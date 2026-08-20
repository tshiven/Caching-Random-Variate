/*
  Name:     rvg_cache.c
  Purpose:  Optional CDF-memoization layer for optimal random variate generation.

  Design:   The public API is a single function, rvg_generate(cdf, prng, dist,
            force_unsafe). Everything about caching -- the hash tables, their
            capacities, per-distribution tuning, and the registry mapping each
            distinct `cdf` seen so far to its own cache -- is private to this
            file. Callers never see or manage a cache object; it is built the
            first time a given `cdf` is used and freed automatically when the
            process exits normally (via atexit), so there is nothing to set
            up or tear down from the outside.

            Two internal-only functions at the bottom (rvg_internal_*) are not
            declared in rvg_cache.h and are not part of the public API -- they
            exist solely so this project's own test/report programs can
            forward-declare them locally and inspect/reset cache state for
            verification purposes.
*/

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "bits.h"
#include "flip.h"
#include "arithmetic.h"
#include "generate.h"

#include "rvg_cache.h"

/* 64 bits (the tree is at max 64 levels deep) */
#define RVG_MAX_LEVELS 64

// Raw-bit float handling

/* Floats are never compared with `==` in this file: a CDF may legitimately
   return -0.0 or a NaN, and those must key distinctly. Everything goes
   through the raw 32-bit pattern, obtained by memcpy to avoid type punning. */

static uint32_t f32_bits(float f) {
    uint32_t u;
    memcpy(&u, &f, sizeof u);
    return u;
}

static float f32_from_bits(uint32_t u) {
    float f;
    memcpy(&f, &u, sizeof f);
    return f;
}

// Cache tuning

#define RVG_HEAD_CAPACITY 65536u  /* head-cache slots, ~1.5 MiB */
#define RVG_TAIL_CAPACITY 16384u  /* tail-cache slots, ~0.6 MiB */
#define RVG_HEAD_MAX_PROBES 6
#define RVG_TAIL_MAX_PROBES 8

/* Max distinct `cdf` functions the registry can manage a cache for at once.
   Fixed, like everything else here -- once full, later new CDFs just fall
   back to an uncached descent instead of erroring out. */
#define RVG_AUTO_REGISTRY_CAPACITY 64u

typedef enum { RVG_STATUS_OK, RVG_STATUS_NOT_RECOMMENDED, RVG_STATUS_UNMEASURED } rvg_status_t;

typedef struct { size_t head_hits, head_misses, tail_hits, tail_misses, head_entries, tail_entries, head_probe_fail; } rvg_stats_t;

// Tables

struct rvg_dist_params {
    unsigned int head_depth;
    int          tail_enabled;
    rvg_status_t status;
};

/* Indexed by rvg_dist_t; order must match the enum in rvg_cache.h. */
static const struct rvg_dist_params RVG_DIST_TABLE[] = {
    /* RVG_DIST_BETA           */ { 28, 0, RVG_STATUS_OK },
    /* RVG_DIST_TDIST          */ { 24, 0, RVG_STATUS_OK },
    /* RVG_DIST_CHISQUARE      */ { 32, 0, RVG_STATUS_OK },
    /* RVG_DIST_FDIST          */ { 24, 0, RVG_STATUS_OK },
    /* RVG_DIST_GAMMA          */ { 28, 0, RVG_STATUS_OK },
    /* RVG_DIST_POISSON        */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_BINOMIAL       */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_NEGBINOMIAL    */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_GEOMETRIC      */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_HYPERGEOMETRIC */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_PASCAL         */ { 32, 1, RVG_STATUS_OK },
    /* RVG_DIST_GAUSSIAN       */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_EXPONENTIAL    */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_LOGNORMAL      */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_CAUCHY         */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_LAPLACE        */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_LOGISTIC       */ { 24, 0, RVG_STATUS_NOT_RECOMMENDED },
    /* RVG_DIST_RAYLEIGH       */ { 24, 0, RVG_STATUS_UNMEASURED },
    /* RVG_DIST_PARETO         */ { 24, 0, RVG_STATUS_UNMEASURED },
    /* RVG_DIST_WEIBULL        */ { 24, 0, RVG_STATUS_UNMEASURED },
    /* RVG_DIST_GUMBEL         */ { 24, 0, RVG_STATUS_UNMEASURED },
};

#define RVG_DIST_COUNT (sizeof(RVG_DIST_TABLE) / sizeof(RVG_DIST_TABLE[0]))

// Cache structures

// Head entry: key (b, l), value = raw bits of the float the CDF returned.
struct head_entry {
    uint64_t b;
    uint32_t l;
    uint32_t val_bits;
    uint8_t  used;
};

// Tail entry: key (b, l, cdf_l, cdf_r, ell), all five compared exactly;
// the two CDF fields are held as raw bit patterns. Value = final double.
struct tail_entry {
    uint64_t b;
    uint32_t l;
    uint32_t cdf_l_bits;
    uint32_t cdf_r_bits;
    uint32_t ell;
    double   value;
    uint8_t  used;
};

struct rvg_cache {
    unsigned int head_depth;
    int          tail_enabled;
    rvg_status_t status;

    struct head_entry *head;
    struct tail_entry *tail;
    size_t head_capacity;   /* fixed at construction, never changes */
    size_t tail_capacity;   /* fixed at construction, never changes */

    rvg_stats_t stats;
};

// Hashing

static uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

static size_t head_hash(uint64_t b, uint32_t l) {
    return (size_t)mix64(mix64(b) ^ ((uint64_t)l * 0x100000001b3ull));
}

static size_t tail_hash(uint64_t b, uint32_t l, uint32_t cl, uint32_t cr, uint32_t ell) {
    uint64_t h = mix64(b);
    h = mix64(h ^ ((uint64_t)l << 32 | (uint64_t)ell));
    h = mix64(h ^ ((uint64_t)cl << 32 | (uint64_t)cr));
    return (size_t)h;
}

// Head cache

/* Return the memoized cdf(d) for (b, l), computing and inserting it on a
   miss. Never resizes and never evicts: after RVG_HEAD_MAX_PROBES occupied,
   non-matching slots the cache gives up and the CDF is called directly. */
static float head_lookup(struct rvg_cache *c, cdf32_t cdf, uint64_t b, uint32_t l, double d) {
    size_t mask = c->head_capacity - 1;
    size_t idx = head_hash(b, l) & mask;

    for (int probe = 0; probe < RVG_HEAD_MAX_PROBES; probe++) {
        struct head_entry *e = &c->head[idx];
        if (!e->used) {
            /* First empty slot in the probe window: miss, then insert here. */
            float v = cdf(d);
            e->b = b;
            e->l = l;
            e->val_bits = f32_bits(v);
            e->used = 1;
            c->stats.head_misses++;
            c->stats.head_entries++;
            return v;
        }
        if (e->b == b && e->l == l) {
            c->stats.head_hits++;
            return f32_from_bits(e->val_bits);
        }
        idx = (idx + 1) & mask;
    }

    /* Probe budget exhausted; fall through to the plain computation. */
    c->stats.head_misses++;
    c->stats.head_probe_fail++;
    return cdf(d);
}


// Tail cache

static int tail_lookup(struct rvg_cache *c, uint64_t b, uint32_t l,
                       uint32_t cl, uint32_t cr, uint32_t ell, double *out) {
    size_t mask = c->tail_capacity - 1;
    size_t idx = tail_hash(b, l, cl, cr, ell) & mask;

    for (int probe = 0; probe < RVG_TAIL_MAX_PROBES; probe++) {
        struct tail_entry *e = &c->tail[idx];
        if (!e->used) { break; }
        if (e->b == b && e->l == l && e->ell == ell
                && e->cdf_l_bits == cl && e->cdf_r_bits == cr) {
            *out = e->value;
            c->stats.tail_hits++;
            return 1;
        }
        idx = (idx + 1) & mask;
    }

    c->stats.tail_misses++;
    return 0;
}

/* Insert if an empty slot is found within the probe budget. An existing entry
   for the same key is left untouched (it holds the same value); a full probe
   window is silently abandoned. Never resizes, never evicts. */
static void tail_insert(struct rvg_cache *c, uint64_t b, uint32_t l,
                        uint32_t cl, uint32_t cr, uint32_t ell, double value) {
    size_t mask = c->tail_capacity - 1;
    size_t idx = tail_hash(b, l, cl, cr, ell) & mask;

    for (int probe = 0; probe < RVG_TAIL_MAX_PROBES; probe++) {
        struct tail_entry *e = &c->tail[idx];
        if (!e->used) {
            e->b = b;
            e->l = l;
            e->cdf_l_bits = cl;
            e->cdf_r_bits = cr;
            e->ell = ell;
            e->value = value;
            e->used = 1;
            c->stats.tail_entries++;
            return;
        }
        if (e->b == b && e->l == l && e->ell == ell
                && e->cdf_l_bits == cl && e->cdf_r_bits == cr) {
            return;
        }
        idx = (idx + 1) & mask;
    }
}

/* Cache lifecycle -- private. Callers never see a cache object.       */

static struct rvg_cache *cache_create(rvg_dist_t dist, int force_unsafe, rvg_status_t *status) {

    if ((unsigned int)dist >= RVG_DIST_COUNT) {
        if (status) { *status = RVG_STATUS_UNMEASURED; }
        return NULL;
    }

    const struct rvg_dist_params *p = &RVG_DIST_TABLE[dist];
    if (status) { *status = p->status; }

    if ((p->status != RVG_STATUS_OK) && !force_unsafe) {
        return NULL;
    }

    struct rvg_cache *c = calloc(1, sizeof *c);
    if (c == NULL) {
        if (status) { *status = RVG_STATUS_UNMEASURED; }
        return NULL;
    }

    c->head_depth   = p->head_depth;
    c->tail_enabled = p->tail_enabled;
    c->status       = p->status;

    c->head_capacity = RVG_HEAD_CAPACITY;
    c->tail_capacity = RVG_TAIL_CAPACITY;

    c->head = calloc(c->head_capacity, sizeof *c->head);
    if (c->head == NULL) {
        free(c);
        if (status) { *status = RVG_STATUS_UNMEASURED; }
        return NULL;
    }

    if (c->tail_enabled) {
        c->tail = calloc(c->tail_capacity, sizeof *c->tail);
        if (c->tail == NULL) {
            free(c->head);
            free(c);
            if (status) { *status = RVG_STATUS_UNMEASURED; }
            return NULL;
        }
    } else {
        c->tail = NULL;
        c->tail_capacity = 0;
    }

    return c;
}

static void cache_free(struct rvg_cache *c) {
    if (c == NULL) { return; }
    free(c->head);
    free(c->tail);
    free(c);
}

/* The descent loop -- a copy of generate_opt's loop (generate.c), so the
   cdf(d) call site can be intercepted; wrapping generate_opt from the
   outside cannot reach that call site. `c` may be NULL for an uncached
   descent, identical to generate_opt. */

struct level_state {
    uint64_t b;
    uint32_t l;
    uint32_t cdf_l_bits;
    uint32_t cdf_r_bits;
    uint32_t ell;
};

static double generate_with_cache(cdf32_t cdf, struct flip_state *prng, struct rvg_cache *c) {

    // Evolving state
    uint64_t b = 0;
    unsigned int ell = 0;
    float cdf_l = 0;
    float cdf_r = 1;

    assert(DBL_SIZE <= RVG_MAX_LEVELS);

    const int tail_on = (c != NULL) && c->tail_enabled && (c->tail != NULL);
    struct level_state levels[RVG_MAX_LEVELS];
    int levels_recorded = 0;
    int last_flip_l = -1;   /* loop index of the most recent flip() call */

    #ifndef NDEBUG
    mpq_t kn0; mpq_init(kn0);
    mpq_t kn1; mpq_init(kn1);
    mpz_t k0; mpz_init(k0);
    mpz_t n0; mpz_init(n0);
    mpz_t k1; mpz_init(k1);
    mpz_t n1; mpz_init(n1);
    #endif

    for (int l = 0; l < DBL_SIZE; l++) {

        uint32_t cdf_l_bits = f32_bits(cdf_l);
        uint32_t cdf_r_bits = f32_bits(cdf_r);

        // Tail cache: the descent from here on consumes no randomness whenever
        // this exact state was previously reached below the last flip().
        if (tail_on) {
            double hit;
            if (tail_lookup(c, b, (uint32_t)l, cdf_l_bits, cdf_r_bits, ell, &hit)) {
                #ifndef NDEBUG
                mpq_clear(kn0); mpq_clear(kn1);
                mpz_clear(k0); mpz_clear(n0);
                mpz_clear(k1); mpz_clear(n1);
                #endif
                return hit;
            }
            levels[levels_recorded].b          = b;
            levels[levels_recorded].l          = (uint32_t)l;
            levels[levels_recorded].cdf_l_bits = cdf_l_bits;
            levels[levels_recorded].cdf_r_bits = cdf_r_bits;
            levels[levels_recorded].ell        = ell;
            levels_recorded++;
        }

        // Compute CDF at midpoint.
        unsigned int m = DBL_SIZE - (l + 1);             // m = n_max - (len(b)+1)
        uint64_t b_lex = (b << (m + 1)) + (1ull << m) - 1; // b+'0' + '1'*m
        uint64_t b_flt = bij64_lex2float(b_lex);
        double d = int2double(b_flt);

        // Head cache: cdf(d) is a pure function of (b, l), memoized shallowly.
        float cdf_m;
        if ((c != NULL) && ((unsigned int)l <= c->head_depth)) {
            cdf_m = head_lookup(c, cdf, b, (uint32_t)l, d);
        } else {
            cdf_m = cdf(d);
        }
        assert(cdf_l <= cdf_m);
        assert(cdf_m <= cdf_r);

        // Compute b+'0' and b+'1'.
        uint64_t b_lex_0 = b << 1;
        uint64_t b_lex_1 = b_lex_0 | 1;

        // Trivial case.
        if (cdf_m == cdf_r) {
            b = b_lex_0;
            cdf_r = cdf_m;
            continue;
        }
        if (cdf_m == cdf_l) {
            b = b_lex_1;
            cdf_l = cdf_m;
            continue;
        }

        // Finite arithmetic case.
        struct subtract_exact_s ss0, ss1;
        subtract_exact(SUB_0, cdf_m, cdf_l, &ss0);
        subtract_exact(SUB_0, cdf_r, cdf_m, &ss1);

        #ifndef NDEBUG
        subtract_gmp(SUB_0, kn0, cdf_m, cdf_l);
        subtract_gmp(SUB_0, kn1, cdf_r, cdf_m);
        mpq_get_num(k0, kn0); mpq_get_den(n0, kn0);
        mpq_get_num(k1, kn1); mpq_get_den(n1, kn1);
        #endif

        if (ell > 0) {
            int a0 = ith_bit_of_exact(&ss0, ell);
            int a1 = ith_bit_of_exact(&ss1, ell);
            #ifndef NDEBUG
            int z0 = ith_bit_of_fraction_gmp(k0, n0, ell);
            int z1 = ith_bit_of_fraction_gmp(k1, n1, ell);
            assert((a0 == z0) && (a1 == z1));
            #endif
            if ((a0 == 1) && (a1 == 0)) {
                b = b_lex_0;
                cdf_r = cdf_m;
                continue;
            }
            if ((a0 == 0) && (a1 == 1)) {
                b = b_lex_1;
                cdf_l = cdf_m;
                continue;
            }
        }

        while (1) {
            ell += 1;
            int a0 = ith_bit_of_exact(&ss0, ell);
            int a1 = ith_bit_of_exact(&ss1, ell);
            #ifndef NDEBUG
            int z0 = ith_bit_of_fraction_gmp(k0, n0, ell);
            int z1 = ith_bit_of_fraction_gmp(k1, n1, ell);
            assert((a0 == z0) && (a1 == z1));
            #endif
            unsigned char x = flip(prng);
            last_flip_l = l;
            if ((x == 0) && (a0 == 1)) {
                b = b_lex_0;
                cdf_r = cdf_m;
                break;
            }
            if ((x == 1) && (a1 == 1)) {
                b = b_lex_1;
                cdf_l = cdf_m;
                break;
            }
        }
    }

    #ifndef NDEBUG
    mpq_clear(kn0); mpq_clear(kn1);
    mpz_clear(k0); mpz_clear(n0);
    mpz_clear(k1); mpz_clear(n1);
    #endif

    b = bij64_lex2float(b);
    double result = int2double(b);

    // Populate the tail cache only now that the sample is complete, and only
    // for levels strictly after the last flip(): from there on the descent was
    // deterministic, so the recorded state alone determines `result`.
    if (tail_on) {
        for (int i = 0; i < levels_recorded; i++) {
            if ((int)levels[i].l <= last_flip_l) { continue; }
            tail_insert(c, levels[i].b, levels[i].l,
                        levels[i].cdf_l_bits, levels[i].cdf_r_bits,
                        levels[i].ell, result);
        }
    }

    return result;
}

/* ------------------------------------------------------------------ */
/* The auto-registry: maps each distinct `cdf` seen so far to its own      */
/* cache (or to NULL, meaning "known -- run this one uncached"). This is   */
/* what makes rvg_generate need no setup/teardown from the caller.         */
/* ------------------------------------------------------------------ */

struct auto_entry {
    cdf32_t cdf;
    struct rvg_cache *cache;
};

static struct auto_entry AUTO_REGISTRY[RVG_AUTO_REGISTRY_CAPACITY];
static size_t AUTO_REGISTRY_COUNT = 0;
static int AUTO_CLEANUP_REGISTERED = 0;

static void auto_registry_cleanup(void) {
    for (size_t i = 0; i < AUTO_REGISTRY_COUNT; i++) {
        cache_free(AUTO_REGISTRY[i].cache);
    }
    AUTO_REGISTRY_COUNT = 0;
}

double rvg_generate(cdf32_t cdf, struct flip_state *prng, rvg_dist_t dist, int force_unsafe) {

    for (size_t i = 0; i < AUTO_REGISTRY_COUNT; i++) {
        if (AUTO_REGISTRY[i].cdf == cdf) {
            return generate_with_cache(cdf, prng, AUTO_REGISTRY[i].cache);
        }
    }

    /* First time this exact cdf has been passed to rvg_generate. */
    struct rvg_cache *c = NULL;
    if (AUTO_REGISTRY_COUNT < RVG_AUTO_REGISTRY_CAPACITY) {
        c = cache_create(dist, force_unsafe, NULL);
        AUTO_REGISTRY[AUTO_REGISTRY_COUNT].cdf = cdf;
        AUTO_REGISTRY[AUTO_REGISTRY_COUNT].cache = c;
        AUTO_REGISTRY_COUNT++;

        if (!AUTO_CLEANUP_REGISTERED) {
            atexit(auto_registry_cleanup);
            AUTO_CLEANUP_REGISTERED = 1;
        }
    }
    /* Registry full: fall through with c == NULL, uncached for this call,
       and not remembered -- every future call for this cdf repeats this
       same linear scan and falls back the same way. No error. */

    return generate_with_cache(cdf, prng, c);
}

// Internal-only helpers (not part of the public API)

rvg_stats_t rvg_internal_stats_for_cdf(cdf32_t cdf) {
    rvg_stats_t empty = {0, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < AUTO_REGISTRY_COUNT; i++) {
        if (AUTO_REGISTRY[i].cdf == cdf) {
            return (AUTO_REGISTRY[i].cache != NULL) ? AUTO_REGISTRY[i].cache->stats : empty;
        }
    }
    return empty;
}

size_t rvg_internal_head_entry_size(void) { return sizeof(struct head_entry); }
size_t rvg_internal_tail_entry_size(void) { return sizeof(struct tail_entry); }

void rvg_internal_reset_for_testing(void) {
    auto_registry_cleanup();
}
