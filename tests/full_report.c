/*
  Name:     full_report.c
  Purpose:  Full accuracy + performance report for rvg_cache across all 11
            supported distributions, each with 3 parameter sets (typical,
            skewed, edge). Writes tests/REPORT.md.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O3 -DNDEBUG -march=native $(gsl-config --cflags) \
                  -I. -Icache -I<gmp include> -L<gmp lib> \
                  tests/full_report.c cache/rvg_cache.c librvg.a \
                  -o tests/full_report.out $(gsl-config --libs-without-cblas) -lgmp
            (run from the repository root).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

#include "generate.h"
#include "flip.h"
#include "rvg_cache.h"

/* ------------------------------------------------------------------ */
/* CDFs: 11 distributions x 3 parameter sets = 33 configs.             */
/* ------------------------------------------------------------------ */

/* --- BETA: typical(5,5), skewed(2,20), edge(0.5,0.5) --- */
MAKE_CDF_P(beta_typical_cdf, gsl_cdf_beta_P, 5, 5)
MAKE_CDF_P(beta_skewed_cdf, gsl_cdf_beta_P, 2, 20)
MAKE_CDF_P(beta_edge_cdf, gsl_cdf_beta_P, 0.5, 0.5)

/* --- TDIST: typical(nu=5), skewed(nu=2, heavy-tailed), edge(nu=1, Cauchy) --- */
MAKE_CDF_P(tdist_typical_cdf, gsl_cdf_tdist_P, 5)
MAKE_CDF_P(tdist_skewed_cdf, gsl_cdf_tdist_P, 2)
MAKE_CDF_P(tdist_edge_cdf, gsl_cdf_tdist_P, 1)

/* --- CHISQUARE: typical(nu=13), skewed(nu=3), edge(nu=0.5) --- */
MAKE_CDF_P(chisq_typical_cdf, gsl_cdf_chisq_P, 13)
MAKE_CDF_P(chisq_skewed_cdf, gsl_cdf_chisq_P, 3)
MAKE_CDF_P(chisq_edge_cdf, gsl_cdf_chisq_P, 0.5)

/* --- FDIST: typical(5,2), skewed(1,10), edge(1,1) --- */
MAKE_CDF_P(fdist_typical_cdf, gsl_cdf_fdist_P, 5, 2)
MAKE_CDF_P(fdist_skewed_cdf, gsl_cdf_fdist_P, 1, 10)
MAKE_CDF_P(fdist_edge_cdf, gsl_cdf_fdist_P, 1, 1)

/* --- GAMMA: typical(shape=0.5,scale=1), skewed(shape=2,scale=3), edge(shape=0.05,scale=1) --- */
MAKE_CDF_P(gamma_typical_cdf, gsl_cdf_gamma_P, 0.5, 1)
MAKE_CDF_P(gamma_skewed_cdf, gsl_cdf_gamma_P, 2, 3)
MAKE_CDF_P(gamma_edge_cdf, gsl_cdf_gamma_P, 0.05, 1)

/* --- POISSON: typical(mu=71), skewed(mu=5), edge(mu=1) --- */
MAKE_CDF_UINT_P(poisson_typical_cdf, gsl_cdf_poisson_P, 71)
MAKE_CDF_UINT_P(poisson_skewed_cdf, gsl_cdf_poisson_P, 5)
MAKE_CDF_UINT_P(poisson_edge_cdf, gsl_cdf_poisson_P, 1)

/* --- BINOMIAL: typical(p=0.2,n=100), skewed(p=0.05,n=200), edge(p=0.5,n=3) --- */
MAKE_CDF_UINT_P(binomial_typical_cdf, gsl_cdf_binomial_P, 0.2, 100)
MAKE_CDF_UINT_P(binomial_skewed_cdf, gsl_cdf_binomial_P, 0.05, 200)
MAKE_CDF_UINT_P(binomial_edge_cdf, gsl_cdf_binomial_P, 0.5, 3)

/* --- NEGBINOMIAL: typical(p=0.71,n=18), skewed(p=0.1,n=5), edge(p=0.9,n=1) --- */
MAKE_CDF_UINT_P(negbinomial_typical_cdf, gsl_cdf_negative_binomial_P, 0.71, 18)
MAKE_CDF_UINT_P(negbinomial_skewed_cdf, gsl_cdf_negative_binomial_P, 0.1, 5)
MAKE_CDF_UINT_P(negbinomial_edge_cdf, gsl_cdf_negative_binomial_P, 0.9, 1)

/* --- GEOMETRIC: typical(p=0.4), skewed(p=0.05), edge(p=0.95) --- */
MAKE_CDF_UINT_P(geometric_typical_cdf, gsl_cdf_geometric_P, 0.4)
MAKE_CDF_UINT_P(geometric_skewed_cdf, gsl_cdf_geometric_P, 0.05)
MAKE_CDF_UINT_P(geometric_edge_cdf, gsl_cdf_geometric_P, 0.95)

/* --- HYPERGEOMETRIC: typical(5,20,7), skewed(2,48,10), edge(1,1,1) --- */
MAKE_CDF_UINT_P(hypergeometric_typical_cdf, gsl_cdf_hypergeometric_P, 5, 20, 7)
MAKE_CDF_UINT_P(hypergeometric_skewed_cdf, gsl_cdf_hypergeometric_P, 2, 48, 10)
MAKE_CDF_UINT_P(hypergeometric_edge_cdf, gsl_cdf_hypergeometric_P, 1, 1, 1)

/* --- PASCAL: typical(p=0.5,n=5), skewed(p=0.1,n=3), edge(p=0.9,n=1) --- */
MAKE_CDF_UINT_P(pascal_typical_cdf, gsl_cdf_pascal_P, 0.5, 5)
MAKE_CDF_UINT_P(pascal_skewed_cdf, gsl_cdf_pascal_P, 0.1, 3)
MAKE_CDF_UINT_P(pascal_edge_cdf, gsl_cdf_pascal_P, 0.9, 1)

/* ------------------------------------------------------------------ */
/* Config tables                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *dist_name;
    const char *param_label;
    const char *variant;    /* "typical", "skewed", "edge" */
    rvg_dist_t  dist;
    cdf32_t     cdf;
} config_t;

static config_t CONFIGS[] = {
    {"BETA",           "Beta(5,5)",         "typical", RVG_DIST_BETA,           beta_typical_cdf},
    {"BETA",           "Beta(2,20)",        "skewed",  RVG_DIST_BETA,           beta_skewed_cdf},
    {"BETA",           "Beta(0.5,0.5)",     "edge",    RVG_DIST_BETA,           beta_edge_cdf},

    {"TDIST",          "Tdist(5)",          "typical", RVG_DIST_TDIST,          tdist_typical_cdf},
    {"TDIST",          "Tdist(2)",          "skewed",  RVG_DIST_TDIST,          tdist_skewed_cdf},
    {"TDIST",          "Tdist(1)",          "edge",    RVG_DIST_TDIST,          tdist_edge_cdf},

    {"CHISQUARE",      "ChiSquare(13)",     "typical", RVG_DIST_CHISQUARE,      chisq_typical_cdf},
    {"CHISQUARE",      "ChiSquare(3)",      "skewed",  RVG_DIST_CHISQUARE,      chisq_skewed_cdf},
    {"CHISQUARE",      "ChiSquare(0.5)",    "edge",    RVG_DIST_CHISQUARE,      chisq_edge_cdf},

    {"FDIST",          "Fdist(5,2)",        "typical", RVG_DIST_FDIST,          fdist_typical_cdf},
    {"FDIST",          "Fdist(1,10)",       "skewed",  RVG_DIST_FDIST,          fdist_skewed_cdf},
    {"FDIST",          "Fdist(1,1)",        "edge",    RVG_DIST_FDIST,          fdist_edge_cdf},

    {"GAMMA",          "Gamma(0.5,1)",      "typical", RVG_DIST_GAMMA,          gamma_typical_cdf},
    {"GAMMA",          "Gamma(2,3)",        "skewed",  RVG_DIST_GAMMA,          gamma_skewed_cdf},
    {"GAMMA",          "Gamma(0.05,1)",     "edge",    RVG_DIST_GAMMA,          gamma_edge_cdf},

    {"POISSON",        "Poisson(71)",       "typical", RVG_DIST_POISSON,        poisson_typical_cdf},
    {"POISSON",        "Poisson(5)",        "skewed",  RVG_DIST_POISSON,        poisson_skewed_cdf},
    {"POISSON",        "Poisson(1)",        "edge",    RVG_DIST_POISSON,        poisson_edge_cdf},

    {"BINOMIAL",       "Binomial(0.2,100)", "typical", RVG_DIST_BINOMIAL,       binomial_typical_cdf},
    {"BINOMIAL",       "Binomial(0.05,200)","skewed",  RVG_DIST_BINOMIAL,       binomial_skewed_cdf},
    {"BINOMIAL",       "Binomial(0.5,3)",   "edge",    RVG_DIST_BINOMIAL,       binomial_edge_cdf},

    {"NEGBINOMIAL",    "NegBinomial(0.71,18)","typical", RVG_DIST_NEGBINOMIAL,  negbinomial_typical_cdf},
    {"NEGBINOMIAL",    "NegBinomial(0.1,5)",  "skewed",  RVG_DIST_NEGBINOMIAL,  negbinomial_skewed_cdf},
    {"NEGBINOMIAL",    "NegBinomial(0.9,1)",  "edge",    RVG_DIST_NEGBINOMIAL,  negbinomial_edge_cdf},

    {"GEOMETRIC",      "Geometric(0.4)",    "typical", RVG_DIST_GEOMETRIC,      geometric_typical_cdf},
    {"GEOMETRIC",      "Geometric(0.05)",   "skewed",  RVG_DIST_GEOMETRIC,      geometric_skewed_cdf},
    {"GEOMETRIC",      "Geometric(0.95)",   "edge",    RVG_DIST_GEOMETRIC,      geometric_edge_cdf},

    {"HYPERGEOMETRIC", "Hypergeometric(5,20,7)", "typical", RVG_DIST_HYPERGEOMETRIC, hypergeometric_typical_cdf},
    {"HYPERGEOMETRIC", "Hypergeometric(2,48,10)","skewed",  RVG_DIST_HYPERGEOMETRIC, hypergeometric_skewed_cdf},
    {"HYPERGEOMETRIC", "Hypergeometric(1,1,1)",  "edge",    RVG_DIST_HYPERGEOMETRIC, hypergeometric_edge_cdf},

    {"PASCAL",         "Pascal(0.5,5)",     "typical", RVG_DIST_PASCAL,         pascal_typical_cdf},
    {"PASCAL",         "Pascal(0.1,3)",     "skewed",  RVG_DIST_PASCAL,         pascal_skewed_cdf},
    {"PASCAL",         "Pascal(0.9,1)",     "edge",    RVG_DIST_PASCAL,         pascal_edge_cdf},
};
#define NUM_CONFIGS ((int)(sizeof(CONFIGS) / sizeof(CONFIGS[0])))
#define SEEDS_PER_CONFIG 20000

typedef struct {
    const char *name;
    rvg_dist_t  dist;
    cdf32_t     cdf;
} base_dist_t;

/* One entry per distribution, using its typical-parameter CDF, for Section 2. */
static base_dist_t BASE_DISTS[] = {
    {"BETA",           RVG_DIST_BETA,           beta_typical_cdf},
    {"TDIST",          RVG_DIST_TDIST,          tdist_typical_cdf},
    {"CHISQUARE",      RVG_DIST_CHISQUARE,      chisq_typical_cdf},
    {"FDIST",          RVG_DIST_FDIST,          fdist_typical_cdf},
    {"GAMMA",          RVG_DIST_GAMMA,          gamma_typical_cdf},
    {"POISSON",        RVG_DIST_POISSON,        poisson_typical_cdf},
    {"BINOMIAL",       RVG_DIST_BINOMIAL,       binomial_typical_cdf},
    {"NEGBINOMIAL",    RVG_DIST_NEGBINOMIAL,    negbinomial_typical_cdf},
    {"GEOMETRIC",      RVG_DIST_GEOMETRIC,      geometric_typical_cdf},
    {"HYPERGEOMETRIC", RVG_DIST_HYPERGEOMETRIC, hypergeometric_typical_cdf},
    {"PASCAL",         RVG_DIST_PASCAL,         pascal_typical_cdf},
};
#define NUM_BASE_DISTS ((int)(sizeof(BASE_DISTS) / sizeof(BASE_DISTS[0])))

static const long N_LIST[] = {100, 250, 500, 1000, 2000, 5000, 10000, 25000, 50000, 100000};
#define NUM_N ((int)(sizeof(N_LIST) / sizeof(N_LIST[0])))

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

static int double_bits_eq(double a, double b) {
    uint64_t ba, bb;
    memcpy(&ba, &a, sizeof ba);
    memcpy(&bb, &b, sizeof bb);
    return ba == bb;
}

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ------------------------------------------------------------------ */
/* Section 1: accuracy                                                  */
/* ------------------------------------------------------------------ */

static long section1_mismatches[NUM_CONFIGS];
static int section1_ran[NUM_CONFIGS];

static long run_accuracy(config_t *cfg) {
    rvg_status_t status;
    rvg_cache_t *cache = rvg_cache_create(cfg->dist, 1, &status);
    if (cache == NULL) {
        printf("  %-16s %-24s : rvg_cache_create returned NULL even with "
               "force_unsafe=1 -- cannot test.\n", cfg->dist_name, cfg->param_label);
        return -1;
    }

    long mismatches = 0;
    for (unsigned long seed = 1; seed <= SEEDS_PER_CONFIG; seed++) {

        gsl_rng *r1 = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(r1, seed);
        struct flip_state p1 = make_flip_state(r1);
        double v1 = generate_opt(cfg->cdf, &p1);
        unsigned long f1 = p1.num_flips;
        gsl_rng_free(r1);

        gsl_rng *r2 = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(r2, seed);
        struct flip_state p2 = make_flip_state(r2);
        double v2 = rvg_generate(cfg->cdf, &p2, cache);
        unsigned long f2 = p2.num_flips;
        gsl_rng_free(r2);

        if (!double_bits_eq(v1, v2) || (f1 != f2)) {
            mismatches++;
            printf("  MISMATCH seed=%lu dist=%s params=%s baseline=%.17g (flips=%lu) "
                   "cached=%.17g (flips=%lu)\n",
                   seed, cfg->dist_name, cfg->param_label, v1, f1, v2, f2);
        }
    }

    rvg_cache_free(cache);
    return mismatches;
}

/* ------------------------------------------------------------------ */
/* Section 2: performance (typical params only)                        */
/* ------------------------------------------------------------------ */

typedef struct {
    double      baseline_sec;
    double      cached_sec;
    double      speedup;
    rvg_stats_t stats;
} perf_point_t;

/* Every (distribution, N) timing is repeated across these REPEATS fixed,
   distinct seeds -- a single wall-clock sample at N=100..2000 is only a few
   milliseconds and easily dominated by scheduling noise, so one number alone
   is not trustworthy. Each repeat uses a brand-new, no-warmup cache. */
#define REPEATS 5
static const unsigned long PERF_SEEDS[REPEATS] = {42, 1337, 271828, 986532, 555555};

static perf_point_t run_perf(cdf32_t cdf, rvg_dist_t dist, long n, unsigned long seed) {
    perf_point_t res;

    gsl_rng *r1 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r1, seed);
    struct flip_state p1 = make_flip_state(r1);
    double t0 = now_seconds();
    for (long i = 0; i < n; i++) { generate_opt(cdf, &p1); }
    double t1 = now_seconds();
    res.baseline_sec = t1 - t0;
    gsl_rng_free(r1);

    rvg_status_t status;
    rvg_cache_t *cache = rvg_cache_create(dist, 1, &status);
    gsl_rng *r2 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r2, seed);
    struct flip_state p2 = make_flip_state(r2);
    double t2 = now_seconds();
    for (long i = 0; i < n; i++) { rvg_generate(cdf, &p2, cache); }
    double t3 = now_seconds();
    res.cached_sec = t3 - t2;
    res.stats = rvg_cache_stats(cache);
    gsl_rng_free(r2);
    rvg_cache_free(cache);

    res.speedup = (res.cached_sec > 0.0) ? (res.baseline_sec / res.cached_sec) : -1.0;
    return res;
}

/* Ascending insertion sort of a small fixed-size double array; REPEATS is
   small enough that this is simpler than pulling in qsort with a comparator. */
static void sort5(double *v) {
    for (int i = 1; i < REPEATS; i++) {
        double key = v[i];
        int j = i - 1;
        while (j >= 0 && v[j] > key) { v[j + 1] = v[j]; j--; }
        v[j + 1] = key;
    }
}

static double median5(const double *raw) {
    double tmp[REPEATS];
    memcpy(tmp, raw, sizeof tmp);
    sort5(tmp);
    return tmp[REPEATS / 2];
}

typedef struct {
    double speedup_median, speedup_min, speedup_max;
    double baseline_median_sec, cached_median_sec;
    double head_hits_median, head_misses_median, head_probe_fail_median, head_entries_median;
    double tail_hits_median, tail_misses_median, tail_entries_median;
} perf_summary_t;

static perf_summary_t summarize_perf(const perf_point_t *raw /* [REPEATS] */) {
    perf_summary_t s;
    double speedups[REPEATS], baselines[REPEATS], cacheds[REPEATS];
    double hh[REPEATS], hm[REPEATS], hpf[REPEATS], he[REPEATS];
    double th[REPEATS], tm[REPEATS], te[REPEATS];

    s.speedup_min = raw[0].speedup;
    s.speedup_max = raw[0].speedup;
    for (int k = 0; k < REPEATS; k++) {
        speedups[k] = raw[k].speedup;
        baselines[k] = raw[k].baseline_sec;
        cacheds[k] = raw[k].cached_sec;
        hh[k] = (double)raw[k].stats.head_hits;
        hm[k] = (double)raw[k].stats.head_misses;
        hpf[k] = (double)raw[k].stats.head_probe_fail;
        he[k] = (double)raw[k].stats.head_entries;
        th[k] = (double)raw[k].stats.tail_hits;
        tm[k] = (double)raw[k].stats.tail_misses;
        te[k] = (double)raw[k].stats.tail_entries;
        if (raw[k].speedup < s.speedup_min) { s.speedup_min = raw[k].speedup; }
        if (raw[k].speedup > s.speedup_max) { s.speedup_max = raw[k].speedup; }
    }

    s.speedup_median = median5(speedups);
    s.baseline_median_sec = median5(baselines);
    s.cached_median_sec = median5(cacheds);
    s.head_hits_median = median5(hh);
    s.head_misses_median = median5(hm);
    s.head_probe_fail_median = median5(hpf);
    s.head_entries_median = median5(he);
    s.tail_hits_median = median5(th);
    s.tail_misses_median = median5(tm);
    s.tail_entries_median = median5(te);
    return s;
}

static perf_point_t section2_raw[NUM_BASE_DISTS][NUM_N][REPEATS];
static perf_summary_t section2_results[NUM_BASE_DISTS][NUM_N];

/* ------------------------------------------------------------------ */
/* Main                                                                  */
/* ------------------------------------------------------------------ */

int main(void) {

    printf("=== SECTION 1: ACCURACY (33 configs x %d seeds, cold-start) ===\n\n", SEEDS_PER_CONFIG);

    long total_mismatches = 0;
    for (int i = 0; i < NUM_CONFIGS; i++) {
        config_t *cfg = &CONFIGS[i];
        long m = run_accuracy(cfg);
        section1_mismatches[i] = m;
        section1_ran[i] = (m >= 0);
        if (m >= 0) {
            printf("%-16s %-24s (%-7s): %ld/%d mismatches\n",
                   cfg->dist_name, cfg->param_label, cfg->variant, m, SEEDS_PER_CONFIG);
            total_mismatches += m;
        } else {
            total_mismatches += 1; /* count the create-failure itself as a hard failure */
        }
    }

    printf("\nSection 1 total mismatches across all configs: %ld\n", total_mismatches);

    if (total_mismatches != 0) {
        printf("\n*** SECTION 1 FAILED -- stopping before Section 2 and before writing REPORT.md. ***\n");
        printf("*** See MISMATCH lines above for exact seed/distribution/parameter/value/flip detail. ***\n");
        return 1;
    }

    printf("\n=== SECTION 2: PERFORMANCE (typical params only, %d distributions x %d N-values x %d repeats) ===\n\n",
           NUM_BASE_DISTS, NUM_N, REPEATS);

    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        base_dist_t *bd = &BASE_DISTS[i];
        for (int j = 0; j < NUM_N; j++) {
            perf_point_t raw[REPEATS];
            for (int k = 0; k < REPEATS; k++) {
                raw[k] = run_perf(bd->cdf, bd->dist, N_LIST[j], PERF_SEEDS[k]);
                section2_raw[i][j][k] = raw[k];
            }
            perf_summary_t summ = summarize_perf(raw);
            section2_results[i][j] = summ;
            printf("%-16s N=%-7ld speedup median=%.3fx range=[%.3fx, %.3fx] "
                   "(baseline med=%.6fs cached med=%.6fs)\n",
                   bd->name, N_LIST[j], summ.speedup_median, summ.speedup_min, summ.speedup_max,
                   summ.baseline_median_sec, summ.cached_median_sec);
        }
        printf("\n");
    }

    /* Breakeven N per distribution: first N such that the median speedup is
       >= 1.0 for that N and every larger N tested. */
    int breakeven_idx[NUM_BASE_DISTS];
    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        breakeven_idx[i] = -1;
        for (int start = 0; start < NUM_N; start++) {
            int ok = 1;
            for (int j = start; j < NUM_N; j++) {
                if (section2_results[i][j].speedup_median < 1.0) { ok = 0; break; }
            }
            if (ok) { breakeven_idx[i] = start; break; }
        }
    }

    printf("=== Breakeven N ===\n\n");
    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        if (breakeven_idx[i] >= 0) {
            printf("%-16s breakeven N = %ld\n", BASE_DISTS[i].name, N_LIST[breakeven_idx[i]]);
        } else {
            printf("%-16s breakeven N = N/A - never breaks even\n", BASE_DISTS[i].name);
        }
    }

    /* ------------------------------------------------------------------ */
    /* SECTION 3: write tests/REPORT.md                                    */
    /* ------------------------------------------------------------------ */

    printf("\n=== SECTION 3: writing tests/REPORT.md ===\n\n");

    FILE *f = fopen("tests/REPORT.md", "w");
    if (f == NULL) {
        fprintf(stderr, "Could not open tests/REPORT.md for writing.\n");
        return 1;
    }

    fprintf(f, "# rvg_cache Full Report\n\n");

    /* --- Methodology --- */
    fprintf(f, "## Methodology\n\n");
    fprintf(f,
        "**Accuracy (Section 1).** 11 distributions x 3 parameter sets each "
        "(typical, skewed, edge) = 33 configs. For each config, one `rvg_cache_t` "
        "is created with `rvg_cache_create(dist, force_unsafe=1, &status)` and "
        "then reused, unmodified, across %d consecutive seeds (1..%d) -- there is "
        "no separate warmup pass, so the cache is cold on seed 1 and matures live "
        "as the loop proceeds, which is how it would be used in practice. For each "
        "seed, `generate_opt` is run against a fresh `gsl_rng` seeded with that "
        "seed, and `rvg_generate` (same cache, same `cdf`) is run against a "
        "*separate* fresh `gsl_rng` seeded identically. The two results are "
        "compared by `memcpy`-ing each `double` to a `uint64_t` and comparing raw "
        "bits (never `==`, so a NaN or signed-zero divergence would still be "
        "caught), and the two `flip_state.num_flips` counters are compared "
        "exactly, which additionally confirms the cached path consumed the exact "
        "same amount of randomness as the uncached path. `force_unsafe=1` is "
        "passed uniformly since this run is deliberately stress-testing skewed "
        "and edge parameterizations, not just recommended ones; all 11 "
        "distributions in this report are tabulated as `RVG_STATUS_OK` "
        "regardless, so this does not change which entries get built.\n\n",
        SEEDS_PER_CONFIG, SEEDS_PER_CONFIG);
    fprintf(f,
        "**Performance (Section 2).** Typical parameter sets only (11 "
        "distributions, not 33 configs). For each of N in {%s}, the timing is "
        "repeated %d times with %d different fixed seeds (%s), and *each* "
        "repeat uses its own fresh, no-warmup cache and its own fresh `gsl_rng` "
        "-- `generate_opt` and `rvg_generate` in a given repeat see the same "
        "seed as each other, but each repeat uses a different seed from the "
        "others. Wall-clock time is `clock_gettime(CLOCK_MONOTONIC)` around "
        "each N-sample loop, and `speedup = baseline_seconds / cached_seconds` "
        "per repeat. What is reported at each N is the **median** speedup "
        "across the %d repeats plus the **[min, max] range**, not a single "
        "number: single-shot wall-clock timing at N=100..2000 is only a few "
        "milliseconds and is easily dominated by scheduler noise, so a lone "
        "sample would misrepresent how stable the effect actually is. "
        "**Breakeven N** is the smallest tested N at which the *median* "
        "speedup is >= 1.0 *and* stays >= 1.0 for every larger N tested (i.e. "
        "the point past which the cache is never seen to fall behind again in "
        "this run, on the median); \"N/A - never breaks even\" means no such N "
        "was found among the values tested. **Memory@100000** is the *median* "
        "(across the %d repeats) of `rvg_cache_stats()` taken on each N=100000 "
        "cache before it is freed: `head_entries * rvg_cache_head_entry_size() "
        "+ tail_entries * rvg_cache_tail_entry_size()`, i.e. the actual "
        "occupied-slot count times the real struct size on this build, not a "
        "hardcoded guess.\n\n",
        "100, 250, 500, 1000, 2000, 5000, 10000, 25000, 50000, 100000",
        REPEATS, REPEATS, "42, 1337, 271828, 986532, 555555", REPEATS, REPEATS);

    fprintf(f,
        "**Head-cache saturation (below).** BETA, CHISQUARE, and GAMMA are all "
        "head-cache-only distributions (`tail_enabled=0` in the internal table) "
        "with a %u-slot fixed-capacity head table (`RVG_HEAD_CAPACITY`). For "
        "these three, the table below Section 2 additionally reports the median "
        "(over the same %d repeats) `head_entries`, `head_hits`, `head_misses`, "
        "and `head_probe_fail` (the subset of misses where the probe window was "
        "full of *other* keys and the cache gave up without inserting, i.e. "
        "wasted probe work with no caching benefit) at each N, to show what "
        "actually happens as the table fills toward capacity.\n\n",
        (unsigned)RVG_HEAD_CAPACITY, REPEATS);

    /* --- Accuracy table --- */
    fprintf(f, "## Accuracy\n\n");
    fprintf(f, "| Distribution | Parameters | Variant | Seeds Tested | Mismatches |\n");
    fprintf(f, "|---|---|---|---:|---:|\n");
    for (int i = 0; i < NUM_CONFIGS; i++) {
        config_t *cfg = &CONFIGS[i];
        fprintf(f, "| %s | %s | %s | %d | %ld |\n",
                cfg->dist_name, cfg->param_label, cfg->variant,
                SEEDS_PER_CONFIG, section1_mismatches[i]);
    }
    fprintf(f, "\n");

    /* --- Performance table --- */
    fprintf(f, "## Performance (typical parameters, speedup = baseline/cached, %d repeats: median [min, max])\n\n", REPEATS);
    fprintf(f, "| Distribution ");
    for (int j = 0; j < NUM_N; j++) { fprintf(f, "| N=%ld ", N_LIST[j]); }
    fprintf(f, "|\n|---");
    for (int j = 0; j < NUM_N; j++) { fprintf(f, "|---:"); }
    fprintf(f, "|\n");
    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        fprintf(f, "| %s ", BASE_DISTS[i].name);
        for (int j = 0; j < NUM_N; j++) {
            perf_summary_t *s = &section2_results[i][j];
            fprintf(f, "| %.3fx [%.3f, %.3f] ", s->speedup_median, s->speedup_min, s->speedup_max);
        }
        fprintf(f, "|\n");
    }
    fprintf(f, "\n");

    /* --- Head-cache saturation table (BETA, CHISQUARE, GAMMA) --- */
    fprintf(f, "## Head-Cache Saturation: BETA, CHISQUARE, GAMMA\n\n");
    fprintf(f, "Head-cache-only distributions (no tail cache); `RVG_HEAD_CAPACITY` = %u slots. "
               "All counts are the median over the %d repeats at that N.\n\n",
            (unsigned)RVG_HEAD_CAPACITY, REPEATS);

    const char *sat_names[] = {"BETA", "CHISQUARE", "GAMMA"};
    long sat_saturation_n[3];       /* -1 if never saturated within tested N */
    double sat_speedup_at_sat[3];   /* median speedup at the saturation N (or N/A) */
    double sat_speedup_at_end[3];   /* median speedup at N=100000 */

    for (int s = 0; s < 3; s++) {
        int di = -1;
        for (int i = 0; i < NUM_BASE_DISTS; i++) {
            if (strcmp(BASE_DISTS[i].name, sat_names[s]) == 0) { di = i; break; }
        }
        sat_saturation_n[s] = -1;
        sat_speedup_at_sat[s] = -1.0;
        sat_speedup_at_end[s] = -1.0;
        if (di < 0) { continue; }

        fprintf(f, "### %s\n\n", sat_names[s]);
        fprintf(f, "| N | head_entries | head_hits | head_misses | head_probe_fail | hit rate | speedup (median) |\n");
        fprintf(f, "|---:|---:|---:|---:|---:|---:|---:|\n");

        for (int j = 0; j < NUM_N; j++) {
            perf_summary_t *sm = &section2_results[di][j];
            double total_lookups = sm->head_hits_median + sm->head_misses_median;
            double hit_rate = (total_lookups > 0.0) ? (sm->head_hits_median / total_lookups) : 0.0;
            fprintf(f, "| %ld | %.0f | %.0f | %.0f | %.0f | %.1f%% | %.3fx |\n",
                    N_LIST[j], sm->head_entries_median, sm->head_hits_median,
                    sm->head_misses_median, sm->head_probe_fail_median,
                    hit_rate * 100.0, sm->speedup_median);
            if (sat_saturation_n[s] < 0 && sm->head_entries_median >= (double)RVG_HEAD_CAPACITY) {
                sat_saturation_n[s] = N_LIST[j];
                sat_speedup_at_sat[s] = sm->speedup_median;
            }
            if (j == NUM_N - 1) { sat_speedup_at_end[s] = sm->speedup_median; }
        }
        fprintf(f, "\n");

        if (sat_saturation_n[s] >= 0) {
            fprintf(f, "**Observed saturation:** the head table for %s first reaches "
                       "its %u-entry capacity at N=%ld among the tested N values "
                       "(median across repeats).\n\n", sat_names[s], (unsigned)RVG_HEAD_CAPACITY, sat_saturation_n[s]);
        } else {
            fprintf(f, "**Observed saturation:** the head table for %s does not reach "
                       "its %u-entry capacity within the tested N range (up to N=%ld).\n\n",
                       sat_names[s], (unsigned)RVG_HEAD_CAPACITY, N_LIST[NUM_N - 1]);
        }
    }

    /* --- Summary table --- */
    fprintf(f, "## Summary\n\n");
    fprintf(f, "| Distribution | Breakeven N | Speedup@1000 (median [min,max]) | Speedup@10000 (median [min,max]) | Speedup@100000 (median [min,max]) | Memory@100000 (median) | Total Mismatches (3 param sets) |\n");
    fprintf(f, "|---|---:|---:|---:|---:|---:|---:|\n");

    /* indices into N_LIST for 1000, 10000, 100000 */
    int idx_1000 = -1, idx_10000 = -1, idx_100000 = -1;
    for (int j = 0; j < NUM_N; j++) {
        if (N_LIST[j] == 1000) { idx_1000 = j; }
        if (N_LIST[j] == 10000) { idx_10000 = j; }
        if (N_LIST[j] == 100000) { idx_100000 = j; }
    }

    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        long total_for_dist = 0;
        for (int k = 0; k < NUM_CONFIGS; k++) {
            if (strcmp(CONFIGS[k].dist_name, BASE_DISTS[i].name) == 0) {
                total_for_dist += section1_mismatches[k];
            }
        }

        perf_summary_t *s1000 = &section2_results[i][idx_1000];
        perf_summary_t *s10000 = &section2_results[i][idx_10000];
        perf_summary_t *s100000 = &section2_results[i][idx_100000];

        double mem_bytes = s100000->head_entries_median * (double)rvg_cache_head_entry_size()
                          + s100000->tail_entries_median * (double)rvg_cache_tail_entry_size();

        char breakeven_str[64];
        if (breakeven_idx[i] >= 0) {
            snprintf(breakeven_str, sizeof breakeven_str, "%ld", N_LIST[breakeven_idx[i]]);
        } else {
            snprintf(breakeven_str, sizeof breakeven_str, "N/A - never breaks even");
        }

        fprintf(f, "| %s | %s | %.3fx [%.3f,%.3f] | %.3fx [%.3f,%.3f] | %.3fx [%.3f,%.3f] | %.1f KB | %ld |\n",
                BASE_DISTS[i].name, breakeven_str,
                s1000->speedup_median, s1000->speedup_min, s1000->speedup_max,
                s10000->speedup_median, s10000->speedup_min, s10000->speedup_max,
                s100000->speedup_median, s100000->speedup_min, s100000->speedup_max,
                mem_bytes / 1024.0,
                total_for_dist);
    }
    fprintf(f, "\n");

    /* --- Limitations --- */
    fprintf(f, "## Limitations\n\n");
    fprintf(f,
        "1. A cache's keys carry no distribution-parameter or CDF identity beyond "
        "the `cdf` function pointer it binds to on first use (see `rvg_cache.h`); "
        "in practice this means **the distribution parameters must be fixed for "
        "the entire lifetime of one cache object**. Changing parameters requires "
        "creating a new cache -- this report creates one cache per config for "
        "exactly that reason.\n");
    fprintf(f,
        "2. Results for GAUSSIAN, EXPONENTIAL, LOGNORMAL, CAUCHY, LAPLACE, and "
        "LOGISTIC are **not included in this report**: they are tabulated as "
        "`RVG_STATUS_NOT_RECOMMENDED`, meaning caching was measured for them and "
        "was not shown to provide a reliable speedup. RAYLEIGH, PARETO, WEIBULL, "
        "and GUMBEL are tabulated as `RVG_STATUS_UNMEASURED` and are only "
        "correctness-tested separately (not for speed here); nothing in this "
        "report should be read as a claim about any of these ten distributions.\n");
    {
        FILE *unamep = popen("uname -a", "r");
        char unamebuf[512] = "unknown";
        if (unamep) {
            if (fgets(unamebuf, sizeof unamebuf, unamep) == NULL) { unamebuf[0] = '\0'; }
            pclose(unamep);
        }
        size_t ulen = strlen(unamebuf);
        if (ulen > 0 && unamebuf[ulen - 1] == '\n') { unamebuf[ulen - 1] = '\0'; }

        FILE *ccp = popen("gcc --version 2>&1 | head -1", "r");
        char ccbuf[256] = "unknown";
        if (ccp) {
            if (fgets(ccbuf, sizeof ccbuf, ccp) == NULL) { ccbuf[0] = '\0'; }
            pclose(ccp);
        }
        size_t clen = strlen(ccbuf);
        if (clen > 0 && ccbuf[clen - 1] == '\n') { ccbuf[clen - 1] = '\0'; }

        fprintf(f,
            "3. This report was generated on a single machine "
            "(`uname -a`: `%s`; `gcc --version`: `%s` -- on macOS this `gcc` is "
            "Apple Clang, not upstream GCC), compiled with "
            "`-O3 -DNDEBUG -march=native`, and has **not been cross-validated on "
            "other hardware, OS, or compiler combinations**. Absolute and "
            "relative timings (especially `-march=native`, which tunes for this "
            "specific CPU) should not be assumed to transfer to other machines. "
            "**Compiler used for `librvg.a` vs. the cache/test code:** both were "
            "built by the exact same `gcc` invocation on this machine -- the root "
            "Makefile's object rule (`%%.o: %%.c %%.h`) calls the literal command "
            "`gcc`, which on this machine resolves via `/usr/bin/gcc` to Apple "
            "Clang %s (confirmed: `which gcc` -> `/usr/bin/gcc`, and no "
            "unversioned `gcc` exists under `/opt/homebrew/bin`, only "
            "`gcc-16`); the manual build command used for this report's test "
            "binaries also invokes plain `gcc`, so it resolves to the identical "
            "compiler. A separate Homebrew GCC (`gcc-16`, real GNU GCC) is "
            "present on this machine and *is* required elsewhere in this "
            "project -- specifically for `examples/main.c`, which defines its "
            "CDFs with `MAKE_CDF_P`/`MAKE_CDF_UINT_P` *inside* `main()`, a GNU C "
            "nested-function definition that Apple Clang rejects outright "
            "(confirmed directly: a minimal nested-function test file fails "
            "under `/usr/bin/gcc` with \"function definition is not allowed "
            "here\" and compiles and runs correctly under `gcc-16`). Neither "
            "`generate.c`/the rest of core `librvg.a`, nor `cache/rvg_cache.c`, "
            "nor `tests/full_report.c` use nested functions -- this report's CDF "
            "macros are deliberately invoked at file scope, not inside `main()`, "
            "specifically so they do not need that extension -- so this project's "
            "Homebrew-GCC requirement does not apply to anything measured here, "
            "and both the baseline and cached code paths in this report were "
            "compiled by the same compiler binary. This also means fairness is "
            "not a cross-compiler question at all in this run: `rvg_generate` "
            "does not reimplement `generate_opt`'s inner arithmetic "
            "(`subtract_exact`, `ith_bit_of_exact`, `bij64_lex2float`, "
            "`int2double`) -- it links directly against the same compiled "
            "definitions from `librvg.a` that `generate_opt` itself calls -- so "
            "even under a hypothetical mixed-compiler build, both paths would "
            "still be exercising identical machine code for that shared "
            "arithmetic; only the outer descent loop (copied into "
            "`rvg_cache.c` because the `cdf(d)` call site had to be "
            "intercepted) is separately compiled, and it is compiled once, by "
            "one compiler, in any given build.\n",
            unamebuf, ccbuf, ccbuf);
    }

    {
        char sat_str[3][32];
        char cmp_str[3][96]; /* "1.396x -> 1.372x (degrading)" or "never saturated in this run" */
        const double TREND_EPS = 0.02; /* 0.02x = below wall-clock repeat noise */

        for (int s = 0; s < 3; s++) {
            if (sat_saturation_n[s] >= 0) {
                snprintf(sat_str[s], sizeof sat_str[s], "N=%ld", sat_saturation_n[s]);

                double delta = sat_speedup_at_end[s] - sat_speedup_at_sat[s];
                const char *trend = (delta > TREND_EPS) ? "still climbing"
                                  : (delta < -TREND_EPS) ? "degrading"
                                  : "flat";
                snprintf(cmp_str[s], sizeof cmp_str[s], "%.3fx -> %.3fx (**%s**)",
                         sat_speedup_at_sat[s], sat_speedup_at_end[s], trend);
            } else {
                snprintf(sat_str[s], sizeof sat_str[s], "not reached by N=%ld", N_LIST[NUM_N - 1]);
                snprintf(cmp_str[s], sizeof cmp_str[s], "n/a (never saturated in this run)");
            }
        }

        fprintf(f,
            "4. **Head-cache saturation (BETA, CHISQUARE, GAMMA).** All three are "
            "head-cache-only (no tail cache) and share the %u-slot fixed "
            "`RVG_HEAD_CAPACITY`. Observed in this run (median `head_entries` "
            "first reaching capacity; see the per-N tables above): BETA "
            "saturates at %s, CHISQUARE at %s, GAMMA at %s. Past that point, new "
            "(b, l) keys can no longer be inserted -- lookups for keys not "
            "already present either still hit (if that exact key was cached "
            "before saturation) or become `head_probe_fail` (6 wasted probes, "
            "then a direct, uncached `cdf()` call). Comparing median speedup at "
            "the saturation N against median speedup at N=100000 (\"before -> "
            "after\", classified as climbing/flat/degrading using a +/-%.2fx "
            "dead zone around zero change, since a difference smaller than that "
            "is within the noise already visible in this section's [min, max] "
            "repeat ranges): BETA %s; CHISQUARE %s; GAMMA %s. This should not be "
            "read as a guarantee for other distributions, other "
            "parameterizations, or N well beyond 100000 -- it is only what this "
            "run observed, and a distribution whose recurring (b, l) prefixes "
            "are less concentrated at shallow levels could plausibly behave "
            "differently once its table saturates.\n",
            (unsigned)RVG_HEAD_CAPACITY,
            sat_str[0], sat_str[1], sat_str[2],
            TREND_EPS,
            cmp_str[0], cmp_str[1], cmp_str[2]);
    }

    fclose(f);
    printf("Wrote tests/REPORT.md\n");

    printf("\nALL SECTION 1 CHECKS PASSED (0 mismatches across all 33 configs)\n");
    return 0;
}
