/*
  Name:     rvg_cache.h
  Purpose:  Optional CDF-memoization layer for optimal random variate generation.
  Note:     This is an add-on. It does not modify the core library; `generate.c`,
            `generate.h` and every other pre-existing file are untouched.
*/

#ifndef RVG_CACHE_H
#define RVG_CACHE_H

#include <stddef.h>

#include "generate.h"
#include "flip.h"

/* ------------------------------------------------------------------ */
/* Overview                                                            */
/* ------------------------------------------------------------------ */
/*
  `rvg_generate` is a drop-in replacement for `generate_opt` that memoizes
  work performed inside the sampling loop. It is sample-for-sample identical
  to `generate_opt` given the same PRNG stream, because both caches are pure
  memoizations of deterministic functions:

    * Head cache. The value `cdf(d)` evaluated at the midpoint of level `l`
      depends only on the pair (b, l) -- the current bit-string prefix and its
      length. Caching that float therefore changes nothing but the number of
      CDF evaluations. It is enabled only for the shallow levels
      `l <= head_depth`, where prefixes repeat often across samples.

    * Tail cache. Once the loop has consumed its final random bit, the
      remainder of the descent is fully determined by the state
      (b, l, cdf_l, cdf_r, ell): no further `flip()` is performed, so the
      final double is a deterministic function of that state. Entries are
      recorded only for the levels strictly below the last `flip()`, and a hit
      returns the stored double without consuming any randomness -- exactly
      as the uncached loop would have done.

  IMPORTANT: a cache is bound to one CDF. The keys do not include the identity
  of the `cdf32_t` function pointer, so passing two different CDFs to the same
  `rvg_cache_t` would otherwise yield incorrect results. This is enforced at
  runtime: the cache binds to whichever `cdf` it first sees in `rvg_generate`,
  and any later call with a different `cdf` pointer calls `abort()` rather
  than risk silently mixing entries from two distributions. Create one cache
  per distribution instance (including per parameter setting).

  Both tables are fixed-capacity, open-addressed and bounded-probe. They NEVER
  resize and NEVER evict. When a key's probe window is full the cache simply
  gives up for that key and falls through to the ordinary computation, so a
  saturated cache degrades to plain `generate_opt` rather than growing without
  bound. Memory use is therefore constant and known at construction time.
*/

/* Capacities are fixed at construction and are these compile-time constants.
   Both must remain powers of two (the index is a bit-mask of the hash). */
#define RVG_HEAD_CAPACITY 65536u  /* head-cache slots, ~1.5 MiB */
#define RVG_TAIL_CAPACITY 16384u  /* tail-cache slots, ~0.6 MiB */

/* Probe budgets for linear probing. On exhaustion the lookup/insert is
   abandoned; no resize and no eviction ever occurs. */
#define RVG_HEAD_MAX_PROBES 6
#define RVG_TAIL_MAX_PROBES 8

/** Opaque cache handle. */
typedef struct rvg_cache rvg_cache_t;

/** Named distributions recognized by the cache-parameter lookup table.
    There is no auto-detection: the caller states which distribution the CDF
    belongs to, and the table supplies the tuning constants. */
typedef enum {
    RVG_DIST_BETA, RVG_DIST_TDIST, RVG_DIST_CHISQUARE, RVG_DIST_FDIST, RVG_DIST_GAMMA,
    RVG_DIST_POISSON, RVG_DIST_BINOMIAL, RVG_DIST_NEGBINOMIAL, RVG_DIST_GEOMETRIC,
    RVG_DIST_HYPERGEOMETRIC, RVG_DIST_PASCAL,
    RVG_DIST_GAUSSIAN, RVG_DIST_EXPONENTIAL, RVG_DIST_LOGNORMAL, RVG_DIST_CAUCHY,
    RVG_DIST_LAPLACE, RVG_DIST_LOGISTIC,
    RVG_DIST_RAYLEIGH, RVG_DIST_PARETO, RVG_DIST_WEIBULL, RVG_DIST_GUMBEL
} rvg_dist_t;

/** Recommendation status attached to each distribution in the table.

    RVG_STATUS_OK
        Caching has been measured for this distribution and the tabulated
        parameters were the configuration used. `rvg_cache_create` builds a
        cache for these without complaint.

    RVG_STATUS_NOT_RECOMMENDED
        Caching has been measured for this distribution and was NOT shown to
        provide a reliable speedup. These CDFs are cheap and/or their descents
        share too few prefixes for memoization to pay for itself; the cache can
        cost more than it saves. `rvg_cache_create` refuses unless
        `force_unsafe` is set.

    RVG_STATUS_UNMEASURED
        No measurement exists for this distribution. Nothing is claimed about
        it either way, and it likewise has NOT been shown to provide a reliable
        speedup. `rvg_cache_create` refuses unless `force_unsafe` is set.

    Only RVG_STATUS_OK entries carry any speedup claim. NOT_RECOMMENDED and
    UNMEASURED entries have not been shown to provide a reliable speedup, and
    forcing them is a benchmarking affordance, not an endorsement. Correctness
    is unaffected in every case: the sampled values are identical with or
    without a cache. */
typedef enum { RVG_STATUS_OK, RVG_STATUS_NOT_RECOMMENDED, RVG_STATUS_UNMEASURED } rvg_status_t;

/** Create a cache configured from the built-in table for `dist`.

    On success returns a new cache and, if `status` is non-NULL, writes the
    distribution's tabulated status to it.

    If the distribution's status is not RVG_STATUS_OK, returns NULL and writes
    that status through `status` (when non-NULL) -- UNLESS `force_unsafe` is
    non-zero, in which case a cache is built anyway using the table's defaults
    for that distribution, and `*status` still reports the tabulated status.

    Returns NULL (with `*status` left as RVG_STATUS_UNMEASURED) if `dist` is
    out of range or allocation fails. */
rvg_cache_t *rvg_cache_create(rvg_dist_t dist, int force_unsafe, rvg_status_t *status);

/** Release a cache. Passing NULL is a no-op. */
void rvg_cache_free(rvg_cache_t *);

/** Generate one variate from `cdf`, using `cache` to memoize the descent.

    Behaves exactly like `generate_opt(cdf, prng)`, including its consumption
    of the `prng` bit stream. Passing a NULL `cache` is legal and performs an
    uncached descent -- an intentional convenience for A/B use, not part of
    the original cache-key spec.

    `cache` binds to the first `cdf` it is called with and aborts the process
    if a later call passes a different `cdf`; see the CDF-binding note above. */
double rvg_generate(cdf32_t cdf, struct flip_state *prng, rvg_cache_t *cache);

/** Cumulative counters, useful for judging whether a cache is earning its
    keep. `*_entries` are the number of occupied slots in each table.
    `head_probe_fail` is the subset of `head_misses` where every slot in the
    key's probe window was occupied by a *different* key -- i.e. the head
    cache gave up without inserting and the CDF was called directly. It rises
    once the head table nears its fixed capacity (`RVG_HEAD_CAPACITY`) and
    distinct (b, l) keys keep arriving with nowhere left to go. */
typedef struct { size_t head_hits, head_misses, tail_hits, tail_misses, head_entries, tail_entries, head_probe_fail; } rvg_stats_t;

/** Snapshot the counters. Returns an all-zero struct for a NULL cache. */
rvg_stats_t rvg_cache_stats(rvg_cache_t *);

/** Byte size of one head-cache / tail-cache slot (occupied or not), for
    turning `rvg_cache_stats()` entry counts into an actual memory estimate:
    `head_entries * rvg_cache_head_entry_size() + tail_entries * rvg_cache_tail_entry_size()`.
    These reflect the real struct layout on the build in use rather than a
    guess baked into caller code. */
size_t rvg_cache_head_entry_size(void);
size_t rvg_cache_tail_entry_size(void);

#endif
