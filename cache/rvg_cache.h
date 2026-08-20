/*
  Purpose:  Optional CDF-memoization layer for optimal random variate generation.
  Note:     This is an add-on. Does not modify the core library; `generate.c`,
            `generate.h` and every other pre-existing file are untouched.
  Thread-safety: If your program only ever does one thing at a time (true for
            most programs), this doesn't apply to you -- use rvg_generate
            normally and skip this note.

            If your program does run multiple parts of itself at the same
            time ("threads"), don't call rvg_generate from more than one of
            those parts. All caching state is managed internally and shared
            by every caller in the process, so concurrent calls from
            different threads are not safe.
*/

#ifndef RVG_CACHE_H
#define RVG_CACHE_H

#include "generate.h"
#include "flip.h"

/* Overview                                                            */
/*
  `rvg_generate` is a drop-in replacement for `generate_opt` that speeds up
  repeated sampling from the same CDF, e.g. drawing many samples in a loop
  when training a model. It is sample-for-sample identical to `generate_opt`
  given the same PRNG stream. Caching never changes which values you get,
  only how fast you get them.

  There is nothing to set up and nothing to tear down. The first call to
  `rvg_generate` for a given `cdf` builds its caching state automatically;
  every later call for that same `cdf` reuses it and gets faster. Different
  `cdf` functions, even two Gammas with different parameters, are kept
  separate, each with its own bounded memory cost (see
  RVG_AUTO_REGISTRY_CAPACITY in rvg_cache.c for the limit on how many
  distinct CDFs can be tracked at once). All of this is freed automatically
  when the program exits normally; nothing needs to be freed by hand.

  Two things are memoized per CDF, both internal to this file: a shallow
  cache of individual CDF evaluations, and, for a handful of discrete
  distributions, a cache of entire finished samples once no more randomness
  remains to be consumed for them. Both are fixed size and never grow
  without bound. Once a cache fills, it simply stops helping further;
  generation still produces correct results either way.
*/

/** Named distributions this speedup is available for. There is no
    auto-detection: you state which distribution `cdf` belongs to, and the
    library uses that to decide how to cache it, and whether caching is
    even likely to help. For some distributions it reliably isn't, and
    `rvg_generate` quietly skips caching for those unless `force_unsafe`
    is set. Auto-detecting the distribution family from an arbitrary `cdf`
    is future work. Questions and contributions welcome at
    tshiven@ucla.edu. */
typedef enum {
    RVG_DIST_BETA, RVG_DIST_TDIST, RVG_DIST_CHISQUARE, RVG_DIST_FDIST, RVG_DIST_GAMMA,
    RVG_DIST_POISSON, RVG_DIST_BINOMIAL, RVG_DIST_NEGBINOMIAL, RVG_DIST_GEOMETRIC,
    RVG_DIST_HYPERGEOMETRIC, RVG_DIST_PASCAL,
    RVG_DIST_GAUSSIAN, RVG_DIST_EXPONENTIAL, RVG_DIST_LOGNORMAL, RVG_DIST_CAUCHY,
    RVG_DIST_LAPLACE, RVG_DIST_LOGISTIC,
    RVG_DIST_RAYLEIGH, RVG_DIST_PARETO, RVG_DIST_WEIBULL, RVG_DIST_GUMBEL
} rvg_dist_t;

/** Generate one variate from `cdf`, exactly like `generate_opt(cdf, prng)`
    (same value, same `prng` bit consumption) -- but faster on repeated calls
    for the same `cdf`, since caching state is built and reused automatically.

    `dist` tells it which distribution `cdf` is, so it knows how to cache it.
    Some distributions (Gaussian, Exponential, and a few others) were
    measured not to reliably benefit from caching; for those, `rvg_generate`
    quietly samples uncached unless `force_unsafe` is nonzero, in which case
    it caches anyway. Correctness is identical either way -- `force_unsafe`
    only affects speed, never the sampled value.

    `dist`/`force_unsafe` are only read the first time a given `cdf` is
    passed in; later calls for that same `cdf` reuse whatever was decided
    then, regardless of what you pass after. */
double rvg_generate(cdf32_t cdf, struct flip_state *prng, rvg_dist_t dist, int force_unsafe);

#endif
