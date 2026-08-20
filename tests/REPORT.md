# rvg_cache Report: Fully-Automatic API

## Methodology

**Accuracy.** For each of the 11 supported distributions, `rvg_generate` is called once per seed from 1 to 10000, with no separate warmup pass -- its internal cache starts cold on seed 1 and matures live as the loop proceeds, exactly how it is actually used. Each call is compared against `generate_opt` run on a separate `gsl_rng` seeded identically: the two `double` results are compared via `memcpy` to `uint64_t` and a raw-bit comparison (never `==`, so NaN/signed-zero divergences would still be caught), and the two `flip_state.num_flips` counters are compared exactly, confirming the cached path consumed exactly the same amount of randomness as the uncached path.

**Interleaved-CDF correctness.** A separate test calls `rvg_generate` for three different CDFs in a round-robin loop over the same `prng` stream -- two different parameterizations of Gamma (shape=2,scale=3 and shape=0.5,scale=1) plus Poisson(71) -- 3000 rounds, 9000 total calls, each compared against `generate_opt` on a parallel reference stream. This specifically exercises the case the old manual API's abort-guard used to protect against (reusing one cache across different CDFs): under the new design there is no cache handle to misuse, since each distinct `cdf` function pointer automatically gets its own independent internal cache.

**Performance.** For each distribution, `generate_opt` and `rvg_generate` are each timed over N in {500, 1000, 2000} samples, seed=42, via `clock_gettime(CLOCK_MONOTONIC)`. Before each `rvg_generate` timing run, an internal (non-public) reset hook clears all auto-managed cache state, so every measurement starts genuinely cold -- matching the old methodology's "brand-new cache, no warmup" runs, just without a manual `rvg_cache_create` call to build that fresh state.

**Memory leak check.** A separate program calls `rvg_generate` ~20,000 times across 5 distributions (2 head-only, 3 head+tail) with zero manual cleanup calls anywhere, then exits normally. Run under macOS's `leaks --atExit` (code-signed with a debug entitlement so the tool can fully introspect the process) to confirm the `atexit`-registered cleanup actually runs and frees everything.

## Accuracy

| Distribution | Seeds Tested | Mismatches |
|---|---:|---:|
| BETA | 10000 | 0 |
| TDIST | 10000 | 0 |
| CHISQUARE | 10000 | 0 |
| FDIST | 10000 | 0 |
| GAMMA | 10000 | 0 |
| POISSON | 10000 | 0 |
| BINOMIAL | 10000 | 0 |
| NEGBINOMIAL | 10000 | 0 |
| GEOMETRIC | 10000 | 0 |
| HYPERGEOMETRIC | 10000 | 0 |
| PASCAL | 10000 | 0 |

**Interleaved-CDF correctness**: 3000 rounds x 3 CDFs (two different Gammas + Poisson) = 9000 calls, **0 mismatches**, no crash.

## Performance (speedup = baseline/cached, single run per N, seed=42)

| Distribution | N=500 | N=1000 | N=2000 |
|---|---:|---:|---:|
| BETA | 1.280x | 1.299x | 1.379x |
| TDIST | 1.279x | 1.278x | 1.279x |
| CHISQUARE | 1.247x | 1.279x | 1.279x |
| FDIST | 1.286x | 1.312x | 1.349x |
| GAMMA | 1.186x | 1.205x | 1.238x |
| POISSON | 6.643x | 9.763x | 16.332x |
| BINOMIAL | 8.488x | 13.671x | 20.629x |
| NEGBINOMIAL | 9.609x | 14.611x | 21.052x |
| GEOMETRIC | 2.074x | 3.012x | 3.866x |
| HYPERGEOMETRIC | 9.393x | 14.996x | 20.015x |
| PASCAL | 7.973x | 11.824x | 17.132x |

Same overall pattern as the previous design: continuous distributions (head-cache only) see modest but real gains (~1.2-1.4x), discrete distributions (head + tail cache) see large gains that grow with N (8-21x already by N=2000), since a tail-cache hit skips the rest of a sample's descent entirely rather than saving a single CDF call.

## Memory leak check

```
leaks Report Version: 4.0, multi-line stacks
Process 26353: 187 nodes malloced for 26 KB
Process 26353: 0 leaks for 0 total leaked bytes.
```
Zero leaks, with zero manual cleanup calls anywhere in the test program -- confirms the `atexit`-registered automatic cleanup runs correctly on normal process exit.

## Limitations

1. **No introspection.** There is no way to ask "how well is caching working right now" from outside `rvg_cache.c` -- `rvg_cache_stats` and the cache handle it operated on no longer exist in the public API. (Internally, `rvg_internal_stats_for_cdf` still exists for this project's own test programs to forward-declare and use, but it is deliberately not declared in `rvg_cache.h` and is not part of the supported API.)
2. **No independent caches for identical parameterizations.** Caching state is keyed on the exact `cdf` function pointer. Two different `cdf` functions for the same distribution+parameters (e.g. two separately-defined wrappers around the same GSL call with the same arguments) would each get their own cache -- there's no way to deliberately force two calls with the *same* `cdf` pointer to use independent caches.
3. **`dist`/`force_unsafe` only matter on the first call for a given `cdf`.** Later calls for that same `cdf` reuse whatever caching decision was made the first time, even if you pass different values afterward.
4. **Distributions not covered here**: GAUSSIAN, EXPONENTIAL, LOGNORMAL, CAUCHY, LAPLACE, LOGISTIC are tabulated `RVG_STATUS_NOT_RECOMMENDED` internally (caching measured, not shown to help); RAYLEIGH, PARETO, WEIBULL, GUMBEL are `RVG_STATUS_UNMEASURED`. `rvg_generate` quietly runs these uncached unless `force_unsafe` is set on their first call.
5. **Not thread-safe.** All caching state (including the registry mapping each `cdf` to its cache) is shared, unsynchronized, global state inside `rvg_cache.c`. Concurrent calls to `rvg_generate` from multiple threads are not supported.
6. **Single machine, single build.** Generated on `Darwin Shivens-Macbook.local 24.6.0 ... RELEASE_ARM64_T8112 arm64`, compiled with Apple Clang 17.0.0 (`gcc` resolves to `/usr/bin/gcc` on this machine) via `-O3 -DNDEBUG -march=native`. Not cross-validated on other hardware, OS, or compilers -- absolute and relative timings should not be assumed to transfer elsewhere.
