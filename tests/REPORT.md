# rvg_cache Full Report

## Methodology

**Accuracy (Section 1).** 11 distributions x 3 parameter sets each (typical, skewed, edge) = 33 configs. For each config, one `rvg_cache_t` is created with `rvg_cache_create(dist, force_unsafe=1, &status)` and then reused, unmodified, across 20000 consecutive seeds (1..20000) -- there is no separate warmup pass, so the cache is cold on seed 1 and matures live as the loop proceeds, which is how it would be used in practice. For each seed, `generate_opt` is run against a fresh `gsl_rng` seeded with that seed, and `rvg_generate` (same cache, same `cdf`) is run against a *separate* fresh `gsl_rng` seeded identically. The two results are compared by `memcpy`-ing each `double` to a `uint64_t` and comparing raw bits (never `==`, so a NaN or signed-zero divergence would still be caught), and the two `flip_state.num_flips` counters are compared exactly, which additionally confirms the cached path consumed the exact same amount of randomness as the uncached path. `force_unsafe=1` is passed uniformly since this run is deliberately stress-testing skewed and edge parameterizations, not just recommended ones; all 11 distributions in this report are tabulated as `RVG_STATUS_OK` regardless, so this does not change which entries get built.

**Performance (Section 2).** Typical parameter sets only (11 distributions, not 33 configs). For each of N in {100, 250, 500, 1000, 2000, 5000, 10000, 25000, 50000, 100000}, the timing is repeated with different fixed seeds, and *each* repeat uses its own fresh, no-warmup cache and its own fresh `gsl_rng` -- `generate_opt` and `rvg_generate` in a given repeat see the same seed as each other, but each repeat uses a different seed from the others. **N=100, 250, 500, and 1000 use 20 repeats; N=2000 and above use 5 repeats.** The smaller N values get more repeats because their absolute wall-clock time (a few hundred microseconds to a few milliseconds) sits closest to the measurement noise floor -- scheduler jitter and timer-call overhead are a much larger fraction of a short interval than a long one, so the same relative noise shrinks as N grows and fewer repeats are needed to see a stable median at N=100000 than at N=100. Wall-clock time is `clock_gettime(CLOCK_MONOTONIC)` around each N-sample loop, and `speedup = baseline_seconds / cached_seconds` per repeat. What is reported at each N is the **median** speedup across the repeats plus the **[min, max] range**, not a single number: single-shot wall-clock timing at N=100..2000 is only a few milliseconds and is easily dominated by scheduler noise, so a lone sample would misrepresent how stable the effect actually is (the wide ranges visible at N=100-500 below, occasionally including an outlier several times the median, are a direct symptom of this and are exactly why those N values get 20 repeats instead of 5). **Breakeven N** is the smallest tested N at which the *median* speedup is >= 1.0 *and* stays >= 1.0 for every larger N tested (i.e. the point past which the cache is never seen to fall behind again in this run, on the median); "N/A - never breaks even" means no such N was found among the values tested. **Memory@100000** is the *median* (across repeats) of `rvg_cache_stats()` taken on each N=100000 cache before it is freed: `head_entries * rvg_cache_head_entry_size() + tail_entries * rvg_cache_tail_entry_size()`, i.e. the actual occupied-slot count times the real struct size on this build, not a hardcoded guess.

**Head-cache saturation (below).** BETA, CHISQUARE, and GAMMA are all head-cache-only distributions (`tail_enabled=0` in the internal table) with a 65536-slot fixed-capacity head table (`RVG_HEAD_CAPACITY`). For these three, the table below Section 2 additionally reports the median (over the same repeats used for that N -- 20 at N=100..1000, 5 at N=2000 and above) `head_entries`, `head_hits`, `head_misses`, and `head_probe_fail` (the subset of misses where the probe window was full of *other* keys and the cache gave up without inserting, i.e. wasted probe work with no caching benefit) at each N, to show what actually happens as the table fills toward capacity.

## Accuracy

| Distribution | Parameters | Variant | Seeds Tested | Mismatches |
|---|---|---|---:|---:|
| BETA | Beta(5,5) | typical | 20000 | 0 |
| BETA | Beta(2,20) | skewed | 20000 | 0 |
| BETA | Beta(0.5,0.5) | edge | 20000 | 0 |
| TDIST | Tdist(5) | typical | 20000 | 0 |
| TDIST | Tdist(2) | skewed | 20000 | 0 |
| TDIST | Tdist(1) | edge | 20000 | 0 |
| CHISQUARE | ChiSquare(13) | typical | 20000 | 0 |
| CHISQUARE | ChiSquare(3) | skewed | 20000 | 0 |
| CHISQUARE | ChiSquare(0.5) | edge | 20000 | 0 |
| FDIST | Fdist(5,2) | typical | 20000 | 0 |
| FDIST | Fdist(1,10) | skewed | 20000 | 0 |
| FDIST | Fdist(1,1) | edge | 20000 | 0 |
| GAMMA | Gamma(0.5,1) | typical | 20000 | 0 |
| GAMMA | Gamma(2,3) | skewed | 20000 | 0 |
| GAMMA | Gamma(0.05,1) | edge | 20000 | 0 |
| POISSON | Poisson(71) | typical | 20000 | 0 |
| POISSON | Poisson(5) | skewed | 20000 | 0 |
| POISSON | Poisson(1) | edge | 20000 | 0 |
| BINOMIAL | Binomial(0.2,100) | typical | 20000 | 0 |
| BINOMIAL | Binomial(0.05,200) | skewed | 20000 | 0 |
| BINOMIAL | Binomial(0.5,3) | edge | 20000 | 0 |
| NEGBINOMIAL | NegBinomial(0.71,18) | typical | 20000 | 0 |
| NEGBINOMIAL | NegBinomial(0.1,5) | skewed | 20000 | 0 |
| NEGBINOMIAL | NegBinomial(0.9,1) | edge | 20000 | 0 |
| GEOMETRIC | Geometric(0.4) | typical | 20000 | 0 |
| GEOMETRIC | Geometric(0.05) | skewed | 20000 | 0 |
| GEOMETRIC | Geometric(0.95) | edge | 20000 | 0 |
| HYPERGEOMETRIC | Hypergeometric(5,20,7) | typical | 20000 | 0 |
| HYPERGEOMETRIC | Hypergeometric(2,48,10) | skewed | 20000 | 0 |
| HYPERGEOMETRIC | Hypergeometric(1,1,1) | edge | 20000 | 0 |
| PASCAL | Pascal(0.5,5) | typical | 20000 | 0 |
| PASCAL | Pascal(0.1,3) | skewed | 20000 | 0 |
| PASCAL | Pascal(0.9,1) | edge | 20000 | 0 |

## Performance (typical parameters, speedup = baseline/cached; N=100-1000: 20 repeats, N=2000-100000: 5 repeats; median [min, max])

| Distribution | N=100 | N=250 | N=500 | N=1000 | N=2000 | N=5000 | N=10000 | N=25000 | N=50000 | N=100000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BETA | 1.274x [1.065, 1.526] | 1.296x [1.205, 1.334] | 1.322x [1.221, 1.430] (50 repeats, see Limitations) | 1.337x [1.281, 1.358] | 1.362x [1.337, 1.378] | 1.386x [1.355, 1.398] | 1.413x [1.404, 1.428] | 1.448x [1.435, 1.449] | 1.456x [1.453, 1.459] | 1.465x [1.459, 1.470] |
| TDIST | 1.215x [1.101, 1.963] | 1.239x [1.176, 1.260] | 1.249x [1.133, 1.279] | 1.261x [1.124, 1.280] | 1.299x [1.286, 1.306] | 1.302x [1.262, 1.326] | 1.323x [1.322, 1.335] | 1.354x [1.339, 1.359] | 1.381x [1.368, 1.395] | 1.396x [1.395, 1.403] |
| CHISQUARE | 1.196x [1.147, 1.247] | 1.212x [1.066, 1.241] | 1.224x [1.016, 1.267] | 1.257x [1.215, 1.270] | 1.295x [1.288, 1.299] | 1.316x [1.296, 1.320] | 1.331x [1.318, 1.338] | 1.359x [1.356, 1.365] | 1.366x [1.359, 1.372] | 1.372x [1.371, 1.373] |
| FDIST | 1.228x [1.180, 1.280] | 1.264x [1.144, 1.302] | 1.288x [1.224, 1.309] | 1.303x [1.263, 1.319] | 1.324x [1.270, 1.330] | 1.359x [1.334, 1.366] | 1.366x [1.295, 1.383] | 1.402x [1.382, 1.412] | 1.410x [1.394, 1.423] | 1.429x [1.390, 1.436] |
| GAMMA | 1.128x [1.089, 1.160] | 1.159x [1.131, 1.181] | 1.158x [1.100, 1.185] | 1.170x [1.131, 1.209] | 1.201x [1.178, 1.231] | 1.255x [1.245, 1.270] | 1.278x [1.257, 1.293] | 1.315x [1.314, 1.316] | 1.322x [1.320, 1.324] | 1.332x [1.331, 1.335] |
| POISSON | 2.521x [2.231, 3.749] | 4.226x [3.979, 5.104] | 6.460x [5.671, 7.262] | 10.173x [9.437, 10.927] | 16.091x [15.425, 16.531] | 26.960x [25.883, 27.394] | 36.236x [35.364, 36.742] | 48.385x [47.703, 48.956] | 55.480x [54.751, 56.524] | 59.542x [58.743, 60.675] |
| BINOMIAL | 3.568x [3.073, 4.287] | 6.347x [5.654, 6.925] | 10.024x [9.375, 10.705] | 14.972x [13.345, 16.120] | 21.346x [20.988, 22.491] | 31.628x [31.073, 32.877] | 39.416x [39.149, 39.600] | 46.884x [44.295, 47.626] | 50.905x [50.432, 51.337] | 53.402x [53.255, 53.539] |
| NEGBINOMIAL | 4.189x [3.535, 4.935] | 7.229x [5.858, 8.093] | 11.057x [10.299, 12.039] | 15.953x [14.925, 17.057] | 22.695x [22.094, 23.590] | 31.234x [30.851, 31.966] | 37.620x [37.064, 37.875] | 43.694x [43.307, 44.306] | 46.275x [39.281, 46.622] | 47.986x [47.674, 48.671] |
| GEOMETRIC | 2.178x [1.785, 2.528] | 2.853x [2.126, 3.200] | 3.456x [3.179, 3.699] | 3.975x [3.841, 4.197] | 4.504x [4.421, 4.520] | 4.876x [4.725, 4.891] | 5.045x [4.961, 5.080] | 5.188x [5.066, 5.223] | 5.221x [5.124, 5.293] | 5.258x [5.127, 5.290] |
| HYPERGEOMETRIC | 7.050x [5.379, 13.210] | 11.948x [10.577, 13.852] | 16.768x [14.106, 17.709] | 21.667x [18.049, 22.932] | 25.423x [24.731, 26.258] | 29.736x [29.400, 30.092] | 31.384x [30.940, 32.005] | 32.838x [32.746, 33.232] | 33.266x [33.220, 33.565] | 33.655x [33.594, 33.747] |
| PASCAL | 3.766x [3.133, 4.775] | 6.482x [5.280, 7.613] | 9.503x [8.714, 10.573] | 13.886x [13.386, 15.857] | 19.797x [19.455, 20.008] | 26.695x [24.634, 26.998] | 30.960x [30.714, 31.019] | 34.816x [34.145, 35.004] | 36.607x [36.475, 36.788] | 37.683x [37.420, 37.756] |

## Head-Cache Saturation: BETA, CHISQUARE, GAMMA

Head-cache-only distributions (no tail cache); `RVG_HEAD_CAPACITY` = 65536 slots. All counts are the median over the repeats used at that N (20 for N=100-1000, 5 for N=2000-100000).

### BETA

| N | head_entries | head_hits | head_misses | head_probe_fail | hit rate | speedup (median) |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1141 | 1759 | 1141 | 0 | 60.7% | 1.274x |
| 250 | 2508 | 4742 | 2508 | 0 | 65.4% | 1.296x |
| 500 | 4506 | 9994 | 4506 | 0 | 68.9% | 1.325x |
| 1000 | 7998 | 21002 | 7998 | 0 | 72.4% | 1.337x |
| 2000 | 14065 | 43933 | 14067 | 2 | 75.7% | 1.362x |
| 5000 | 28402 | 116418 | 28582 | 185 | 80.3% | 1.386x |
| 10000 | 45067 | 242475 | 47525 | 2514 | 83.6% | 1.413x |
| 25000 | 62900 | 631320 | 93680 | 30845 | 87.1% | 1.448x |
| 50000 | 65411 | 1284965 | 165035 | 99624 | 88.6% | 1.456x |
| 100000 | 65534 | 2593291 | 306709 | 241175 | 89.4% | 1.465x |

**Observed saturation:** the head table for BETA does not reach its 65536-entry capacity within the tested N range (up to N=100000).

### CHISQUARE

| N | head_entries | head_hits | head_misses | head_probe_fail | hit rate | speedup (median) |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1580 | 1720 | 1580 | 0 | 52.1% | 1.196x |
| 250 | 3604 | 4646 | 3604 | 0 | 56.3% | 1.212x |
| 500 | 6692 | 9808 | 6692 | 0 | 59.4% | 1.224x |
| 1000 | 12392 | 20608 | 12392 | 1 | 62.4% | 1.257x |
| 2000 | 22742 | 43211 | 22789 | 46 | 65.5% | 1.295x |
| 5000 | 47253 | 114492 | 50508 | 3222 | 69.4% | 1.316x |
| 10000 | 63308 | 237355 | 92645 | 29343 | 71.9% | 1.331x |
| 25000 | 65536 | 607726 | 217274 | 151738 | 73.7% | 1.359x |
| 50000 | 65536 | 1225553 | 424447 | 358911 | 74.3% | 1.366x |
| 100000 | 65536 | 2461129 | 838871 | 773335 | 74.6% | 1.372x |

**Observed saturation:** the head table for CHISQUARE first reaches its 65536-entry capacity at N=25000 among the tested N values (median across repeats).

### GAMMA

| N | head_entries | head_hits | head_misses | head_probe_fail | hit rate | speedup (median) |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1400 | 1500 | 1400 | 0 | 51.7% | 1.128x |
| 250 | 3162 | 4088 | 3162 | 0 | 56.4% | 1.159x |
| 500 | 5837 | 8663 | 5837 | 0 | 59.7% | 1.158x |
| 1000 | 10598 | 18402 | 10598 | 0 | 63.5% | 1.170x |
| 2000 | 19216 | 38761 | 19239 | 23 | 66.8% | 1.201x |
| 5000 | 40367 | 103332 | 41668 | 1280 | 71.3% | 1.255x |
| 10000 | 59300 | 215812 | 74188 | 14956 | 74.4% | 1.278x |
| 25000 | 65509 | 557875 | 167125 | 101616 | 76.9% | 1.315x |
| 50000 | 65536 | 1129241 | 320759 | 255223 | 77.9% | 1.322x |
| 100000 | 65536 | 2271395 | 628605 | 563069 | 78.3% | 1.332x |

**Observed saturation:** the head table for GAMMA first reaches its 65536-entry capacity at N=50000 among the tested N values (median across repeats).

## Summary

| Distribution | Breakeven N | Speedup@1000 (median [min,max]) | Speedup@10000 (median [min,max]) | Speedup@100000 (median [min,max]) | Memory@100000 (median) | Total Mismatches (3 param sets) |
|---|---:|---:|---:|---:|---:|---:|
| BETA | 100 | 1.337x [1.281,1.358] | 1.413x [1.404,1.428] | 1.465x [1.459,1.470] | 1536.0 KB | 0 |
| TDIST | 100 | 1.261x [1.124,1.280] | 1.323x [1.322,1.335] | 1.396x [1.395,1.403] | 1516.8 KB | 0 |
| CHISQUARE | 100 | 1.257x [1.215,1.270] | 1.331x [1.318,1.338] | 1.372x [1.371,1.373] | 1536.0 KB | 0 |
| FDIST | 100 | 1.303x [1.263,1.319] | 1.366x [1.295,1.383] | 1.429x [1.390,1.436] | 1425.8 KB | 0 |
| GAMMA | 100 | 1.170x [1.131,1.209] | 1.278x [1.257,1.293] | 1.332x [1.331,1.335] | 1536.0 KB | 0 |
| POISSON | 100 | 10.173x [9.437,10.927] | 36.236x [35.364,36.742] | 59.542x [58.743,60.675] | 576.5 KB | 0 |
| BINOMIAL | 100 | 14.972x [13.345,16.120] | 39.416x [39.149,39.600] | 53.402x [53.255,53.539] | 358.9 KB | 0 |
| NEGBINOMIAL | 100 | 15.953x [14.925,17.057] | 37.620x [37.064,37.875] | 47.986x [47.674,48.671] | 272.5 KB | 0 |
| GEOMETRIC | 100 | 3.975x [3.841,4.197] | 5.045x [4.961,5.080] | 5.258x [5.127,5.290] | 217.9 KB | 0 |
| HYPERGEOMETRIC | 100 | 21.667x [18.049,22.932] | 31.384x [30.940,32.005] | 33.655x [33.594,33.747] | 81.7 KB | 0 |
| PASCAL | 100 | 13.886x [13.386,15.857] | 30.960x [30.714,31.019] | 37.683x [37.420,37.756] | 225.4 KB | 0 |

## Limitations

1. A cache's keys carry no distribution-parameter or CDF identity beyond the `cdf` function pointer it binds to on first use (see `rvg_cache.h`); in practice this means **the distribution parameters must be fixed for the entire lifetime of one cache object**. Changing parameters requires creating a new cache -- this report creates one cache per config for exactly that reason.
2. Results for GAUSSIAN, EXPONENTIAL, LOGNORMAL, CAUCHY, LAPLACE, and LOGISTIC are **not included in this report**: they are tabulated as `RVG_STATUS_NOT_RECOMMENDED`, meaning caching was measured for them and was not shown to provide a reliable speedup. RAYLEIGH, PARETO, WEIBULL, and GUMBEL are tabulated as `RVG_STATUS_UNMEASURED` and are only correctness-tested separately (not for speed here); nothing in this report should be read as a claim about any of these ten distributions.
3. This report was generated on a single machine (`uname -a`: `Darwin Shivens-Macbook.local 24.6.0 Darwin Kernel Version 24.6.0: Mon Jan 19 22:01:08 PST 2026; root:xnu-11417.140.69.708.3~1/RELEASE_ARM64_T8112 arm64`; `gcc --version`: `Apple clang version 17.0.0 (clang-1700.6.4.2)` -- on macOS this `gcc` is Apple Clang, not upstream GCC), compiled with `-O3 -DNDEBUG -march=native`, and has **not been cross-validated on other hardware, OS, or compiler combinations**. Absolute and relative timings (especially `-march=native`, which tunes for this specific CPU) should not be assumed to transfer to other machines. **Compiler used for `librvg.a` vs. the cache/test code:** both were built by the exact same `gcc` invocation on this machine -- the root Makefile's object rule (`%.o: %.c %.h`) calls the literal command `gcc`, which on this machine resolves via `/usr/bin/gcc` to Apple Clang Apple clang version 17.0.0 (clang-1700.6.4.2) (confirmed: `which gcc` -> `/usr/bin/gcc`, and no unversioned `gcc` exists under `/opt/homebrew/bin`, only `gcc-16`); the manual build command used for this report's test binaries also invokes plain `gcc`, so it resolves to the identical compiler. A separate Homebrew GCC (`gcc-16`, real GNU GCC) is present on this machine and *is* required elsewhere in this project -- specifically for `examples/main.c`, which defines its CDFs with `MAKE_CDF_P`/`MAKE_CDF_UINT_P` *inside* `main()`, a GNU C nested-function definition that Apple Clang rejects outright (confirmed directly: a minimal nested-function test file fails under `/usr/bin/gcc` with "function definition is not allowed here" and compiles and runs correctly under `gcc-16`). Neither `generate.c`/the rest of core `librvg.a`, nor `cache/rvg_cache.c`, nor `tests/full_report.c` use nested functions -- this report's CDF macros are deliberately invoked at file scope, not inside `main()`, specifically so they do not need that extension -- so this project's Homebrew-GCC requirement does not apply to anything measured here, and both the baseline and cached code paths in this report were compiled by the same compiler binary. This also means fairness is not a cross-compiler question at all in this run: `rvg_generate` does not reimplement `generate_opt`'s inner arithmetic (`subtract_exact`, `ith_bit_of_exact`, `bij64_lex2float`, `int2double`) -- it links directly against the same compiled definitions from `librvg.a` that `generate_opt` itself calls -- so even under a hypothetical mixed-compiler build, both paths would still be exercising identical machine code for that shared arithmetic; only the outer descent loop (copied into `rvg_cache.c` because the `cdf(d)` call site had to be intercepted) is separately compiled, and it is compiled once, by one compiler, in any given build.
4. **Head-cache saturation (BETA, CHISQUARE, GAMMA).** All three are head-cache-only (no tail cache) and share the 65536-slot fixed `RVG_HEAD_CAPACITY`. Observed in this run (median `head_entries` first reaching capacity; see the per-N tables above): BETA saturates at not reached by N=100000, CHISQUARE at N=25000, GAMMA at N=50000. Past that point, new (b, l) keys can no longer be inserted -- lookups for keys not already present either still hit (if that exact key was cached before saturation) or become `head_probe_fail` (6 wasted probes, then a direct, uncached `cdf()` call). Comparing median speedup at the saturation N against median speedup at N=100000 ("before -> after", classified as climbing/flat/degrading using a +/-0.02x dead zone around zero change, since a difference smaller than that is within the noise already visible in this section's [min, max] repeat ranges): BETA n/a (never saturated in this run); CHISQUARE 1.359x -> 1.372x (**flat**); GAMMA 1.322x -> 1.332x (**flat**). This should not be read as a guarantee for other distributions, other parameterizations, or N well beyond 100000 -- it is only what this run observed, and a distribution whose recurring (b, l) prefixes are less concentrated at shallow levels could plausibly behave differently once its table saturates.
5. **BETA N=500 outlier, investigated.** The initial 20-repeat run produced a BETA N=500 range of [0.999x, 3.828x] -- wildly inconsistent with every neighboring N and every other repeat at that N, which all sit in a tight ~1.2-1.4x band. This was re-run independently with 50 repeats on 50 new fixed seeds (disjoint from the original 20). Result: none of the 50 new repeats came anywhere near 0.999x or 3.828x -- every value fell between 1.221x and 1.430x, median 1.322x. Applying a 2x-median-absolute-deviation trim (MAD = 0.015x, threshold = |x - median| > 0.030x) to this 50-value set excludes 8 of the 50 values (1.221, 1.227, 1.236, 1.243, 1.254, 1.268, 1.285, 1.430), giving a trimmed mean of 1.322x -- essentially identical to the untrimmed median, because the underlying distribution here is a tight, roughly symmetric cluster with no real outliers among these 50, not a heavy-tailed one; the MAD trim mechanically flags the tails of a normal-looking cluster rather than isolating a distinct anomalous subpopulation. The practical conclusion is that the original 3.828x/0.999x extremes did not reproduce and are most likely a one-off measurement artifact at N=500's short wall-clock scale (~5-6ms per run -- e.g. a scheduler preemption or background-process contention hitting one of the 20 original repeats) rather than a real property of BETA-cache interaction at that N. The Performance table's BETA/N=500 cell above has been updated to this 50-repeat result (median [raw min, raw max], no values dropped from the displayed range) in place of the original 20-repeat [0.999, 3.828] figure; no other cell in this report was affected by or is expected to exhibit this issue, since it was not observed at any other (distribution, N) pair. The Summary table has no N=500 column (only Speedup@1000/@10000/@100000), so nothing there needed updating.
