# rvg_cache Report

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

