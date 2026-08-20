/*
  Name:     breakeven_smalln.c
  Purpose:  Find the smallest N at which rvg_generate starts outperforming
            generate_opt, down to N=1. A naive version of this measurement
            has a cold-start artifact: the very first generate_opt call for
            a given cdf in the whole process pays a one-time cost (page
            faults, branch-predictor warmup, symbol resolution) that
            inflates the first reading at the smallest N tested, biasing the
            reported breakeven point upward.

            Fix: prime each distribution's code paths (both generate_opt
            and rvg_generate) with a throwaway run on a dedicated RNG/seed
            before any timed measurement for that distribution, then reset
            the auto-managed cache so the real (timed) rvg_generate calls
            still start genuinely cold.

            Also reports raw absolute seconds (not just the ratio) for a
            couple of representative N values, to make explicit how small
            these measurements are and where measurement noise dominates.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O3 -DNDEBUG -march=native $(gsl-config --cflags) \
                  -I. -Icache tests/breakeven_smalln.c cache/rvg_cache.c librvg.a \
                  -o tests/breakeven_smalln.out $(gsl-config --libs-without-cblas) -lgmp
            (run from the repository root; add -I/-L flags for gmp.h if it
            isn't already on your include path, e.g. via `brew --prefix gmp`).
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

extern void rvg_internal_reset_for_testing(void);

MAKE_CDF_P(beta_cdf, gsl_cdf_beta_P, 5, 5)
MAKE_CDF_P(tdist_cdf, gsl_cdf_tdist_P, 5)
MAKE_CDF_P(chisq_cdf, gsl_cdf_chisq_P, 13)
MAKE_CDF_P(fdist_cdf, gsl_cdf_fdist_P, 5, 2)
MAKE_CDF_P(gamma_cdf, gsl_cdf_gamma_P, 0.5, 1)

MAKE_CDF_UINT_P(poisson_cdf, gsl_cdf_poisson_P, 71)
MAKE_CDF_UINT_P(binomial_cdf, gsl_cdf_binomial_P, 0.2, 100)
MAKE_CDF_UINT_P(negbinomial_cdf, gsl_cdf_negative_binomial_P, 0.71, 18)
MAKE_CDF_UINT_P(geometric_cdf, gsl_cdf_geometric_P, 0.4)
MAKE_CDF_UINT_P(hypergeometric_cdf, gsl_cdf_hypergeometric_P, 5, 20, 7)
MAKE_CDF_UINT_P(pascal_cdf, gsl_cdf_pascal_P, 0.5, 5)

typedef struct {
    const char *name;
    const char *param_label;
    rvg_dist_t  dist;
    cdf32_t     cdf;
} dist_entry_t;

static dist_entry_t DISTS[] = {
    {"BETA",           "Beta(5,5)",          RVG_DIST_BETA,           beta_cdf},
    {"TDIST",          "Tdist(5)",           RVG_DIST_TDIST,          tdist_cdf},
    {"CHISQUARE",      "ChiSquare(13)",      RVG_DIST_CHISQUARE,      chisq_cdf},
    {"FDIST",          "Fdist(5,2)",         RVG_DIST_FDIST,          fdist_cdf},
    {"GAMMA",          "Gamma(0.5,1)",       RVG_DIST_GAMMA,          gamma_cdf},
    {"POISSON",        "Poisson(71)",        RVG_DIST_POISSON,        poisson_cdf},
    {"BINOMIAL",       "Binomial(0.2,100)",  RVG_DIST_BINOMIAL,       binomial_cdf},
    {"NEGBINOMIAL",    "NegBinomial(0.71,18)", RVG_DIST_NEGBINOMIAL,  negbinomial_cdf},
    {"GEOMETRIC",      "Geometric(0.4)",     RVG_DIST_GEOMETRIC,      geometric_cdf},
    {"HYPERGEOMETRIC", "Hypergeometric(5,20,7)", RVG_DIST_HYPERGEOMETRIC, hypergeometric_cdf},
    {"PASCAL",         "Pascal(0.5,5)",      RVG_DIST_PASCAL,         pascal_cdf},
};
#define NUM_DISTS ((int)(sizeof(DISTS) / sizeof(DISTS[0])))

/* Below 50, plus a few checkpoints up to 2000 to connect with the previous
   report's range. */
static const int SWEEP_N[] = {1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 75, 100, 200, 500, 2000};
#define NUM_SWEEP_N ((int)(sizeof(SWEEP_N) / sizeof(SWEEP_N[0])))

static const unsigned long SEEDS[] = {42, 1337, 271828, 986532, 555555};
#define NUM_SEEDS ((int)(sizeof(SEEDS) / sizeof(SEEDS[0])))

#define PRIME_CALLS 5000
#define PRIME_SEED 999999999UL /* never used as a real measurement seed */

/* ------------------------------------------------------------------ */

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Warm up code paths (page faults, branch predictor, symbol resolution,
   CPU frequency ramp) for this cdf on both the baseline and cached paths,
   with a seed never used for real measurements, so the very first *timed*
   call for this distribution isn't paying a one-time cold-start cost.
   Cache state is reset afterward so timed rvg_generate calls still start
   genuinely cold, unaffected by this warmup. */
static void prime(cdf32_t cdf, rvg_dist_t dist) {
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, PRIME_SEED);
    struct flip_state p = make_flip_state(r);
    for (int i = 0; i < PRIME_CALLS; i++) {
        generate_opt(cdf, &p);
        rvg_generate(cdf, &p, dist, 0);
    }
    gsl_rng_free(r);
    rvg_internal_reset_for_testing();
}

static double run_baseline(cdf32_t cdf, unsigned long seed, int n) {
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, seed);
    struct flip_state p = make_flip_state(r);

    double t0 = now_seconds();
    for (int i = 0; i < n; i++) {
        generate_opt(cdf, &p);
    }
    double t1 = now_seconds();

    gsl_rng_free(r);
    return t1 - t0;
}

static double run_cached(cdf32_t cdf, rvg_dist_t dist, unsigned long seed, int n) {
    rvg_internal_reset_for_testing();

    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, seed);
    struct flip_state p = make_flip_state(r);

    double t0 = now_seconds();
    for (int i = 0; i < n; i++) {
        rvg_generate(cdf, &p, dist, 0);
    }
    double t1 = now_seconds();

    gsl_rng_free(r);
    return t1 - t0;
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

static double median_n(const double *v, int n) {
    double tmp[NUM_SEEDS];
    memcpy(tmp, v, (size_t)n * sizeof(double));
    qsort(tmp, (size_t)n, sizeof(double), cmp_double);
    return n % 2 == 1 ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}

int main(void) {
    static double speedups[NUM_DISTS][NUM_SWEEP_N][NUM_SEEDS];
    static double baseline_secs[NUM_DISTS][NUM_SWEEP_N][NUM_SEEDS];
    static double cached_secs[NUM_DISTS][NUM_SWEEP_N][NUM_SEEDS];

    for (int d = 0; d < NUM_DISTS; d++) {
        prime(DISTS[d].cdf, DISTS[d].dist);
        for (int ni = 0; ni < NUM_SWEEP_N; ni++) {
            int n = SWEEP_N[ni];
            for (int si = 0; si < NUM_SEEDS; si++) {
                unsigned long seed = SEEDS[si];
                double tb = run_baseline(DISTS[d].cdf, seed, n);
                double tc = run_cached(DISTS[d].cdf, DISTS[d].dist, seed, n);
                baseline_secs[d][ni][si] = tb;
                cached_secs[d][ni][si] = tc;
                speedups[d][ni][si] = tb / tc;
            }
        }
    }

    printf("## Raw per-N speedups, primed (baseline_seconds / cached_seconds, 5 repeats)\n\n");
    for (int d = 0; d < NUM_DISTS; d++) {
        printf("### %s -- %s\n\n", DISTS[d].name, DISTS[d].param_label);
        printf("| N | seed=42 | seed=1337 | seed=271828 | seed=986532 | seed=555555 | median | min | max |\n");
        printf("|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        for (int ni = 0; ni < NUM_SWEEP_N; ni++) {
            double *v = speedups[d][ni];
            double med = median_n(v, NUM_SEEDS);
            double mn = v[0], mx = v[0];
            for (int si = 1; si < NUM_SEEDS; si++) {
                if (v[si] < mn) mn = v[si];
                if (v[si] > mx) mx = v[si];
            }
            printf("| %d | %.3fx | %.3fx | %.3fx | %.3fx | %.3fx | %.3fx | %.3fx | %.3fx |\n",
                   SWEEP_N[ni], v[0], v[1], v[2], v[3], v[4], med, mn, mx);
        }
        printf("\n");
    }

    /* Absolute times (microseconds) at a few representative N, to show the
       scale of what's being measured and where noise dominates signal. */
    printf("## Absolute timings in microseconds (median of 5 seeds), primed\n\n");
    printf("| Distribution | N | baseline us | cached us | delta us |\n");
    printf("|---|---:|---:|---:|---:|\n");
    static const int SHOW_N[] = {1, 5, 10, 20, 50, 100};
    for (int d = 0; d < NUM_DISTS; d++) {
        for (int si2 = 0; si2 < (int)(sizeof(SHOW_N)/sizeof(SHOW_N[0])); si2++) {
            int want = SHOW_N[si2];
            for (int ni = 0; ni < NUM_SWEEP_N; ni++) {
                if (SWEEP_N[ni] != want) continue;
                double bmed = median_n(baseline_secs[d][ni], NUM_SEEDS) * 1e6;
                double cmed = median_n(cached_secs[d][ni], NUM_SEEDS) * 1e6;
                printf("| %s | %d | %.2f | %.2f | %.2f |\n", DISTS[d].name, want, bmed, cmed, bmed - cmed);
            }
        }
    }
    printf("\n");

    printf("## Breakeven summary (primed, N from 1 upward)\n\n");
    printf("| Distribution | Median breakeven N | Conservative breakeven N |\n");
    printf("|---|---:|---:|\n");
    for (int d = 0; d < NUM_DISTS; d++) {
        int median_breakeven = -1;
        int conservative_breakeven = -1;

        int median_ok_from_here = 1;
        int conservative_ok_from_here = 1;
        for (int ni = NUM_SWEEP_N - 1; ni >= 0; ni--) {
            double *v = speedups[d][ni];
            double med = median_n(v, NUM_SEEDS);
            double mn = v[0];
            for (int si = 1; si < NUM_SEEDS; si++) if (v[si] < mn) mn = v[si];

            if (med >= 1.0 && median_ok_from_here) {
                median_breakeven = SWEEP_N[ni];
            } else {
                median_ok_from_here = 0;
            }

            if (mn >= 1.0 && conservative_ok_from_here) {
                conservative_breakeven = SWEEP_N[ni];
            } else {
                conservative_ok_from_here = 0;
            }
        }

        if (median_breakeven < 0) {
            printf("| %s | never (<1.0x at all tested N) | ", DISTS[d].name);
        } else {
            printf("| %s | %d | ", DISTS[d].name, median_breakeven);
        }
        if (conservative_breakeven < 0) {
            printf("never (<1.0x at all tested N) |\n");
        } else {
            printf("%d |\n", conservative_breakeven);
        }
    }

    return 0;
}
