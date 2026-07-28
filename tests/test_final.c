/*
  Name:     test_final.c
  Purpose:  Correctness (cold-start, no warmup) and performance validation
            for the rvg_cache layer, across all 11 distributions it supports.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O3 -DNDEBUG -march=native $(gsl-config --cflags) \
                  -I.. tests/test_final.c cache/rvg_cache.c librvg.a \
                  -o tests/test_final.out $(gsl-config --libs-without-cblas) -lgmp
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
/* CDFs under test, matching examples/main.c's MAKE_CDF_* convention.  */
/* ------------------------------------------------------------------ */

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
    rvg_dist_t  dist;
    cdf32_t     cdf;
} dist_entry_t;

static dist_entry_t DISTS[] = {
    {"BETA",           RVG_DIST_BETA,           beta_cdf},
    {"TDIST",          RVG_DIST_TDIST,          tdist_cdf},
    {"CHISQUARE",      RVG_DIST_CHISQUARE,      chisq_cdf},
    {"FDIST",          RVG_DIST_FDIST,          fdist_cdf},
    {"GAMMA",          RVG_DIST_GAMMA,          gamma_cdf},
    {"POISSON",        RVG_DIST_POISSON,        poisson_cdf},
    {"BINOMIAL",       RVG_DIST_BINOMIAL,       binomial_cdf},
    {"NEGBINOMIAL",    RVG_DIST_NEGBINOMIAL,    negbinomial_cdf},
    {"GEOMETRIC",      RVG_DIST_GEOMETRIC,      geometric_cdf},
    {"HYPERGEOMETRIC", RVG_DIST_HYPERGEOMETRIC, hypergeometric_cdf},
    {"PASCAL",         RVG_DIST_PASCAL,         pascal_cdf},
};

#define NUM_DISTS ((int)(sizeof(DISTS) / sizeof(DISTS[0])))
#define NUM_SEEDS 10000

static const int SWEEP_N[] = {500, 1000, 2000};
#define NUM_SWEEP_N ((int)(sizeof(SWEEP_N) / sizeof(SWEEP_N[0])))

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
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
    rvg_status_t status;
    rvg_cache_t *cache = rvg_cache_create(dist, 0, &status);
    if (cache == NULL) { return -1.0; }

    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, seed);
    struct flip_state p = make_flip_state(r);

    double t0 = now_seconds();
    for (int i = 0; i < n; i++) {
        rvg_generate(cdf, &p, cache);
    }
    double t1 = now_seconds();

    gsl_rng_free(r);
    rvg_cache_free(cache);
    return t1 - t0;
}

/* ------------------------------------------------------------------ */
/* Main                                                                 */
/* ------------------------------------------------------------------ */

int main(void) {

    long cold_mismatches[NUM_DISTS];
    double speedups[NUM_DISTS][NUM_SWEEP_N];
    int usable[NUM_DISTS];

    for (int i = 0; i < NUM_DISTS; i++) {
        cold_mismatches[i] = -1;
        usable[i] = 0;
        for (int j = 0; j < NUM_SWEEP_N; j++) { speedups[i][j] = -1.0; }
    }

    printf("=== PART 1: Cold-start exactness ===\n\n");

    for (int i = 0; i < NUM_DISTS; i++) {
        dist_entry_t *d = &DISTS[i];

        rvg_status_t status;
        rvg_cache_t *cache = rvg_cache_create(d->dist, 0, &status);
        if (cache == NULL) {
            printf("%s: rvg_cache_create did NOT return OK (status=%d) -- skipping.\n",
                   d->name, (int)status);
            continue;
        }
        usable[i] = 1;

        long mismatches = 0;
        for (unsigned long seed = 1; seed <= NUM_SEEDS; seed++) {

            gsl_rng *r1 = gsl_rng_alloc(gsl_rng_default);
            gsl_rng_set(r1, seed);
            struct flip_state p1 = make_flip_state(r1);
            double v1 = generate_opt(d->cdf, &p1);
            unsigned long f1 = p1.num_flips;
            gsl_rng_free(r1);

            gsl_rng *r2 = gsl_rng_alloc(gsl_rng_default);
            gsl_rng_set(r2, seed);
            struct flip_state p2 = make_flip_state(r2);
            double v2 = rvg_generate(d->cdf, &p2, cache);
            unsigned long f2 = p2.num_flips;
            gsl_rng_free(r2);

            if (!double_bits_eq(v1, v2) || (f1 != f2)) {
                mismatches++;
                printf("  MISMATCH seed=%lu dist=%s baseline=%.17g (flips=%lu) "
                       "cached=%.17g (flips=%lu)\n",
                       seed, d->name, v1, f1, v2, f2);
            }
        }

        rvg_cache_free(cache);
        cold_mismatches[i] = mismatches;
        printf("%s COLD: %ld/%d mismatches\n", d->name, mismatches, NUM_SEEDS);
    }

    printf("\n=== PART 2: Realistic-usage performance sweep ===\n\n");

    for (int i = 0; i < NUM_DISTS; i++) {
        if (!usable[i]) { continue; }
        dist_entry_t *d = &DISTS[i];

        for (int j = 0; j < NUM_SWEEP_N; j++) {
            int n = SWEEP_N[j];

            double baseline = run_baseline(d->cdf, 42, n);
            double cached = run_cached(d->cdf, d->dist, 42, n);
            double speedup = (cached > 0.0) ? (baseline / cached) : -1.0;
            speedups[i][j] = speedup;

            printf("%s N=%d: baseline=%.6fs cached=%.6fs speedup=%.3fx\n",
                   d->name, n, baseline, cached, speedup);
        }
    }

    printf("\n=== PART 3: Summary ===\n\n");
    printf("%-16s | %-18s | %10s | %10s | %10s\n",
           "distribution", "cold mismatches", "N=500", "N=1000", "N=2000");
    printf("-----------------------------------------------------------------------\n");

    int all_zero = 1;
    for (int i = 0; i < NUM_DISTS; i++) {
        dist_entry_t *d = &DISTS[i];
        if (!usable[i]) {
            printf("%-16s | %-18s | %10s | %10s | %10s\n",
                   d->name, "N/A (not OK)", "-", "-", "-");
            all_zero = 0;
            continue;
        }
        if (cold_mismatches[i] != 0) { all_zero = 0; }
        printf("%-16s | %18ld | %9.3fx | %9.3fx | %9.3fx\n",
               d->name, cold_mismatches[i],
               speedups[i][0], speedups[i][1], speedups[i][2]);
    }

    printf("\n");
    if (all_zero) {
        printf("ALL TESTS PASSED\n");
        return 0;
    } else {
        printf("TESTS FAILED\n");
        return 1;
    }
}
