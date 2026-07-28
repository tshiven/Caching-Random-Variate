/*
  Name:     full_report_lowN.c
  Purpose:  Recompute only N in {100, 250, 500, 1000} with 20 repeats (instead
            of the 5 repeats used for N>=2000 in full_report.c), for tighter
            confidence where wall-clock time is closest to the measurement
            noise floor. Prints values to be spliced into tests/REPORT.md by
            hand -- this file does not write REPORT.md itself.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O3 -DNDEBUG -march=native $(gsl-config --cflags) \
                  -I. -Icache -I<gmp include> -L<gmp lib> \
                  tests/full_report_lowN.c cache/rvg_cache.c librvg.a \
                  -o tests/full_report_lowN.out $(gsl-config --libs-without-cblas) -lgmp
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

/* Same 11 typical-parameter CDFs as full_report.c / test_final.c. */
MAKE_CDF_P(beta_typical_cdf, gsl_cdf_beta_P, 5, 5)
MAKE_CDF_P(tdist_typical_cdf, gsl_cdf_tdist_P, 5)
MAKE_CDF_P(chisq_typical_cdf, gsl_cdf_chisq_P, 13)
MAKE_CDF_P(fdist_typical_cdf, gsl_cdf_fdist_P, 5, 2)
MAKE_CDF_P(gamma_typical_cdf, gsl_cdf_gamma_P, 0.5, 1)
MAKE_CDF_UINT_P(poisson_typical_cdf, gsl_cdf_poisson_P, 71)
MAKE_CDF_UINT_P(binomial_typical_cdf, gsl_cdf_binomial_P, 0.2, 100)
MAKE_CDF_UINT_P(negbinomial_typical_cdf, gsl_cdf_negative_binomial_P, 0.71, 18)
MAKE_CDF_UINT_P(geometric_typical_cdf, gsl_cdf_geometric_P, 0.4)
MAKE_CDF_UINT_P(hypergeometric_typical_cdf, gsl_cdf_hypergeometric_P, 5, 20, 7)
MAKE_CDF_UINT_P(pascal_typical_cdf, gsl_cdf_pascal_P, 0.5, 5)

typedef struct {
    const char *name;
    rvg_dist_t  dist;
    cdf32_t     cdf;
} base_dist_t;

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

static const long N_LIST[] = {100, 250, 500, 1000};
#define NUM_N ((int)(sizeof(N_LIST) / sizeof(N_LIST[0])))

#define REPEATS 20
static const unsigned long PERF_SEEDS[REPEATS] = {
    42, 1337, 271828, 986532, 555555,
    7, 77, 777, 7777, 77777,
    3, 33, 333, 3333, 33333,
    9, 99, 999, 9999, 99999
};

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

typedef struct {
    double      baseline_sec;
    double      cached_sec;
    double      speedup;
    rvg_stats_t stats;
} perf_point_t;

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

static void sortN(double *v, int n) {
    for (int i = 1; i < n; i++) {
        double key = v[i];
        int j = i - 1;
        while (j >= 0 && v[j] > key) { v[j + 1] = v[j]; j--; }
        v[j + 1] = key;
    }
}

static double medianN(const double *raw, int n) {
    double tmp[REPEATS];
    memcpy(tmp, raw, (size_t)n * sizeof(double));
    sortN(tmp, n);
    return (n % 2 == 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}

typedef struct {
    double speedup_median, speedup_min, speedup_max;
    double head_hits_median, head_misses_median, head_probe_fail_median, head_entries_median;
} perf_summary_t;

static perf_summary_t summarize_perf(const perf_point_t *raw) {
    perf_summary_t s;
    double speedups[REPEATS], hh[REPEATS], hm[REPEATS], hpf[REPEATS], he[REPEATS];

    s.speedup_min = raw[0].speedup;
    s.speedup_max = raw[0].speedup;
    for (int k = 0; k < REPEATS; k++) {
        speedups[k] = raw[k].speedup;
        hh[k] = (double)raw[k].stats.head_hits;
        hm[k] = (double)raw[k].stats.head_misses;
        hpf[k] = (double)raw[k].stats.head_probe_fail;
        he[k] = (double)raw[k].stats.head_entries;
        if (raw[k].speedup < s.speedup_min) { s.speedup_min = raw[k].speedup; }
        if (raw[k].speedup > s.speedup_max) { s.speedup_max = raw[k].speedup; }
    }

    s.speedup_median = medianN(speedups, REPEATS);
    s.head_hits_median = medianN(hh, REPEATS);
    s.head_misses_median = medianN(hm, REPEATS);
    s.head_probe_fail_median = medianN(hpf, REPEATS);
    s.head_entries_median = medianN(he, REPEATS);
    return s;
}

int main(void) {
    printf("=== full_report_lowN: N in {100,250,500,1000}, %d repeats ===\n\n", REPEATS);

    perf_summary_t results[NUM_BASE_DISTS][NUM_N];

    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        base_dist_t *bd = &BASE_DISTS[i];
        for (int j = 0; j < NUM_N; j++) {
            perf_point_t raw[REPEATS];
            for (int k = 0; k < REPEATS; k++) {
                raw[k] = run_perf(bd->cdf, bd->dist, N_LIST[j], PERF_SEEDS[k]);
            }
            perf_summary_t summ = summarize_perf(raw);
            results[i][j] = summ;
            double total = summ.head_hits_median + summ.head_misses_median;
            double hit_rate = (total > 0.0) ? (summ.head_hits_median / total) : 0.0;
            printf("PERF %-16s N=%-6ld speedup_median=%.3f speedup_min=%.3f speedup_max=%.3f "
                   "head_entries=%.0f head_hits=%.0f head_misses=%.0f head_probe_fail=%.0f hit_rate=%.1f\n",
                   bd->name, N_LIST[j], summ.speedup_median, summ.speedup_min, summ.speedup_max,
                   summ.head_entries_median, summ.head_hits_median, summ.head_misses_median,
                   summ.head_probe_fail_median, hit_rate * 100.0);
        }
    }

    printf("\n=== Performance table cells (Distribution | N=100 | N=250 | N=500 | N=1000) ===\n\n");
    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        printf("%-16s", BASE_DISTS[i].name);
        for (int j = 0; j < NUM_N; j++) {
            perf_summary_t *s = &results[i][j];
            printf(" | %.3fx [%.3f, %.3f]", s->speedup_median, s->speedup_min, s->speedup_max);
        }
        printf("\n");
    }

    printf("\n=== Summary Speedup@1000 column ===\n\n");
    for (int i = 0; i < NUM_BASE_DISTS; i++) {
        perf_summary_t *s = &results[i][3]; /* N_LIST[3] == 1000 */
        printf("%-16s %.3fx [%.3f,%.3f]\n", BASE_DISTS[i].name, s->speedup_median, s->speedup_min, s->speedup_max);
    }

    return 0;
}
