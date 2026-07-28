/*
  Name:     beta_n500_investigate.c
  Purpose:  Investigate the BETA N=500 outlier seen in the 20-repeat run
            (range [0.999, 3.828] versus neighboring N values all sitting in
            a tight ~1.3x band). Re-runs BETA N=500 with 50 repeats and 50
            new fixed seeds, prints every individual speedup value sorted,
            and computes both the median and a MAD-trimmed mean.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O3 -DNDEBUG -march=native $(gsl-config --cflags) \
                  -I. -Icache -I<gmp include> -L<gmp lib> \
                  tests/beta_n500_investigate.c cache/rvg_cache.c librvg.a \
                  -o tests/beta_n500_investigate.out $(gsl-config --libs-without-cblas) -lgmp
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

MAKE_CDF_P(beta_typical_cdf, gsl_cdf_beta_P, 5, 5)

#define N_SAMPLES 500
#define REPEATS 50

/* 50 new fixed seeds, distinct from the {42,1337,271828,986532,555555,
   7,77,777,7777,77777,3,33,333,3333,33333,9,99,999,9999,99999} used in the
   earlier 20-repeat run, so this is genuinely independent evidence. */
static const unsigned long SEEDS[REPEATS] = {
    101, 202, 303, 404, 505, 606, 707, 808, 909, 1010,
    1111, 2222, 3333, 4444, 6666, 8888, 12121, 13131, 14141, 15151,
    21212, 23232, 25252, 27272, 29292, 31313, 33333333, 41414, 43434, 45454,
    51515, 53535, 55555555, 61616, 63636, 65656, 71717, 73737, 75757, 77777777,
    81818, 83838, 85858, 91919, 93939, 95959, 111213, 141516, 171819, 202122
};

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double run_one(unsigned long seed) {
    gsl_rng *r1 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r1, seed);
    struct flip_state p1 = make_flip_state(r1);
    double t0 = now_seconds();
    for (int i = 0; i < N_SAMPLES; i++) { generate_opt(beta_typical_cdf, &p1); }
    double t1 = now_seconds();
    double baseline_sec = t1 - t0;
    gsl_rng_free(r1);

    rvg_status_t status;
    rvg_cache_t *cache = rvg_cache_create(RVG_DIST_BETA, 1, &status);
    gsl_rng *r2 = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r2, seed);
    struct flip_state p2 = make_flip_state(r2);
    double t2 = now_seconds();
    for (int i = 0; i < N_SAMPLES; i++) { rvg_generate(beta_typical_cdf, &p2, cache); }
    double t3 = now_seconds();
    double cached_sec = t3 - t2;
    gsl_rng_free(r2);
    rvg_cache_free(cache);

    return (cached_sec > 0.0) ? (baseline_sec / cached_sec) : -1.0;
}

static void sortD(double *v, int n) {
    for (int i = 1; i < n; i++) {
        double key = v[i];
        int j = i - 1;
        while (j >= 0 && v[j] > key) { v[j + 1] = v[j]; j--; }
        v[j + 1] = key;
    }
}

static double medianD(const double *v, int n) {
    double tmp[REPEATS];
    memcpy(tmp, v, (size_t)n * sizeof(double));
    sortD(tmp, n);
    return (n % 2 == 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}

int main(void) {
    double speedups[REPEATS];

    printf("=== BETA N=%d, %d repeats, 50 new fixed seeds ===\n\n", N_SAMPLES, REPEATS);

    for (int k = 0; k < REPEATS; k++) {
        speedups[k] = run_one(SEEDS[k]);
        printf("seed=%-10lu speedup=%.6f\n", SEEDS[k], speedups[k]);
    }

    double sorted[REPEATS];
    memcpy(sorted, speedups, sizeof sorted);
    sortD(sorted, REPEATS);

    printf("\n=== Sorted speedup values (all %d) ===\n\n", REPEATS);
    for (int i = 0; i < REPEATS; i++) {
        printf("%2d: %.6f\n", i + 1, sorted[i]);
    }

    double median = medianD(speedups, REPEATS);

    double dev[REPEATS];
    for (int i = 0; i < REPEATS; i++) { dev[i] = sorted[i] - median; if (dev[i] < 0) { dev[i] = -dev[i]; } }
    double mad = medianD(dev, REPEATS);

    double threshold = 2.0 * mad;
    double sum = 0.0;
    int kept = 0;
    int excluded_count = 0;
    double excluded_vals[REPEATS];

    printf("\n=== Outlier analysis ===\n\n");
    printf("median = %.6f\n", median);
    printf("MAD (median absolute deviation) = %.6f\n", mad);
    printf("exclusion threshold = 2 * MAD = %.6f (exclude |x - median| > threshold)\n\n", threshold);

    for (int i = 0; i < REPEATS; i++) {
        double d = sorted[i] - median;
        if (d < 0) { d = -d; }
        if (d > threshold) {
            excluded_vals[excluded_count++] = sorted[i];
            printf("EXCLUDED: %.6f (|dev|=%.6f > %.6f)\n", sorted[i], d, threshold);
        } else {
            sum += sorted[i];
            kept++;
        }
    }

    double trimmed_mean = (kept > 0) ? (sum / kept) : 0.0;

    printf("\nkept %d of %d values for trimmed mean\n", kept, REPEATS);
    if (excluded_count == 0) {
        printf("no values excluded\n");
    } else {
        printf("excluded %d value(s):", excluded_count);
        for (int i = 0; i < excluded_count; i++) { printf(" %.6f", excluded_vals[i]); }
        printf("\n");
    }
    printf("\nmedian        = %.6f\n", median);
    printf("trimmed mean  = %.6f\n", trimmed_mean);
    printf("raw min       = %.6f\n", sorted[0]);
    printf("raw max       = %.6f\n", sorted[REPEATS - 1]);

    return 0;
}
