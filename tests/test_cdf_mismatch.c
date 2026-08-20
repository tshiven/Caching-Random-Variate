/*
  Name:     test_cdf_mismatch.c
  Purpose:  Originally tested that reusing one manually-created rvg_cache_t
            across two different CDFs triggered a deliberate abort(). That
            scenario no longer exists: the public API is now a single
            rvg_generate(cdf, prng, dist, force_unsafe) with no cache handle
            for a caller to misuse -- caching is entirely internal, keyed on
            the exact `cdf` function pointer automatically.

            This file now proves the *positive* version of that same
            guarantee: using several different CDFs -- including two
            different parameterizations of the SAME distribution family,
            interleaved in the same loop -- never crashes and never mixes up
            results between them, because each distinct `cdf` function
            automatically gets its own independent internal cache.
  Note:     New file; does not modify any existing file. Compile manually:
              gcc -O0 -g $(gsl-config --cflags) -I. -Icache -I<gmp include> \
                  -L<gmp lib> tests/test_cdf_mismatch.c cache/rvg_cache.c \
                  librvg.a -o tests/test_cdf_mismatch.out \
                  $(gsl-config --libs-without-cblas) -lgmp
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

#include "generate.h"
#include "flip.h"
#include "rvg_cache.h"

/* Two DIFFERENT Gamma parameterizations -- same rvg_dist_t (RVG_DIST_GAMMA),
   but different cdf32_t function pointers. Under the old manual-cache API,
   sharing one rvg_cache_t between these would deliberately abort(). Under
   the new automatic API, there is no shared cache to misuse: each gets its
   own automatically. */
MAKE_CDF_P(gamma_cdf_a, gsl_cdf_gamma_P, 2, 3)     /* Gamma(shape=2, scale=3) */
MAKE_CDF_P(gamma_cdf_b, gsl_cdf_gamma_P, 0.5, 1)   /* Gamma(shape=0.5, scale=1) */
MAKE_CDF_UINT_P(poisson_cdf, gsl_cdf_poisson_P, 71)

static int double_bits_eq(double a, double b) {
    uint64_t ba, bb;
    memcpy(&ba, &a, sizeof ba);
    memcpy(&bb, &b, sizeof bb);
    return ba == bb;
}

int main(void) {
    printf("=== Interleaved-CDF test: multiple different CDFs, no crash, no cross-contamination ===\n\n");

    long mismatches = 0;
    const int ROUNDS = 3000;

    /* One PRNG stream shared across all three CDFs, exactly as a real
       program interleaving several distributions would do. */
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 42);
    struct flip_state prng = make_flip_state(rng);

    /* A second, independent stream to compute the ground truth with
       generate_opt, seeded identically so the two streams stay comparable. */
    gsl_rng *rng_ref = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng_ref, 42);
    struct flip_state prng_ref = make_flip_state(rng_ref);

    for (int i = 0; i < ROUNDS; i++) {
        double ref_a = generate_opt(gamma_cdf_a, &prng_ref);
        double got_a = rvg_generate(gamma_cdf_a, &prng, RVG_DIST_GAMMA, 0);
        if (!double_bits_eq(ref_a, got_a)) {
            mismatches++;
            printf("  MISMATCH round=%d gamma_a: expected=%.17g got=%.17g\n", i, ref_a, got_a);
        }

        double ref_b = generate_opt(gamma_cdf_b, &prng_ref);
        double got_b = rvg_generate(gamma_cdf_b, &prng, RVG_DIST_GAMMA, 0);
        if (!double_bits_eq(ref_b, got_b)) {
            mismatches++;
            printf("  MISMATCH round=%d gamma_b: expected=%.17g got=%.17g\n", i, ref_b, got_b);
        }

        double ref_p = generate_opt(poisson_cdf, &prng_ref);
        double got_p = rvg_generate(poisson_cdf, &prng, RVG_DIST_POISSON, 0);
        if (!double_bits_eq(ref_p, got_p)) {
            mismatches++;
            printf("  MISMATCH round=%d poisson: expected=%.17g got=%.17g\n", i, ref_p, got_p);
        }
    }

    gsl_rng_free(rng);
    gsl_rng_free(rng_ref);

    printf("%d rounds x 3 interleaved CDFs (two different Gammas + one Poisson) = %d calls\n",
           ROUNDS, ROUNDS * 3);
    printf("mismatches=%ld\n", mismatches);

    if (mismatches == 0) {
        printf("RESULT: PASS -- no crash, and every interleaved CDF stayed independently correct.\n");
        return 0;
    } else {
        printf("RESULT: FAIL -- see MISMATCH lines above.\n");
        return 1;
    }
}
