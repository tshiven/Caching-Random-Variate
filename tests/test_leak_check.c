/*
  Name:     test_leak_check.c
  Purpose:  Exercise rvg_generate across several distributions (head-only and
            head+tail) -- with NO manual create/free of any kind, since the
            public API no longer exposes a cache object at all -- so `leaks`
            has real allocation activity to inspect, and to prove the
            automatic atexit-registered cleanup actually runs and frees
            everything when the process exits normally.
  Note:     New file; does not modify any existing file. Compile with debug
            symbols:
              gcc -g -O0 $(gsl-config --cflags) -I. -Icache -I<gmp include> \
                  -L<gmp lib> tests/test_leak_check.c cache/rvg_cache.c \
                  librvg.a -o tests/test_leak_check.out \
                  $(gsl-config --libs-without-cblas) -lgmp
            Then run under: leaks --atExit -- ./tests/test_leak_check.out
*/

#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

#include "generate.h"
#include "flip.h"
#include "rvg_cache.h"

MAKE_CDF_P(beta_cdf, gsl_cdf_beta_P, 5, 5)
MAKE_CDF_P(gamma_cdf, gsl_cdf_gamma_P, 0.5, 1)
MAKE_CDF_UINT_P(poisson_cdf, gsl_cdf_poisson_P, 71)
MAKE_CDF_UINT_P(binomial_cdf, gsl_cdf_binomial_P, 0.2, 100)
MAKE_CDF_UINT_P(hypergeometric_cdf, gsl_cdf_hypergeometric_P, 5, 20, 7)

typedef struct {
    const char *name;
    rvg_dist_t  dist;
    cdf32_t     cdf;
} entry_t;

static entry_t ENTRIES[] = {
    {"BETA (head-only)",           RVG_DIST_BETA,           beta_cdf},
    {"GAMMA (head-only)",          RVG_DIST_GAMMA,          gamma_cdf},
    {"POISSON (head+tail)",        RVG_DIST_POISSON,        poisson_cdf},
    {"BINOMIAL (head+tail)",       RVG_DIST_BINOMIAL,       binomial_cdf},
    {"HYPERGEOMETRIC (head+tail)", RVG_DIST_HYPERGEOMETRIC, hypergeometric_cdf},
};
#define NUM_ENTRIES ((int)(sizeof(ENTRIES) / sizeof(ENTRIES[0])))

int main(void) {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 42);
    struct flip_state prng = make_flip_state(rng);

    /* Just call rvg_generate, repeatedly, for several distributions -- no
       create, no free, nothing to manage. The first call for each distinct
       cdf silently builds its internal cache; nothing here ever tears one
       down manually. If this leaks, it will show up in `leaks --atExit`. */
    for (int round = 0; round < 2; round++) {
        for (int i = 0; i < NUM_ENTRIES; i++) {
            entry_t *e = &ENTRIES[i];
            for (int k = 0; k < 2000; k++) {
                rvg_generate(e->cdf, &prng, e->dist, 0);
            }
            printf("round=%d %-28s: 2000 samples drawn\n", round, e->name);
        }
    }

    gsl_rng_free(rng);
    printf("done -- process exiting normally now; atexit cleanup should run here\n");
    return 0;
}
