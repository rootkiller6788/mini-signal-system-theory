/* test_stability.c - Tests for stability analysis */
#include "laplace_z_transform_core.h"
#include "stability_analysis.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-8

int main(void) {
    printf("test_stability: Starting...\n");
    fflush(stdout);

    /* Routh-Hurwitz: s^2 + 2s + 1 = (s+1)^2, stable */
    double poly1[] = {1.0, 2.0, 1.0};
    assert(stability_routh_necessary_condition(poly1, 2) == 1);
    assert(stability_routh_is_stable(poly1, 2) == 1);
    int rhp1 = stability_routh_rhp_roots(poly1, 2);
    assert(rhp1 == 0);
    printf("  Routh-Hurwitz basic: PASS\n");
    fflush(stdout);

    /* Routh-Hurwitz: s^2 - s + 1, unstable */
    double poly2[] = {1.0, -1.0, 1.0};
    assert(stability_routh_necessary_condition(poly2, 2) == 0);
    printf("  Routh-Hurwitz unstable: PASS\n");
    fflush(stdout);

    /* Jury: z^2 - 0.5z + 0.06, stable */
    double jpoly1[] = {0.06, -0.5, 1.0};
    assert(stability_jury_is_stable(jpoly1, 2) == 1);
    printf("  Jury stable: PASS\n");
    fflush(stdout);

    /* Jury: unstable */
    double jpoly2[] = {2.0, 0.0, 1.0};
    assert(stability_jury_is_stable(jpoly2, 2) == 0);
    printf("  Jury unstable: PASS\n");
    fflush(stdout);

    /* Pole-based stability (continuous) */
    complex_t poles_stable[] = {
        complex_make(-1.0, 0.0), complex_make(-0.5, 1.0), complex_make(-0.5, -1.0)
    };
    assert(stability_from_poles(poles_stable, 3, 0) == 1);
    printf("  Pole stability: PASS\n");
    fflush(stdout);

    /* Discrete-time pole stability */
    complex_t dpoles_stable[] = { complex_make(0.5, 0.0), complex_make(0.0, 0.8) };
    assert(stability_from_poles(dpoles_stable, 2, 1) == 1);
    printf("  Discrete pole stability: PASS\n");
    fflush(stdout);

    /* Stability margin */
    double margin_stable = stability_margin(poles_stable, 3, 0);
    assert(margin_stable > 0);
    printf("  Stability margin: PASS\n");
    fflush(stdout);

    printf("test_stability: All tests passed!\n");
    return 0;
}
