/* test_z.c - Tests for Z-transform operations */
#include "laplace_z_transform_core.h"
#include "z_transform.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-8

int main(void) {
    printf("test_z: Starting...\n");

    rational_func_t X_imp = zt_impulse();
    assert(rational_is_proper(&X_imp));

    rational_func_t X_step = zt_step();
    double d0 = polynomial_eval(&X_step.den, 1.0);
    assert(fabs(d0) < EPS);

    /* Geometric a=0.5: 1/(1-0.5z^{-1}), DC = 2 */
    rational_func_t X_geo = zt_geometric(0.5);
    double dc = rational_mag_at(&X_geo, 0.0, 1);
    assert(fabs(dc - 2.0) < EPS);

    /* cos(pi/4*n) */
    rational_func_t X_cos = zt_cos(M_PI/4.0);
    double dc_cos = rational_mag_at(&X_cos, 0.0, 1);
    assert(dc_cos > 0);

    /* Stability */
    assert(zt_is_stable(&X_geo) == 1);
    rational_func_t X_unstable = zt_geometric(2.0);
    assert(zt_is_stable(&X_unstable) == 0);
    assert(zt_roc_includes_unit_circle(&X_geo) == 1);
    assert(zt_roc_includes_unit_circle(&X_unstable) == 0);

    /* Initial value: x[0] = 1 */
    double iv = zt_initial_value(&X_geo);
    assert(fabs(iv - 1.0) < EPS);

    /* Convolution: deg1 * deg1 = deg2 */
    rational_func_t X_conv = zt_convolution(&X_geo, &X_geo);
    assert(X_conv.den.degree == 2);

    /* Time shift right by 2 */
    rational_func_t X_shift = zt_time_shift_right(&X_geo, 2);
    assert(X_shift.num.degree == X_geo.num.degree + 2);

    rational_free(&X_imp); rational_free(&X_step); rational_free(&X_geo);
    rational_free(&X_cos); rational_free(&X_unstable);
    rational_free(&X_conv); rational_free(&X_shift);

    printf("test_z: All tests passed!\n");
    return 0;
}
