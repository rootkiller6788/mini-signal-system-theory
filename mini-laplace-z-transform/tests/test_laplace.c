/* test_laplace.c - Tests for Laplace transform operations */
#include "laplace_z_transform_core.h"
#include "laplace_transform.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-8

int main(void) {
    printf("test_laplace: Starting...\n");

    /* Test Laplace transform pairs */
    rational_func_t H_step = laplace_step();
    assert(rational_is_strictly_proper(&H_step));

    double d0 = polynomial_eval(&H_step.den, 0.0);
    assert(fabs(d0) < EPS);

    /* exp(-2t): 1/(s+2), DC gain = 1/2 */
    rational_func_t H_exp = laplace_exp_decay(2.0);
    double dc = rational_mag_at(&H_exp, 0.0, 0);
    assert(fabs(dc - 0.5) < EPS);

    /* cos(t): s/(s^2+1) */
    rational_func_t H_cos = laplace_cos(1.0);
    complex_t val0 = rational_eval(&H_cos, complex_make(0.0, 0.0));
    assert(fabs(val0.real) < EPS);

    /* sin(t): 1/(s^2+1), DC gain = 1 */
    rational_func_t H_sin = laplace_sin(1.0);
    double dcs = rational_mag_at(&H_sin, 0.0, 0);
    assert(fabs(dcs - 1.0) < EPS);

    /* Second order: wn=2, zeta=0.5 */
    rational_func_t H2 = laplace_second_order(1.0, 2.0, 0.5);
    double dc2 = rational_mag_at(&H2, 0.0, 0);
    assert(fabs(dc2 - 1.0) < EPS);

    /* Convolution: step*step = ramp (both are 1/s^2) */
    rational_func_t H_r1 = laplace_step();
    rational_func_t H_conv = laplace_convolution(&H_r1, &H_r1);
    rational_func_t H_ramp = laplace_ramp();
    assert(H_conv.den.degree == H_ramp.den.degree);

    /* Derivative property at non-zero frequency */
    rational_func_t H_deriv = laplace_derivative(&H_step);
    double dcd = rational_mag_at(&H_deriv, 0.01, 0);
    assert(fabs(dcd - 1.0) < 0.02);

    /* Frequency shifting */
    rational_func_t H_shift = laplace_freq_shift(&H_cos, 3.0);
    complex_t vs = rational_eval(&H_shift, complex_make(0.0, 0.0));
    assert(vs.real > 0);

    /* Initial value theorem: for H(s)=1/(s+2), h(0+)=1 */
    double iv = laplace_initial_value(&H_exp);
    assert(fabs(iv - 1.0) < EPS);

    /* Final value theorem: for h(t)=e^{-2t}, h(inf)=0 */
    double fv = laplace_final_value(&H_exp);
    assert(fabs(fv) < EPS || isnan(fv));

    /* Stability: 1/(s+2) stable */
    assert(laplace_is_stable(&H_exp) == 1);

    /* Fourier at DC */
    complex_t ft = laplace_to_fourier(&H_exp, 0.0);
    assert(fabs(ft.real - 0.5) < EPS);

    /* Parseval energy (non-negative) */
    double energy = laplace_parseval_energy(&H_exp, 10.0, 100);
    assert(energy >= 0.0);

    rational_free(&H_step); rational_free(&H_exp); rational_free(&H_cos);
    rational_free(&H_sin); rational_free(&H2);
    rational_free(&H_r1); rational_free(&H_conv); rational_free(&H_ramp);
    rational_free(&H_deriv); rational_free(&H_shift);

    printf("test_laplace: All tests passed!\n");
    return 0;
}
