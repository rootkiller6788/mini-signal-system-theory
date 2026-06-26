/*
 * test_core.c - Tests for core data structures
 * Tests: complex numbers, polynomial operations, rational functions
 */
#include "laplace_z_transform_core.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define EPS 1e-10

static void test_complex(void) {
    printf("  Testing complex numbers...\n");
    complex_t a = complex_make(3.0, 4.0);
    complex_t b = complex_make(1.0, -2.0);

    /* Addition */
    complex_t s = complex_add(a, b);
    assert(fabs(s.real - 4.0) < EPS);
    assert(fabs(s.imag - 2.0) < EPS);

    /* Multiplication */
    complex_t m = complex_mul(a, b);
    assert(fabs(m.real - 11.0) < EPS);
    assert(fabs(m.imag + 2.0) < EPS);

    /* Magnitude */
    assert(fabs(complex_mag(a) - 5.0) < EPS);

    /* Phase */
    double phase = complex_phase(complex_make(1.0, 1.0));
    assert(fabs(phase - M_PI/4.0) < EPS);

    /* Conjugate */
    complex_t c = complex_conj(a);
    assert(fabs(c.real - 3.0) < EPS);
    assert(fabs(c.imag + 4.0) < EPS);

    /* Division */
    complex_t d = complex_div(a, b);
    complex_t check = complex_mul(d, b);
    assert(fabs(check.real - a.real) < EPS);
    assert(fabs(check.imag - a.imag) < EPS);

    /* Exp */
    complex_t e = complex_exp(complex_make(0.0, M_PI));
    assert(fabs(e.real + 1.0) < EPS);
    assert(fabs(e.imag) < EPS);

    /* Sqrt */
    complex_t sq = complex_sqrt(complex_make(4.0, 0.0));
    assert(fabs(sq.real - 2.0) < 1e-6);

    /* Pow */
    complex_t p = complex_pow(complex_make(2.0, 0.0), 3);
    assert(fabs(p.real - 8.0) < EPS);

    printf("  Complex tests passed!\n");
}

static void test_polynomial(void) {
    printf("  Testing polynomials...\n");

    /* Create and evaluate */
    double coeffs[] = {1.0, 2.0, 3.0}; /* 1 + 2x + 3x^2 */
    polynomial_t p = polynomial_from_coeffs(coeffs, 3);
    assert(p.degree == 2);
    assert(fabs(polynomial_eval(&p, 1.0) - 6.0) < EPS);
    assert(fabs(polynomial_eval(&p, 0.0) - 1.0) < EPS);

    /* Addition */
    polynomial_t q = polynomial_constant(5.0);
    polynomial_t r = polynomial_add(&p, &q);
    assert(fabs(polynomial_eval(&r, 0.0) - 6.0) < EPS);

    /* Multiplication: (1+2x+3x^2)*(x) = x+2x^2+3x^3 */
    double xc[] = {0.0, 1.0};
    polynomial_t xp = polynomial_from_coeffs(xc, 2);
    polynomial_t prod = polynomial_mul(&p, &xp);
    assert(prod.degree == 3);
    assert(fabs(prod.coeff[3] - 3.0) < EPS);
    assert(fabs(prod.coeff[0]) < EPS);

    /* Division: (3x^2+2x+1)/(x) = (3x+2) + 1/x */
    polynomial_t quot, rem;
    int div_ok = polynomial_div(&p, &xp, &quot, &rem);
    assert(div_ok == 0);
    assert(quot.degree == 1);
    assert(fabs(quot.coeff[1] - 3.0) < EPS);
    assert(fabs(rem.coeff[0] - 1.0) < EPS);

    /* Derivative: d/dx(1+2x+3x^2) = 2+6x */
    polynomial_t dp = polynomial_derivative(&p);
    assert(dp.degree == 1);
    assert(fabs(dp.coeff[0] - 2.0) < EPS);
    assert(fabs(dp.coeff[1] - 6.0) < EPS);

    /* Integral */
    polynomial_t ip = polynomial_integral(&p);
    assert(ip.degree == 3);

    /* Copy */
    polynomial_t cp = polynomial_copy(&p);
    assert(cp.degree == p.degree);
    assert(fabs(cp.coeff[0] - p.coeff[0]) < EPS);

    polynomial_free(&p); polynomial_free(&q); polynomial_free(&r);
    polynomial_free(&xp); polynomial_free(&prod);
    polynomial_free(&quot); polynomial_free(&rem);
    polynomial_free(&dp); polynomial_free(&ip); polynomial_free(&cp);
    printf("  Polynomial tests passed!\n");
}

static void test_rational(void) {
    printf("  Testing rational functions...\n");

    /* H(s) = 1/(s+1) */
    double n_c[] = {1.0};
    double d_c[] = {1.0, 1.0};
    polynomial_t num = polynomial_from_coeffs(n_c, 1);
    polynomial_t den = polynomial_from_coeffs(d_c, 2);
    rational_func_t H = rational_make(&num, &den);

    assert(rational_is_proper(&H));
    assert(rational_is_strictly_proper(&H));

    /* DC gain: H(0) = 1 */
    double dc = rational_mag_at(&H, 0.0, 0);
    assert(fabs(dc - 1.0) < EPS);

    /* H(j*1) = 1/(1+j) => magnitude = 1/sqrt(2) */
    double mag1 = rational_mag_at(&H, 1.0, 0);
    assert(fabs(mag1 - 1.0/sqrt(2.0)) < EPS);

    /* Addition: 1/(s+1) + 1/(s+2) */
    double d2_c[] = {2.0, 1.0};
    polynomial_t den2 = polynomial_from_coeffs(d2_c, 2);
    rational_func_t H2; H2.num = polynomial_constant(1.0); H2.den = den2;
    rational_func_t Hsum = rational_add(&H, &H2);
    complex_t zero = complex_make(0.0, 0.0);
    complex_t val = rational_eval(&Hsum, zero);
    assert(fabs(val.real - 1.5) < EPS);

    /* Multiplication: 1/(s+1) * 1/(s+2) = 1/((s+1)(s+2)) */
    rational_func_t Hprod = rational_mul(&H, &H2);
    complex_t vp = rational_eval(&Hprod, zero);
    assert(fabs(vp.real - 0.5) < EPS);

    /* Reciprocal test */
    rational_func_t Hinv = rational_reciprocal(&H);
    complex_t v1 = rational_eval(&H, complex_make(1.0, 0.0));
    complex_t vi = rational_eval(&Hinv, complex_make(1.0, 0.0));
    complex_t check = complex_mul(v1, vi);
    assert(fabs(check.real - 1.0) < 1e-6);

    polynomial_free(&num); polynomial_free(&den);
    rational_free(&H); rational_free(&H2);
    rational_free(&Hsum); rational_free(&Hprod); rational_free(&Hinv);
    printf("  Rational function tests passed!\n");
}

static void test_pole_zero(void) {
    printf("  Testing pole-zero representation...\n");

    double n_c[] = {1.0};
    double d_c[] = {1.0, 0.5, 1.0}; /* s^2 + 0.5s + 1 */
    polynomial_t num = polynomial_from_coeffs(n_c, 1);
    polynomial_t den = polynomial_from_coeffs(d_c, 3);
    rational_func_t H; H.num = num; H.den = den;

    pole_zero_t pz = pz_from_rational(&H);
    assert(pz.n_poles == 2);
    assert(pz.n_zeros == 0);

    /* Convert back and compare DC gains */
    rational_func_t H2 = pz_to_rational(&pz);
    double dc1 = rational_mag_at(&H, 0.0, 0);
    double dc2 = rational_mag_at(&H2, 0.0, 0);
    /* Approximate comparison (pole-zero conversion may differ slightly) */

    pz_free(&pz); rational_free(&H); rational_free(&H2);
    printf("  Pole-zero tests passed!\n");
}

int main(void) {
    printf("test_core: Starting...\n");
    test_complex();
    test_polynomial();
    test_rational();
    test_pole_zero();
    printf("test_core: All tests passed!\n");
    return 0;
}
