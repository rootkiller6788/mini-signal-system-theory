/*
 * laplace_z_transform_core.h - Core Data Structures for Laplace and Z Transforms
 *
 * Provides foundational complex number, polynomial, and rational function
 * types used throughout the Laplace/Z transform module.
 *
 * Knowledge Coverage:
 *   L1: complex numbers (s-domain, z-domain variables)
 *   L3: polynomial algebra, rational function representation
 *
 * Reference: Oppenheim and Willsky, Signals and Systems (1997), Ch.9
 *            Proakis and Manolakis, Digital Signal Processing (2007), Ch.3
 */

#ifndef LAPLACE_Z_TRANSFORM_CORE_H
#define LAPLACE_Z_TRANSFORM_CORE_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Complex Number Type - foundation for s-plane and z-plane analysis
 * ========================================================================== */

typedef struct {
    double real;
    double imag;
} complex_t;

complex_t complex_make(double real, double imag);
complex_t complex_add(complex_t a, complex_t b);
complex_t complex_sub(complex_t a, complex_t b);
complex_t complex_mul(complex_t a, complex_t b);
complex_t complex_div(complex_t a, complex_t b);
double complex_mag(complex_t z);
double complex_phase(complex_t z);
complex_t complex_conj(complex_t z);
complex_t complex_exp(complex_t z);
complex_t complex_sqrt(complex_t z);
complex_t complex_pow(complex_t z, int n);
complex_t complex_inv(complex_t z);
complex_t complex_neg(complex_t z);
int complex_is_real(complex_t z, double eps);
const char *complex_str(complex_t z);

/* ==========================================================================
 * L3: Polynomial Type
 * ========================================================================== */

typedef struct {
    int degree;
    double *coeff;
} polynomial_t;

polynomial_t polynomial_zero(void);
polynomial_t polynomial_constant(double c);
polynomial_t polynomial_from_coeffs(const double *coeffs, int n);
polynomial_t polynomial_copy(const polynomial_t *p);
void polynomial_free(polynomial_t *p);
double polynomial_eval(const polynomial_t *p, double x);
complex_t polynomial_eval_complex(const polynomial_t *p, complex_t x);
polynomial_t polynomial_add(const polynomial_t *p, const polynomial_t *q);
polynomial_t polynomial_sub(const polynomial_t *p, const polynomial_t *q);
polynomial_t polynomial_mul(const polynomial_t *p, const polynomial_t *q);
int polynomial_div(const polynomial_t *p, const polynomial_t *q,
                   polynomial_t *quotient, polynomial_t *remainder);
polynomial_t polynomial_derivative(const polynomial_t *p);
polynomial_t polynomial_integral(const polynomial_t *p);
void polynomial_trim(polynomial_t *p);
void polynomial_print(const polynomial_t *p, const char *varname);

/* ==========================================================================
 * L3: Rational Function
 * ========================================================================== */

typedef struct {
    polynomial_t num;
    polynomial_t den;
} rational_func_t;

rational_func_t rational_make(const polynomial_t *num, const polynomial_t *den);
rational_func_t rational_copy(const rational_func_t *rf);
void rational_free(rational_func_t *rf);
int rational_is_proper(const rational_func_t *rf);
int rational_is_strictly_proper(const rational_func_t *rf);
complex_t rational_eval(const rational_func_t *rf, complex_t x);
double rational_mag_at(const rational_func_t *rf, double omega, int is_z_domain);
double rational_phase_at(const rational_func_t *rf, double omega, int is_z_domain);
rational_func_t rational_add(const rational_func_t *r1, const rational_func_t *r2);
rational_func_t rational_sub(const rational_func_t *r1, const rational_func_t *r2);
rational_func_t rational_mul(const rational_func_t *r1, const rational_func_t *r2);
rational_func_t rational_div(const rational_func_t *r1, const rational_func_t *r2);
rational_func_t rational_reciprocal(const rational_func_t *rf);
rational_func_t rational_reduce(const rational_func_t *rf);
void rational_print(const rational_func_t *rf, const char *varname);

/* ==========================================================================
 * L1: Pole-Zero Representation
 * ========================================================================== */

typedef struct {
    int n_poles;
    int n_zeros;
    complex_t *poles;
    complex_t *zeros;
    double gain;
} pole_zero_t;

pole_zero_t pz_create(int n_poles, int n_zeros);
void pz_free(pole_zero_t *pz);
pole_zero_t pz_from_rational(const rational_func_t *rf);
rational_func_t pz_to_rational(const pole_zero_t *pz);
complex_t pz_eval(const pole_zero_t *pz, complex_t x);
double pz_mag_at(const pole_zero_t *pz, double omega, int is_z_domain);
double pz_phase_at(const pole_zero_t *pz, double omega, int is_z_domain);
void pz_print(const pole_zero_t *pz);

/* ==========================================================================
 * L1: Partial Fraction
 * ========================================================================== */

typedef struct {
    complex_t residue;
    complex_t pole;
    int multiplicity;
} pf_term_t;

typedef struct {
    int n_terms;
    pf_term_t *terms;
    polynomial_t direct;
} partial_fraction_t;

partial_fraction_t pf_create(int n_terms);
void pf_free(partial_fraction_t *pf);
void pf_print(const partial_fraction_t *pf, int is_z_domain);

#ifdef __cplusplus
}
#endif

#endif /* LAPLACE_Z_TRANSFORM_CORE_H */
