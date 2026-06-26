#include "laplace_z_transform_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

complex_t complex_make(double real, double imag) {
    complex_t z; z.real = real; z.imag = imag; return z;
}
complex_t complex_add(complex_t a, complex_t b) { return complex_make(a.real+b.real, a.imag+b.imag); }
complex_t complex_sub(complex_t a, complex_t b) { return complex_make(a.real-b.real, a.imag-b.imag); }
complex_t complex_mul(complex_t a, complex_t b) { return complex_make(a.real*b.real - a.imag*b.imag, a.real*b.imag + a.imag*b.real); }
complex_t complex_div(complex_t a, complex_t b) {
    double d = b.real*b.real + b.imag*b.imag;
    if (d < 1e-30) return complex_make(1e300, 0.0);
    return complex_make((a.real*b.real + a.imag*b.imag)/d, (a.imag*b.real - a.real*b.imag)/d);
}
double complex_mag(complex_t z) { return sqrt(z.real*z.real + z.imag*z.imag); }
double complex_phase(complex_t z) { return atan2(z.imag, z.real); }
complex_t complex_conj(complex_t z) { return complex_make(z.real, -z.imag); }
complex_t complex_exp(complex_t z) { double m = exp(z.real); return complex_make(m*cos(z.imag), m*sin(z.imag)); }
complex_t complex_sqrt(complex_t z) {
    double r = complex_mag(z), theta = atan2(z.imag, z.real);
    double sr = sqrt(r);
    return complex_make(sr*cos(theta/2.0), sr*sin(theta/2.0));
}
complex_t complex_pow(complex_t z, int n) {
    if (n == 0) return complex_make(1.0, 0.0);
    if (n < 0) { complex_t inv = complex_inv(z); return complex_pow(inv, -n); }
    complex_t result = complex_make(1.0, 0.0), base = z;
    int exp = n;
    while (exp > 0) { if (exp & 1) result = complex_mul(result, base); base = complex_mul(base, base); exp >>= 1; }
    return result;
}
complex_t complex_inv(complex_t z) {
    double d = z.real*z.real + z.imag*z.imag;
    if (d < 1e-30) return complex_make(1e300, 0.0);
    return complex_make(z.real/d, -z.imag/d);
}
complex_t complex_neg(complex_t z) { return complex_make(-z.real, -z.imag); }
int complex_is_real(complex_t z, double eps) { return fabs(z.imag) < eps; }
const char *complex_str(complex_t z) {
    static char buf[128];
    if (fabs(z.imag) < 1e-12) snprintf(buf, sizeof(buf), "%.6g", z.real);
    else if (z.imag > 0) snprintf(buf, sizeof(buf), "%.6g + j%.6g", z.real, z.imag);
    else snprintf(buf, sizeof(buf), "%.6g - j%.6g", z.real, -z.imag);
    return buf;
}

/* Polynomial operations (L3: Polynomial Algebra) */
polynomial_t polynomial_zero(void) { polynomial_t p; p.degree = -1; p.coeff = NULL; return p; }

polynomial_t polynomial_constant(double c) {
    polynomial_t p;
    if (fabs(c) < 1e-30) { p.degree = -1; p.coeff = NULL; }
    else { p.degree = 0; p.coeff = (double*)malloc(sizeof(double)); if (p.coeff) p.coeff[0] = c; }
    return p;
}

polynomial_t polynomial_from_coeffs(const double *coeffs, int n) {
    polynomial_t p;
    if (!coeffs || n <= 0) return polynomial_zero();
    int deg = n - 1;
    while (deg >= 0 && fabs(coeffs[deg]) < 1e-30) deg--;
    if (deg < 0) return polynomial_zero();
    p.degree = deg;
    p.coeff = (double*)malloc((deg + 1) * sizeof(double));
    if (p.coeff) memcpy(p.coeff, coeffs, (deg + 1) * sizeof(double));
    return p;
}

polynomial_t polynomial_copy(const polynomial_t *p) {
    polynomial_t q;
    if (!p || p->degree < 0) return polynomial_zero();
    q.degree = p->degree;
    q.coeff = (double*)malloc((q.degree + 1) * sizeof(double));
    if (q.coeff) memcpy(q.coeff, p->coeff, (q.degree + 1) * sizeof(double));
    return q;
}

void polynomial_free(polynomial_t *p) {
    if (p && p->coeff) { free(p->coeff); p->coeff = NULL; p->degree = -1; }
}

double polynomial_eval(const polynomial_t *p, double x) {
    /* Horner method, O(n), numerically stable */
    if (!p || p->degree < 0) return 0.0;
    double result = 0.0;
    for (int i = p->degree; i >= 0; i--) result = result * x + p->coeff[i];
    return result;
}

complex_t polynomial_eval_complex(const polynomial_t *p, complex_t x) {
    if (!p || p->degree < 0) return complex_make(0.0, 0.0);
    complex_t result = complex_make(0.0, 0.0);
    for (int i = p->degree; i >= 0; i--) {
        complex_t term = complex_mul(result, x);
        result = complex_make(term.real + p->coeff[i], term.imag);
    }
    return result;
}

polynomial_t polynomial_add(const polynomial_t *p, const polynomial_t *q) {
    if (!p || p->degree < 0) return polynomial_copy(q);
    if (!q || q->degree < 0) return polynomial_copy(p);
    int max_deg = (p->degree > q->degree) ? p->degree : q->degree;
    double *coeffs = (double*)calloc(max_deg + 1, sizeof(double));
    if (!coeffs) return polynomial_zero();
    for (int i = 0; i <= p->degree; i++) coeffs[i] += p->coeff[i];
    for (int i = 0; i <= q->degree; i++) coeffs[i] += q->coeff[i];
    polynomial_t result = polynomial_from_coeffs(coeffs, max_deg + 1);
    free(coeffs);
    return result;
}

polynomial_t polynomial_sub(const polynomial_t *p, const polynomial_t *q) {
    if (!p || p->degree < 0) {
        if (!q || q->degree < 0) return polynomial_zero();
        polynomial_t neg = polynomial_copy(q);
        for (int i = 0; i <= neg.degree; i++) neg.coeff[i] = -neg.coeff[i];
        return neg;
    }
    if (!q || q->degree < 0) return polynomial_copy(p);
    int max_deg = (p->degree > q->degree) ? p->degree : q->degree;
    double *coeffs = (double*)calloc(max_deg + 1, sizeof(double));
    if (!coeffs) return polynomial_zero();
    for (int i = 0; i <= p->degree; i++) coeffs[i] += p->coeff[i];
    for (int i = 0; i <= q->degree; i++) coeffs[i] -= q->coeff[i];
    polynomial_t result = polynomial_from_coeffs(coeffs, max_deg + 1);
    free(coeffs);
    return result;
}

polynomial_t polynomial_mul(const polynomial_t *p, const polynomial_t *q) {
    /* Convolution of coefficient arrays, O(n*m) */
    if (!p || p->degree < 0 || !q || q->degree < 0) return polynomial_zero();
    int result_deg = p->degree + q->degree;
    double *coeffs = (double*)calloc(result_deg + 1, sizeof(double));
    if (!coeffs) return polynomial_zero();
    for (int i = 0; i <= p->degree; i++) {
        if (fabs(p->coeff[i]) < 1e-30) continue;
        for (int j = 0; j <= q->degree; j++) coeffs[i + j] += p->coeff[i] * q->coeff[j];
    }
    polynomial_t result = polynomial_from_coeffs(coeffs, result_deg + 1);
    free(coeffs);
    return result;
}

int polynomial_div(const polynomial_t *p, const polynomial_t *q,
                   polynomial_t *quotient, polynomial_t *remainder) {
    /* Polynomial long division, O(n^2) */
    if (!p || !q || q->degree < 0) return -1;
    if (!quotient || !remainder) return -1;
    if (p->degree < 0) { *quotient = polynomial_zero(); *remainder = polynomial_zero(); return 0; }
    if (p->degree < q->degree) { *quotient = polynomial_zero(); *remainder = polynomial_copy(p); return 0; }
    int r_deg = p->degree;
    double *r = (double*)malloc((r_deg + 1) * sizeof(double));
    if (!r) return -1;
    memcpy(r, p->coeff, (r_deg + 1) * sizeof(double));
    int q_deg = q->degree, quot_deg = r_deg - q_deg;
    double *quot = (double*)calloc(quot_deg + 1, sizeof(double));
    if (!quot) { free(r); return -1; }
    double lead_q = q->coeff[q_deg];
    for (int k = quot_deg; k >= 0; k--) {
        double factor = r[q_deg + k] / lead_q;
        quot[k] = factor;
        for (int j = 0; j <= q_deg; j++) r[k + j] -= factor * q->coeff[j];
    }
    while (r_deg >= 0 && fabs(r[r_deg]) < 1e-30) r_deg--;
    *quotient = polynomial_from_coeffs(quot, quot_deg + 1);
    *remainder = polynomial_from_coeffs(r, r_deg + 1);
    free(r); free(quot);
    return 0;
}

polynomial_t polynomial_derivative(const polynomial_t *p) {
    if (!p || p->degree <= 0) return polynomial_zero();
    int new_deg = p->degree - 1;
    double *coeffs = (double*)malloc((new_deg + 1) * sizeof(double));
    if (!coeffs) return polynomial_zero();
    for (int i = 1; i <= p->degree; i++) coeffs[i - 1] = i * p->coeff[i];
    polynomial_t result = polynomial_from_coeffs(coeffs, new_deg + 1);
    free(coeffs);
    return result;
}

polynomial_t polynomial_integral(const polynomial_t *p) {
    if (!p || p->degree < 0) return polynomial_zero();
    int new_deg = p->degree + 1;
    double *coeffs = (double*)calloc(new_deg + 1, sizeof(double));
    if (!coeffs) return polynomial_zero();
    coeffs[0] = 0.0;
    for (int i = 0; i <= p->degree; i++) coeffs[i + 1] = p->coeff[i] / (i + 1.0);
    polynomial_t result = polynomial_from_coeffs(coeffs, new_deg + 1);
    free(coeffs);
    return result;
}

void polynomial_trim(polynomial_t *p) {
    if (!p || p->degree < 0) return;
    while (p->degree >= 0 && fabs(p->coeff[p->degree]) < 1e-30) p->degree--;
    if (p->degree < 0) { free(p->coeff); p->coeff = NULL; }
}

void polynomial_print(const polynomial_t *p, const char *varname) {
    if (!p || p->degree < 0) { printf("0\n"); return; }
    int first = 1;
    for (int i = p->degree; i >= 0; i--) {
        double c = p->coeff[i];
        if (fabs(c) < 1e-15) continue;
        if (!first) printf(c > 0 ? " + " : " - ");
        else if (c < 0) printf("-");
        double abs_c = fabs(c);
        if (i == 0) printf("%.4g", abs_c);
        else if (i == 1) { if (fabs(abs_c - 1.0) > 1e-12) printf("%.4g*", abs_c); printf("%s", varname); }
        else { if (fabs(abs_c - 1.0) > 1e-12) printf("%.4g*", abs_c); printf("%s^%d", varname, i); }
        first = 0;
    }
    if (first) printf("0");
    printf("\n");
}

/* ========================================================================
 * Rational Function Operations (L3: Rational Function Algebra)
 * ======================================================================== */

rational_func_t rational_make(const polynomial_t *num, const polynomial_t *den) {
    rational_func_t rf; rf.num = polynomial_copy(num); rf.den = polynomial_copy(den); return rf;
}
rational_func_t rational_copy(const rational_func_t *rf) {
    rational_func_t c; c.num = polynomial_copy(&rf->num); c.den = polynomial_copy(&rf->den); return c;
}
void rational_free(rational_func_t *rf) { if(rf){polynomial_free(&rf->num);polynomial_free(&rf->den);} }
int rational_is_proper(const rational_func_t *rf) { return rf && rf->den.degree >= 0 && rf->num.degree <= rf->den.degree; }
int rational_is_strictly_proper(const rational_func_t *rf) { return rf && rf->den.degree >= 0 && rf->num.degree < rf->den.degree; }

complex_t rational_eval(const rational_func_t *rf, complex_t x) {
    complex_t nv = polynomial_eval_complex(&rf->num, x);
    complex_t dv = polynomial_eval_complex(&rf->den, x);
    return complex_div(nv, dv);
}

double rational_mag_at(const rational_func_t *rf, double omega, int is_z_domain) {
    complex_t s = is_z_domain ? complex_make(cos(omega), sin(omega)) : complex_make(0.0, omega);
    return complex_mag(rational_eval(rf, s));
}

double rational_phase_at(const rational_func_t *rf, double omega, int is_z_domain) {
    complex_t s = is_z_domain ? complex_make(cos(omega), sin(omega)) : complex_make(0.0, omega);
    return complex_phase(rational_eval(rf, s));
}

rational_func_t rational_add(const rational_func_t *r1, const rational_func_t *r2) {
    polynomial_t n1d2 = polynomial_mul(&r1->num, &r2->den);
    polynomial_t n2d1 = polynomial_mul(&r2->num, &r1->den);
    polynomial_t num = polynomial_add(&n1d2, &n2d1);
    polynomial_t den = polynomial_mul(&r1->den, &r2->den);
    rational_func_t result; result.num = num; result.den = den;
    polynomial_free(&n1d2);polynomial_free(&n2d1);return result;
}

rational_func_t rational_sub(const rational_func_t *r1, const rational_func_t *r2) {
    polynomial_t n1d2 = polynomial_mul(&r1->num, &r2->den);
    polynomial_t n2d1 = polynomial_mul(&r2->num, &r1->den);
    polynomial_t num = polynomial_sub(&n1d2, &n2d1);
    polynomial_t den = polynomial_mul(&r1->den, &r2->den);
    rational_func_t result; result.num = num; result.den = den;
    polynomial_free(&n1d2);polynomial_free(&n2d1);return result;
}

rational_func_t rational_mul(const rational_func_t *r1, const rational_func_t *r2) {
    rational_func_t result;
    result.num = polynomial_mul(&r1->num, &r2->num);
    result.den = polynomial_mul(&r1->den, &r2->den);
    return result;
}

rational_func_t rational_div(const rational_func_t *r1, const rational_func_t *r2) {
    rational_func_t result;
    result.num = polynomial_mul(&r1->num, &r2->den);
    result.den = polynomial_mul(&r1->den, &r2->num);
    return result;
}

rational_func_t rational_reciprocal(const rational_func_t *rf) {
    rational_func_t result;
    result.num = polynomial_copy(&rf->den);
    result.den = polynomial_copy(&rf->num);
    return result;
}

rational_func_t rational_reduce(const rational_func_t *rf) {
    polynomial_t a = polynomial_copy(&rf->num);
    polynomial_t b = polynomial_copy(&rf->den);
    int max_iter = 100;
    while (b.degree >= 0 && max_iter-- > 0) {
        polynomial_t q, r;
        if (polynomial_div(&a, &b, &q, &r) != 0) break;
        polynomial_free(&a); a = b; b = r; polynomial_free(&q);
    }
    rational_func_t result;
    if (a.degree < 0 || a.degree == 0) {
        result.num = polynomial_copy(&rf->num);
        result.den = polynomial_copy(&rf->den);
    } else {
        polynomial_t q1, r1, q2, r2;
        polynomial_div(&rf->num, &a, &q1, &r1);
        polynomial_div(&rf->den, &a, &q2, &r2);
        result.num = q1; result.den = q2;
        polynomial_free(&r1); polynomial_free(&r2);
    }
    polynomial_free(&a); polynomial_free(&b);
    return result;
}

void rational_print(const rational_func_t *rf, const char *varname) {
    printf("N(%s) = ", varname); polynomial_print(&rf->num, varname);
    printf("D(%s) = ", varname); polynomial_print(&rf->den, varname);
}

/* ========================================================================
 * Pole-Zero Operations (L1: Pole-Zero Representation)
 * ======================================================================== */

pole_zero_t pz_create(int n_poles, int n_zeros) {
    pole_zero_t pz; pz.n_poles = n_poles; pz.n_zeros = n_zeros; pz.gain = 1.0;
    pz.poles = n_poles>0 ? (complex_t*)calloc(n_poles, sizeof(complex_t)) : NULL;
    pz.zeros = n_zeros>0 ? (complex_t*)calloc(n_zeros, sizeof(complex_t)) : NULL;
    return pz;
}

void pz_free(pole_zero_t *pz) {
    if(pz){free(pz->poles);free(pz->zeros);pz->poles=NULL;pz->zeros=NULL;pz->n_poles=0;pz->n_zeros=0;}
}

pole_zero_t pz_from_rational(const rational_func_t *rf) {
    pole_zero_t pz;
    int nz = rf->num.degree, np = rf->den.degree;
    pz.n_zeros = nz; pz.n_poles = np; pz.gain = 1.0;
    if (nz >= 0 && np >= 0) pz.gain = rf->num.coeff[nz] / rf->den.coeff[np];
    else if (nz >= 0) pz.gain = rf->num.coeff[nz];
    else if (np >= 0) pz.gain = 1.0 / rf->den.coeff[np];
    pz.zeros = nz>0 ? (complex_t*)calloc(nz, sizeof(complex_t)) : NULL;
    pz.poles = np>0 ? (complex_t*)calloc(np, sizeof(complex_t)) : NULL;

    for (int i = 0; i < nz; i++) {
        double angle = 2.0*M_PI*i/(nz>0?nz:1);
        complex_t x = complex_make(cos(angle), sin(angle));
        for (int iter = 0; iter < 50; iter++) {
            complex_t fx = polynomial_eval_complex(&rf->num, x);
            polynomial_t deriv = polynomial_derivative(&rf->num);
            complex_t fpx = polynomial_eval_complex(&deriv, x);
            polynomial_free(&deriv);
            if (complex_mag(fpx) < 1e-15) break;
            complex_t dx = complex_div(fx, fpx);
            x = complex_sub(x, dx);
            if (complex_mag(dx) < 1e-12) break;
        }
        pz.zeros[i] = x;
    }
    for (int i = 0; i < np; i++) {
        double angle = 2.0*M_PI*i/(np>0?np:1);
        complex_t x = complex_make(cos(angle), sin(angle));
        for (int iter = 0; iter < 50; iter++) {
            complex_t fx = polynomial_eval_complex(&rf->den, x);
            polynomial_t deriv = polynomial_derivative(&rf->den);
            complex_t fpx = polynomial_eval_complex(&deriv, x);
            polynomial_free(&deriv);
            if (complex_mag(fpx) < 1e-15) break;
            complex_t dx = complex_div(fx, fpx);
            x = complex_sub(x, dx);
            if (complex_mag(dx) < 1e-12) break;
        }
        pz.poles[i] = x;
    }
    return pz;
}

rational_func_t pz_to_rational(const pole_zero_t *pz) {
    polynomial_t num = polynomial_constant(1.0);
    polynomial_t den = polynomial_constant(1.0);
    for (int i = 0; i < pz->n_zeros; i++) {
        double c[2] = {-pz->zeros[i].real, 1.0};
        polynomial_t f = polynomial_from_coeffs(c,2);
        polynomial_t nn = polynomial_mul(&num,&f);
        polynomial_free(&num);polynomial_free(&f);num=nn;
    }
    for (int i = 0; i < pz->n_poles; i++) {
        double c[2] = {-pz->poles[i].real, 1.0};
        polynomial_t f = polynomial_from_coeffs(c,2);
        polynomial_t nd = polynomial_mul(&den,&f);
        polynomial_free(&den);polynomial_free(&f);den=nd;
    }
    for (int i = 0; i <= num.degree; i++) num.coeff[i] *= pz->gain;
    rational_func_t r; r.num=num; r.den=den; return r;
}

complex_t pz_eval(const pole_zero_t *pz, complex_t x) {
    complex_t num = complex_make(pz->gain, 0.0);
    complex_t den = complex_make(1.0, 0.0);
    for (int i = 0; i < pz->n_zeros; i++) { complex_t d=complex_sub(x,pz->zeros[i]); num=complex_mul(num,d); }
    for (int i = 0; i < pz->n_poles; i++) { complex_t d=complex_sub(x,pz->poles[i]); den=complex_mul(den,d); }
    return complex_div(num, den);
}

double pz_mag_at(const pole_zero_t *pz, double omega, int is_z_domain) {
    complex_t x = is_z_domain ? complex_make(cos(omega), sin(omega)) : complex_make(0.0, omega);
    return complex_mag(pz_eval(pz, x));
}

double pz_phase_at(const pole_zero_t *pz, double omega, int is_z_domain) {
    complex_t x = is_z_domain ? complex_make(cos(omega), sin(omega)) : complex_make(0.0, omega);
    return complex_phase(pz_eval(pz, x));
}

void pz_print(const pole_zero_t *pz) {
    printf("Pole-Zero Plot:\nGain K = %.6g\nZeros (%d):\n", pz->gain, pz->n_zeros);
    for (int i=0;i<pz->n_zeros;i++) printf("  z%d = %s\n",i+1,complex_str(pz->zeros[i]));
    printf("Poles (%d):\n",pz->n_poles);
    for (int i=0;i<pz->n_poles;i++) printf("  p%d = %s\n",i+1,complex_str(pz->poles[i]));
}

/* ========================================================================
 * Partial Fraction Data Structures (L1)
 * ======================================================================== */

partial_fraction_t pf_create(int n_terms) {
    partial_fraction_t pf; pf.n_terms = n_terms;
    pf.terms = n_terms>0 ? (pf_term_t*)calloc(n_terms, sizeof(pf_term_t)) : NULL;
    pf.direct = polynomial_zero();
    return pf;
}

void pf_free(partial_fraction_t *pf) {
    if(pf){free(pf->terms);polynomial_free(&pf->direct);pf->terms=NULL;pf->n_terms=0;}
}

void pf_print(const partial_fraction_t *pf, int is_z_domain) {
    if(pf->direct.degree >= 0){printf("Direct term: ");polynomial_print(&pf->direct,"t");}
    for(int i=0;i<pf->n_terms;i++){
        pf_term_t *t=&pf->terms[i];
        if(is_z_domain)
            printf("  %s / (1 - %s*z^{-1})^%d\n",complex_str(t->residue),complex_str(t->pole),t->multiplicity);
        else
            printf("  %s / (s - %s)^%d\n",complex_str(t->residue),complex_str(t->pole),t->multiplicity);
    }
}
