/* partial_fraction.c - Partial Fraction Expansion */
#include "laplace_z_transform_core.h"
#include "partial_fraction.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static complex_t find_one_root_newton(const polynomial_t *p, complex_t x0, int max_iter) {
    complex_t x = x0;
    for (int iter = 0; iter < max_iter; iter++) {
        complex_t fx = polynomial_eval_complex(p, x);
        polynomial_t dp = polynomial_derivative(p);
        complex_t fpx = polynomial_eval_complex(&dp, x);
        polynomial_free(&dp);
        if (complex_mag(fpx) < 1e-15) break;
        complex_t dx = complex_div(fx, fpx);
        x = complex_sub(x, dx);
        if (complex_mag(dx) < 1e-12) break;
    }
    return x;
}

static polynomial_t deflate_polynomial(const polynomial_t *p, complex_t root) {
    if (p->degree <= 0) return polynomial_copy(p);
    int n = p->degree;
    double *qcoeffs = (double*)calloc(n, sizeof(double));
    if (fabs(root.imag) < 1e-12) {
        qcoeffs[n - 1] = p->coeff[n];
        for (int k = n - 1; k >= 1; k--) qcoeffs[k - 1] = p->coeff[k] + root.real * qcoeffs[k];
    } else {
        double b1 = -2.0 * root.real, b0 = root.real * root.real + root.imag * root.imag;
        if (n >= 2) {
            qcoeffs[n - 1] = p->coeff[n];
            qcoeffs[n - 2] = p->coeff[n - 1] - b1 * qcoeffs[n - 1];
            for (int k = n - 2; k >= 1; k--)
                qcoeffs[k - 1] = p->coeff[k] - b1 * qcoeffs[k] - b0 * qcoeffs[k + 1];
        }
    }
    polynomial_t result = polynomial_from_coeffs(qcoeffs, n);
    free(qcoeffs);
    return result;
}

int pf_find_poles(const polynomial_t *den, complex_t *poles, int *mult, int *n_distinct) {
    if (!den || den->degree < 0 || !poles || !mult || !n_distinct) return -1;
    int deg = den->degree;
    if (deg == 0) { *n_distinct = 0; return 0; }
    polynomial_t working = polynomial_copy(den);
    int found = 0;
    for (int i = 0; i < deg && working.degree >= 0; i++) {
        double angle = 2.0 * M_PI * i / (deg > 0 ? deg : 1);
        complex_t guess = complex_make(cos(angle), sin(angle));
        complex_t root = find_one_root_newton(&working, guess, 100);
        int is_new = 1;
        for (int j = 0; j < found; j++) {
            if (complex_mag(complex_sub(root, poles[j])) < 1e-6) { is_new = 0; mult[j]++; break; }
        }
        if (is_new) { poles[found] = root; mult[found] = 1; found++; }
        polynomial_t new_working = deflate_polynomial(&working, root);
        polynomial_free(&working);
        working = new_working;
    }
    polynomial_free(&working);
    *n_distinct = found;
    return 0;
}

int pf_residues_simple(const polynomial_t *N, const polynomial_t *D, const complex_t *poles, int n, complex_t *residues) {
    if (!N || !D || !poles || !residues || n <= 0) return -1;
    polynomial_t Dp = polynomial_derivative(D);
    for (int i = 0; i < n; i++) {
        complex_t Npk = polynomial_eval_complex(N, poles[i]);
        complex_t Dppk = polynomial_eval_complex(&Dp, poles[i]);
        residues[i] = (complex_mag(Dppk) < 1e-15) ? complex_make(NAN, NAN) : complex_div(Npk, Dppk);
    }
    polynomial_free(&Dp);
    return 0;
}

int pf_residues_repeated(const rational_func_t *F, complex_t pole, int m, complex_t *residues) {
    if (!F || m <= 0 || !residues) return -1;
    polynomial_t factor = polynomial_constant(1.0);
    for (int i = 0; i < m; i++) {
        double coeffs[2] = {-pole.real, 1.0};
        polynomial_t sp = polynomial_from_coeffs(coeffs, 2);
        polynomial_t nf = polynomial_mul(&factor, &sp);
        polynomial_free(&factor); polynomial_free(&sp);
        factor = nf;
    }
    polynomial_t quot, rem;
    polynomial_div(&F->den, &factor, &quot, &rem);
    polynomial_free(&rem);
    for (int k = 0; k < m; k++) {
        complex_t Np = polynomial_eval_complex(&F->num, pole);
        complex_t Dp = polynomial_eval_complex(&quot, pole);
        if (complex_mag(Dp) < 1e-15) { residues[k] = complex_make(NAN, NAN); continue; }
        complex_t Gp = complex_div(Np, Dp);
        if (k == 0) { residues[k] = Gp; }
        else {
            double h = 1e-6;
            complex_t hp = complex_make(h, 0.0);
            complex_t Gp_plus = complex_div(
                polynomial_eval_complex(&F->num, complex_add(pole, hp)),
                polynomial_eval_complex(&quot, complex_add(pole, hp)));
            complex_t deriv = complex_div(complex_sub(Gp_plus, Gp), hp);
            residues[k] = complex_make(deriv.real / k, deriv.imag / k);
        }
    }
    polynomial_free(&factor); polynomial_free(&quot);
    return 0;
}

int pf_decompose_laplace(const rational_func_t *F, partial_fraction_t *pf) {
    if (!F || !pf) return -1;
    if (F->den.degree < 0) return -1;
    rational_func_t proper_part;
    polynomial_t direct_term;
    if (F->num.degree >= F->den.degree) {
        polynomial_t quot, rem;
        if (polynomial_div(&F->num, &F->den, &quot, &rem) != 0) return -1;
        direct_term = quot;
        proper_part.num = rem; proper_part.den = polynomial_copy(&F->den);
    } else {
        direct_term = polynomial_zero();
        proper_part.num = polynomial_copy(&F->num);
        proper_part.den = polynomial_copy(&F->den);
    }
    int max_poles = F->den.degree;
    complex_t *poles_arr = (complex_t*)calloc(max_poles, sizeof(complex_t));
    int *mult_arr = (int*)calloc(max_poles, sizeof(int));
    int n_distinct = 0;
    if (pf_find_poles(&proper_part.den, poles_arr, mult_arr, &n_distinct) != 0) {
        free(poles_arr); free(mult_arr);
        rational_free(&proper_part); polynomial_free(&direct_term);
        return -1;
    }
    int total_terms = 0;
    for (int i = 0; i < n_distinct; i++) total_terms += mult_arr[i];
    pf->n_terms = total_terms;
    pf->terms = (pf_term_t*)calloc(total_terms, sizeof(pf_term_t));
    pf->direct = direct_term;
    int term_idx = 0;
    for (int i = 0; i < n_distinct; i++) {
        int m = mult_arr[i];
        if (m == 1) {
            complex_t Np = polynomial_eval_complex(&proper_part.num, poles_arr[i]);
            polynomial_t Dp = polynomial_derivative(&proper_part.den);
            complex_t Dpp = polynomial_eval_complex(&Dp, poles_arr[i]);
            polynomial_free(&Dp);
            pf->terms[term_idx].pole = poles_arr[i];
            pf->terms[term_idx].residue = complex_div(Np, Dpp);
            pf->terms[term_idx].multiplicity = 1;
            term_idx++;
        } else {
            complex_t *res = (complex_t*)calloc(m, sizeof(complex_t));
            pf_residues_repeated(&proper_part, poles_arr[i], m, res);
            for (int k = m - 1; k >= 0; k--) {
                pf->terms[term_idx].pole = poles_arr[i];
                pf->terms[term_idx].residue = res[k];
                pf->terms[term_idx].multiplicity = m - k;
                term_idx++;
            }
            free(res);
        }
    }
    free(poles_arr); free(mult_arr);
    rational_free(&proper_part);
    return 0;
}

int pf_decompose_z(const rational_func_t *X, partial_fraction_t *pf) {
    if (!X || !pf) return -1;
    if (X->den.degree < 0) return -1;
    rational_func_t proper_part;
    polynomial_t direct_term;
    if (X->num.degree >= X->den.degree) {
        polynomial_t quot, rem;
        if (polynomial_div(&X->num, &X->den, &quot, &rem) != 0) return -1;
        direct_term = quot;
        proper_part.num = rem; proper_part.den = polynomial_copy(&X->den);
    } else {
        direct_term = polynomial_zero();
        proper_part.num = polynomial_copy(&X->num);
        proper_part.den = polynomial_copy(&X->den);
    }
    int max_poles = X->den.degree;
    complex_t *poles_arr = (complex_t*)calloc(max_poles, sizeof(complex_t));
    int *mult_arr = (int*)calloc(max_poles, sizeof(int));
    int n_distinct = 0;
    if (pf_find_poles(&proper_part.den, poles_arr, mult_arr, &n_distinct) != 0) {
        free(poles_arr); free(mult_arr);
        rational_free(&proper_part); polynomial_free(&direct_term);
        return -1;
    }
    int total_terms = 0;
    for (int i = 0; i < n_distinct; i++) total_terms += mult_arr[i];
    pf->n_terms = total_terms;
    pf->terms = (pf_term_t*)calloc(total_terms, sizeof(pf_term_t));
    pf->direct = direct_term;
    int term_idx = 0;
    for (int i = 0; i < n_distinct; i++) {
        int m = mult_arr[i];
        if (m == 1) {
            complex_t Np = polynomial_eval_complex(&proper_part.num, poles_arr[i]);
            polynomial_t Dp = polynomial_derivative(&proper_part.den);
            complex_t Dpp = polynomial_eval_complex(&Dp, poles_arr[i]);
            polynomial_free(&Dp);
            pf->terms[term_idx].pole = poles_arr[i];
            pf->terms[term_idx].residue = complex_div(Np, Dpp);
            pf->terms[term_idx].multiplicity = 1;
            term_idx++;
        } else {
            complex_t *res = (complex_t*)calloc(m, sizeof(complex_t));
            pf_residues_repeated(&proper_part, poles_arr[i], m, res);
            for (int k = m - 1; k >= 0; k--) {
                pf->terms[term_idx].pole = poles_arr[i];
                pf->terms[term_idx].residue = res[k];
                pf->terms[term_idx].multiplicity = m - k;
                term_idx++;
            }
            free(res);
        }
    }
    free(poles_arr); free(mult_arr);
    rational_free(&proper_part);
    return 0;
}

rational_func_t pf_to_rational_func(const partial_fraction_t *pf) {
    rational_func_t result;
    if (pf->direct.degree >= 0) {
        result.num = polynomial_copy(&pf->direct);
    } else {
        result.num = polynomial_constant(0.0);
    }
    result.den = polynomial_constant(1.0);
    for (int i = 0; i < pf->n_terms; i++) {
        pf_term_t *t = &pf->terms[i];
        polynomial_t term_num = polynomial_constant(t->residue.real);
        polynomial_t term_den = polynomial_constant(1.0);
        for (int j = 0; j < t->multiplicity; j++) {
            double coeffs[2] = {-t->pole.real, 1.0};
            polynomial_t factor = polynomial_from_coeffs(coeffs, 2);
            polynomial_t nd = polynomial_mul(&term_den, &factor);
            polynomial_free(&term_den); polynomial_free(&factor);
            term_den = nd;
        }
        rational_func_t trf; trf.num = term_num; trf.den = term_den;
        rational_func_t nr = rational_add(&result, &trf);
        rational_free(&result); rational_free(&trf);
        result = nr;
    }
    return result;
}