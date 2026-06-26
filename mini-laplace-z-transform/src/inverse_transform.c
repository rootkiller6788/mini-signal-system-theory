/*
 * inverse_transform.c - Inverse Laplace and Z-Transform Computation
 *
 * Implements inverse transforms using partial fraction expansion
 * and table lookup, plus numerical methods for cases where analytic
 * inversion is not practical.
 *
 * Knowledge Coverage:
 *   L5: Inverse Laplace transform via partial fractions + table lookup
 *   L5: Inverse Z-transform via power series, partial fractions, residue method
 *   L8: Numerical inverse Laplace transform (Gaver-Stehfest method)
 *
 * Reference: Oppenheim & Willsky, "Signals and Systems" (1997), Ch.9
 *            Oppenheim & Schafer, "DSP" (2010), Ch.3
 *            Cohen, "Numerical Methods for Laplace Transform Inversion" (2007)
 */

#include "laplace_z_transform_core.h"
#include "partial_fraction.h"
#include "transfer_function.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * ========================================================================
 * L5: Analytic Inverse Laplace Transform
 * ========================================================================
 *
 * Given F(s) = N(s)/D(s), compute f(t) = L^{-1}{F(s)} at time t.
 *
 * Method:
 *   1. Partial fraction expansion: F(s) = sum r_k/(s-p_k)^m_k + direct term
 *   2. Table lookup:
 *      - r/(s-p)     <-> r*exp(p*t)
 *      - r/(s-p)^2   <-> r*t*exp(p*t)
 *      - r/(s-p)^m   <-> r*t^{m-1}/(m-1)! * exp(p*t)
 *      - Direct term <-> impulses and derivatives (at t=0 only)
 *
 * For complex poles, the result is always real (conjugate pairs sum to real).
 */

double inverse_laplace_analytic(const rational_func_t *F, double t) {
    /*
     * Compute f(t) analytically for a rational Laplace transform F(s).
     *
     * Algorithm: partial fraction expansion + table lookup.
     * O(n) where n = number of partial fraction terms.
     *
     * Reference: Oppenheim & Willsky (1997), Table 9.2
     */
    return tf_impulse_response(F, t);
}

/*
 * Compute f(t) at multiple time points for plotting/analysis.
 *
 * Fills array 'values' with f(t_i) for each time point in 'times'.
 */
void inverse_laplace_vector(const rational_func_t *F, const double *times,
                             double *values, int n_points) {
    for (int i = 0; i < n_points; i++) {
        values[i] = inverse_laplace_analytic(F, times[i]);
    }
}

/*
 * ========================================================================
 * L5: Analytic Inverse Z-Transform
 * ========================================================================
 *
 * Given X(z^{-1}) = N(z^{-1})/D(z^{-1}), compute x[n] = Z^{-1}{X(z)}.
 *
 * Methods:
 *   1. Partial fraction expansion (preferred for analytic solutions)
 *   2. Power series (long division) - good for first few samples
 *   3. Residue method (contour integration)
 *
 * For causal sequences:
 *   X(z^{-1}) = sum r_k/(1 - p_k*z^{-1})
 *   x[n] = sum r_k * p_k^n * u[n]
 */

double inverse_z_analytic(const rational_func_t *X, int n) {
    /*
     * Compute x[n] analytically using partial fraction expansion.
     *
     * Reference: Oppenheim & Schafer (2010), Section 3.3.2
     */
    return tf_discrete_impulse_response(X, n);
}

void inverse_z_vector(const rational_func_t *X, int n_max,
                       double *values) {
    /*
     * Compute first n_max+1 samples of the sequence x[n].
     */
    for (int n = 0; n <= n_max; n++) {
        values[n] = inverse_z_analytic(X, n);
    }
}

/*
 * Inverse Z-transform via power series (long division).
 *
 * This method directly computes the first N terms of x[n] by
 * dividing the numerator by the denominator polynomial in z^{-1}.
 *
 * Algorithm: synthetic division / deconvolution.
 * O(N * deg(den)) for N samples.
 *
 * This is particularly useful when the partial fraction expansion
 * is difficult (e.g., high-order systems).
 *
 * Reference: Oppenheim & Schafer (2010), Section 3.3.1
 */
int inverse_z_power_series(const rational_func_t *X, double *x, int N) {
    /*
     * Compute first N samples of x[n] using power series expansion.
     *
     * X(z^{-1}) = (b0 + b1*z^{-1} + ... + bM*z^{-M}) / (1 + a1*z^{-1} + ... + aN*z^{-N})
     *
     * Requires denominator constant term = 1 (normalized form).
     *
     * Returns 0 on success, -1 if denominator constant term is zero.
     */
    if (!X || !x || N <= 0) return -1;
    if (X->den.degree < 0 || X->num.degree < 0) return -1;

    /* Check denominator normalization */
    double a0 = X->den.coeff[0];
    if (fabs(a0) < 1e-15) return -1;

    int M = X->num.degree;
    int N_den = X->den.degree;

    /* Allocate normalized coefficient arrays */
    double *b = (double*)calloc(M + 1, sizeof(double));
    double *a = (double*)calloc(N_den + 1, sizeof(double));

    if (!b || !a) { free(b); free(a); return -1; }

    for (int i = 0; i <= M; i++) b[i] = X->num.coeff[i] / a0;
    for (int i = 0; i <= N_den; i++) a[i] = X->den.coeff[i] / a0;

    /* Long division (deconvolution):
     * x[n] = b[n] - sum_{k=1}^{min(n, N_den)} a[k] * x[n-k] */
    for (int n = 0; n < N; n++) {
        double sum = (n <= M) ? b[n] : 0.0;
        for (int k = 1; k <= N_den && k <= n; k++) {
            sum -= a[k] * x[n - k];
        }
        x[n] = sum;
    }

    free(b); free(a);
    return 0;
}

/*
 * Inverse Z-transform via the residue method (contour integration).
 *
 * x[n] = sum of residues of X(z) * z^{n-1} at poles inside the ROC.
 *
 * For a rational X(z) with simple poles at z = p_k:
 * x[n] = sum_k Res[X(z)*z^{n-1}, z=p_k]
 *      = sum_k (z-p_k)*X(z)*z^{n-1} |_{z=p_k}
 *
 * This is the formal complex analysis method.
 *
 * Reference: Churchill & Brown, "Complex Variables" (2009), Ch.7
 */
double inverse_z_residue(const rational_func_t *X, int n) {
    /*
     * Compute x[n] using the residue method.
     *
     * For X(z^{-1}) with simple poles:
     * Convert to X(z) in positive powers of z, find residues at poles.
     */
    if (!X || n < 0) return 0.0;
    if (X->den.degree < 0) return 0.0;

    /* The residue method gives the same result as partial fractions.
     * We delegate to the partial fraction method for robustness. */
    return tf_discrete_impulse_response(X, n);
}

/*
 * ========================================================================
 * L8: Numerical Inverse Laplace Transform (Gaver-Stehfest Method)
 * ========================================================================
 *
 * The Gaver-Stehfest method is a numerical algorithm for computing
 * the inverse Laplace transform when an analytic solution is impractical
 * (e.g., for non-rational F(s) or very high-order systems).
 *
 * Formula:
 *   f(t) approx (ln(2)/t) * sum_{k=1}^{2M} V_k * F(k*ln(2)/t)
 *
 * where V_k are the Stehfest coefficients:
 *   V_k = (-1)^{M+k} * sum_{j=floor((k+1)/2)}^{min(k,M)}
 *         j^{M+1}/(M!) * C(M,j) * C(2j,j) * C(j,k-j)
 *
 * M is the number of terms (typically 6-18, even).
 * Larger M gives better accuracy but requires higher precision.
 *
 * This is an L8 (advanced) method for numerical transform inversion.
 *
 * Reference: Stehfest, "Numerical Inversion of Laplace Transforms",
 *            Communications of the ACM (1970)
 *            Cohen, "Numerical Methods for Laplace Transform Inversion" (2007)
 */

/* Compute binomial coefficient C(n, k) */
static double binomial_coeff(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    /* Use multiplicative formula to avoid overflow */
    double result = 1.0;
    for (int i = 1; i <= k; i++) {
        result *= (n - k + i) / (double)i;
    }
    return result;
}

/* Compute factorial */
static double factorial(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; i++) result *= i;
    return result;
}

int inverse_laplace_gaver_stehfest_coeffs(int M, double *V) {
    /*
     * Compute Gaver-Stehfest coefficients V_k for k = 1..2M.
     *
     * M must be even, typically 8-16.
     * V_k = (-1)^{M+k} * sum_{j=(k+1)/2}^{min(k,M)}
     *       j^{M+1} / M! * C(M,j) * C(2j,j) * C(j,k-j)
     *
     * Returns 0 on success, -1 on error.
     */
    if (!V || M < 2 || M % 2 != 0) return -1;

    double M_fact = factorial(M);

    for (int k = 1; k <= 2*M; k++) {
        double sum = 0.0;
        int j_start = (k + 1) / 2;
        int j_end = (k < M) ? k : M;

        for (int j = j_start; j <= j_end; j++) {
            double term = pow((double)j, M + 1) / M_fact;
            term *= binomial_coeff(M, j);
            term *= binomial_coeff(2*j, j);
            term *= binomial_coeff(j, k - j);
            sum += term;
        }

        V[k - 1] = ((M + k) % 2 == 0 ? 1.0 : -1.0) * sum;
    }

    return 0;
}

/*
 * User-provided function pointer type for F(s) evaluation.
 * Allows the Gaver-Stehfest method to work with arbitrary F(s).
 */
typedef complex_t (*laplace_func_t)(complex_t s, void *user_data);

double inverse_laplace_gaver_stehfest(laplace_func_t F, void *user_data,
                                       double t, int M) {
    /*
     * Compute f(t) using the Gaver-Stehfest numerical inversion method.
     *
     * @param F         Function evaluating F(s) at complex s
     * @param user_data User data passed to F
     * @param t         Time point (t > 0)
     * @param M         Number of terms (even, typically 8-16)
     * @return          f(t) estimate
     *
     * Formula:
     *   f(t) = (ln(2)/t) * sum_{k=1}^{2M} V_k * Re{F(k*ln(2)/t)}
     *
     * Note: Only the real part of F(s) at real s is needed for the
     * standard Gaver-Stehfest formula. This is an advantage over methods
     * requiring complex arithmetic.
     */
    if (!F || t <= 0.0 || M < 2 || M % 2 != 0) return NAN;

    double *V = (double*)malloc(2 * M * sizeof(double));
    if (!V) return NAN;

    if (inverse_laplace_gaver_stehfest_coeffs(M, V) != 0) {
        free(V);
        return NAN;
    }

    double ln2_over_t = log(2.0) / t;
    double sum = 0.0;

    for (int k = 1; k <= 2*M; k++) {
        double s_real = k * ln2_over_t;
        complex_t s = complex_make(s_real, 0.0);
        complex_t Fs = F(s, user_data);
        sum += V[k - 1] * Fs.real;
    }

    free(V);
    return ln2_over_t * sum;
}

/*
 * Convenience wrapper for rational function F(s) using Gaver-Stehfest.
 */
static complex_t rational_F_wrapper(complex_t s, void *user_data) {
    rational_func_t *F = (rational_func_t*)user_data;
    return rational_eval(F, s);
}

double inverse_laplace_numerical(const rational_func_t *F, double t, int M) {
    /*
     * Numerical inverse Laplace for a rational function using Gaver-Stehfest.
     *
     * This provides an alternative to analytic inversion, useful for
     * verification or when the partial fraction decomposition fails.
     */
    return inverse_laplace_gaver_stehfest(rational_F_wrapper, (void*)F, t, M);
}

/*
 * ========================================================================
 * L8: Chirp Z-Transform (CZT)
 * ========================================================================
 *
 * The Chirp Z-Transform (CZT) evaluates the Z-transform on a spiral
 * contour in the z-plane, enabling zoomed frequency analysis.
 *
 * X(z_k) = sum_{n=0}^{N-1} x[n] * z_k^{-n}
 * where z_k = A * W^{-k}, k = 0, 1, ..., M-1
 *
 * A = A0 * exp(j*theta0): starting point
 * W = W0 * exp(-j*phi0): ratio between successive points
 *
 * The CZT generalizes the DFT (which evaluates on the unit circle only).
 * It's an L8 advanced topic.
 *
 * Reference: Rabiner, Schafer, Rader, "The Chirp Z-Transform Algorithm"
 *            IEEE Trans. Audio Electroacoustics (1969)
 */

void chirp_z_transform(const double *x, int N,
                        double A0, double theta0,
                        double W0, double phi0,
                        int M, complex_t *X) {
    /*
     * Compute the Chirp Z-Transform of sequence x[n].
     *
     * @param x       Input sequence (length N)
     * @param N       Sequence length
     * @param A0      Starting radius
     * @param theta0  Starting angle (radians)
     * @param W0      Spiral expansion factor (>1: expanding, <1: contracting)
     * @param phi0    Angular step (radians)
     * @param M       Number of output points
     * @param X       Output array (length M, complex)
     *
     * Direct computation: O(N*M). Can be accelerated to O((N+M)*log(N+M))
     * using FFT convolution (not implemented here).
     */
    if (!x || !X || N <= 0 || M <= 0) return;

    complex_t A = complex_make(A0 * cos(theta0), A0 * sin(theta0));
    complex_t W = complex_make(W0 * cos(-phi0), W0 * sin(-phi0));

    for (int k = 0; k < M; k++) {
        /* z_k = A * W^{-k} */
        complex_t Wk = complex_pow(W, -k);
        complex_t zk = complex_mul(A, Wk);

        /* X(z_k) = sum x[n] * z_k^{-n} */
        complex_t sum = complex_make(0.0, 0.0);
        complex_t zk_pow = complex_make(1.0, 0.0);  /* z_k^0 = 1 */
        complex_t zk_inv = complex_inv(zk);         /* z_k^{-1} */

        for (int n = 0; n < N; n++) {
            /* sum += x[n] * zk_pow^{-n} = x[n] * (zk^{-1})^n */
            complex_t term = complex_make(x[n] * zk_pow.real, x[n] * zk_pow.imag);
            sum = complex_add(sum, term);
            zk_pow = complex_mul(zk_pow, zk_inv);
        }

        X[k] = sum;
    }
}

/*
 * ========================================================================
 * L6: Difference Equation Solution via Z-Transform
 * ========================================================================
 *
 * Solve a linear constant-coefficient difference equation:
 *   sum_{k=0}^{N} a_k * y[n-k] = sum_{k=0}^{M} b_k * x[n-k]
 *
 * Using the Z-transform:
 *   Y(z) * A(z) = X(z) * B(z)
 *   => Y(z) = H(z) * X(z)  where H(z) = B(z)/A(z)
 *
 * For zero initial conditions and a specific input, compute y[n].
 */

/**
 * Compute the output sequence y[n] for a given input x[n] using the
 * difference equation directly (time-domain recursion).
 *
 * This implements the canonical Direct Form I realization.
 *
 * Reference: Oppenheim & Schafer (2010), Section 6.2
 */
void difference_equation_solve(const double *b, int M,
                                const double *a, int N,
                                const double *x, int input_len,
                                double *y) {
    /*
     * y[n] = sum_{k=0}^{M} b_k*x[n-k] - sum_{k=1}^{N} a_k*y[n-k]
     *
     * @param b         Numerator coefficients (feedforward)
     * @param M         Numerator order
     * @param a         Denominator coefficients (feedback), a[0] = 1
     * @param N         Denominator order
     * @param x         Input sequence
     * @param input_len Length of input
     * @param y         Output sequence (same length)
     *
     * Initial rest conditions assumed (y[n]=0, x[n]=0 for n<0).
     */
    if (!b || !a || !x || !y || input_len <= 0) return;

    /* Normalize: ensure a[0] = 1 */
    double a0 = a[0];
    if (fabs(a0) < 1e-15) return;

    for (int n = 0; n < input_len; n++) {
        double sum_b = 0.0;
        for (int k = 0; k <= M && k <= n; k++) {
            sum_b += b[k] * x[n - k] / a0;
        }

        double sum_a = 0.0;
        for (int k = 1; k <= N && k <= n; k++) {
            sum_a += a[k] * y[n - k] / a0;
        }

        y[n] = sum_b - sum_a;
    }
}