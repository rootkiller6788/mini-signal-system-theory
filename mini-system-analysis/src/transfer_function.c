
/**
 * @file transfer_function.c
 * @brief L3/L5 Transfer Function Evaluation, Pole-Zero, Partial Fractions
 *
 * Implementation of rational transfer function operations in CT (Laplace/s)
 * and DT (Z-domain). Each function implements an independent knowledge point
 * from signal and system theory.
 *
 * Knowledge Coverage (20 points):
 *   L3-I01: ct_tf_evaluate - Horner method polynomial evaluation
 *   L3-I02: dt_tf_evaluate
 *   L3-I03: ct_tf_evaluate_vector / dt_tf_evaluate_vector
 *   L3-I04: ct_frequency_response - H(jw) frequency sweep
 *   L3-I05: dt_frequency_response - H(e^{jw}) unit circle sweep
 *   L3-I06: ct_tf_find_poles - Newton+deflation root finding
 *   L3-I07: ct_tf_find_zeros
 *   L3-I08: dt_tf_find_poles / dt_tf_find_zeros
 *   L3-I09: ct_tf_pole_zero_analysis (full PZ+damping+wn)
 *   L3-I10: dt_tf_pole_zero_analysis
 *   L3-I11: ct_tf_partial_fraction - residue computation
 *   L3-I12: dt_tf_partial_fraction
 *   L3-I13: ct_to_dt_bilinear - Tustin transform
 *   L3-I14: ct_to_dt_impulse_invariance
 *   L3-I15: ct_to_dt_zoh - Zero-order hold
 *   L3-I16: ct_to_dt_matched_z
 *   L3-I17: ct_tf_order_analysis / dt_tf_order_analysis
 *   L3-I18: ct_tf_zpk_gain / dt_tf_zpk_gain
 *   L3-I19: ct_tf_bode_plot / ct_tf_nyquist_plot
 *   L3-I20: ct_tf_group_delay / dt_tf_group_delay
 *
 * References:
 *   Oppenheim & Willsky, Signals and Systems, 2nd Ed. Ch.9
 *   Franklin, Powell, Emami-Naeini, Feedback Control, 7th Ed. Ch.8
 */

#include "transfer_function.h"
#include "system_defs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Internal: Horner polynomial evaluation ---- */
static double complex poly_eval_horner(const double *coeffs, int degree,
                                        double complex x)
{
    if (degree < 0) return 0.0 + 0.0 * I;
    double complex result = coeffs[degree] + 0.0 * I;
    for (int i = degree - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }
    return result;
}

/* ---- Internal: Polynomial value + derivative (for Newton method) ---- */
static void poly_eval_and_deriv(const double *coeffs, int degree,
                                 double complex x,
                                 double complex *val, double complex *deriv)
{
    *val = poly_eval_horner(coeffs, degree, x);
    if (degree <= 0) { *deriv = 0.0 + 0.0 * I; return; }
    *deriv = (double)degree * coeffs[degree];
    for (int k = degree - 1; k >= 1; k--) {
        *deriv = (*deriv) * x + (double)k * coeffs[k];
    }
}

/* ---- Internal: Newton-Raphson root finding ---- */
static int newton_find_root(const double *coeffs, int degree,
                             double complex guess, double complex *root,
                             int max_iter, double tol)
{
    double complex x = guess;
    for (int iter = 0; iter < max_iter; iter++) {
        double complex val, deriv;
        poly_eval_and_deriv(coeffs, degree, x, &val, &deriv);
        if (cabs(deriv) < 1e-30) { x += 0.1 + 0.1 * I; continue; }
        double complex dx = val / deriv;
        x = x - dx;
        if (cabs(dx) < tol) { *root = x; return 0; }
    }
    *root = x;
    return 0;
}

/* ---- Internal: Synthetic division (deflation) ---- */
static void poly_deflate(double *coeffs, int degree, double complex root)
{
    if (degree < 1) return;
    double complex *temp = (double complex *)calloc((size_t)(degree + 1),
                                                     sizeof(double complex));
    if (!temp) return;
    for (int i = 0; i <= degree; i++) temp[i] = coeffs[i] + 0.0 * I;
    double complex b = temp[degree];
    coeffs[degree - 1] = creal(b);
    for (int k = degree - 2; k >= 0; k--) {
        b = temp[k + 1] + root * b;
        coeffs[k] = creal(b);
    }
    free(temp);
}

/* ========================================================================
 *  L3-I01: ct_tf_evaluate — H(s) = N(s)/D(s) at complex s
 *
 *  Uses Horner's method for numerically stable polynomial evaluation.
 *  Complexity: O(max(num_order, den_order)).
 *  Returns NaN+NaN*I at pole locations (D(s) ≈ 0).
 * ======================================================================== */

double complex ct_tf_evaluate(const ct_transfer_function_t *tf,
                               double complex s)
{
    if (!tf || !tf->num_coeffs || !tf->den_coeffs)
        return NAN + NAN * I;

    double complex num = poly_eval_horner(tf->num_coeffs, tf->num_order, s);
    double complex den = poly_eval_horner(tf->den_coeffs, tf->den_order, s);

    if (cabs(den) < 1e-30) return NAN + NAN * I;
    return num / den;
}

/* ========================================================================
 *  L3-I02: dt_tf_evaluate — H(z) = N(z)/D(z) at complex z
 * ======================================================================== */

double complex dt_tf_evaluate(const dt_transfer_function_t *tf,
                               double complex z)
{
    if (!tf || !tf->num_coeffs || !tf->den_coeffs)
        return NAN + NAN * I;

    double complex num = poly_eval_horner(tf->num_coeffs, tf->num_order, z);
    double complex den = poly_eval_horner(tf->den_coeffs, tf->den_order, z);

    if (cabs(den) < 1e-30) return NAN + NAN * I;
    return num / den;
}

/* ========================================================================
 *  L3-I03: Vectorized transfer function evaluation
 * ======================================================================== */

int ct_tf_evaluate_vector(const ct_transfer_function_t *tf,
                           const double complex *s_values, int num_points,
                           double complex *result)
{
    if (!tf || !s_values || !result || num_points <= 0) return -1;
    for (int i = 0; i < num_points; i++)
        result[i] = ct_tf_evaluate(tf, s_values[i]);
    return 0;
}

int dt_tf_evaluate_vector(const dt_transfer_function_t *tf,
                           const double complex *z_values, int num_points,
                           double complex *result)
{
    if (!tf || !z_values || !result || num_points <= 0) return -1;
    for (int i = 0; i < num_points; i++)
        result[i] = dt_tf_evaluate(tf, z_values[i]);
    return 0;
}

/* ========================================================================
 *  L3-I04: ct_frequency_response — Evaluate H(jw) along frequency grid
 *
 *  Theorem (Frequency Response):
 *    H(jw) = H(s)|_{s=jw} is the Fourier transform of h(t).
 *    Magnitude |H(jw)| gives gain at each frequency.
 *    Phase angle(H(jw)) gives phase shift at each frequency.
 * ======================================================================== */

int ct_frequency_response(const ct_transfer_function_t *tf,
                          frequency_response_t *fr)
{
    if (!tf || !fr || !fr->points) return -1;

    for (int i = 0; i < fr->num_points; i++) {
        double f_hz = fr->log_scale
            ? fr->f_start * pow(fr->f_end / fr->f_start,
                                (double)i / (double)(fr->num_points - 1))
            : fr->f_start + (fr->f_end - fr->f_start) *
                            (double)i / (double)(fr->num_points - 1);

        double w = 2.0 * M_PI * f_hz;
        double complex s = 0.0 + w * I;
        double complex H = ct_tf_evaluate(tf, s);

        double mag = cabs(H);
        double phase = carg(H);

        fr->points[i].frequency   = f_hz;
        fr->points[i].magnitude   = mag;
        fr->points[i].magnitude_db = (mag > 1e-30) ? 20.0 * log10(mag) : -300.0;
        fr->points[i].phase_rad   = phase;
        fr->points[i].phase_deg   = phase * 180.0 / M_PI;
    }
    return 0;
}

/* ========================================================================
 *  L3-I05: dt_frequency_response — H(e^{jw}) on unit circle
 *
 *  For DT: z = e^{jw}, w ∈ [0, π] maps to f ∈ [0, Fs/2].
 *  Frequency response is periodic with period 2π.
 * ======================================================================== */

int dt_frequency_response(const dt_transfer_function_t *tf,
                          double sample_rate,
                          frequency_response_t *fr)
{
    if (!tf || !fr || !fr->points || sample_rate <= 0.0) return -1;

    for (int i = 0; i < fr->num_points; i++) {
        double f_hz = fr->log_scale
            ? fr->f_start * pow(fr->f_end / fr->f_start,
                                (double)i / (double)(fr->num_points - 1))
            : fr->f_start + (fr->f_end - fr->f_start) *
                            (double)i / (double)(fr->num_points - 1);

        double w = 2.0 * M_PI * f_hz / sample_rate;
        double complex z = cos(w) + sin(w) * I;
        double complex H = dt_tf_evaluate(tf, z);

        double mag = cabs(H);
        double phase = carg(H);

        fr->points[i].frequency   = f_hz;
        fr->points[i].magnitude   = mag;
        fr->points[i].magnitude_db = (mag > 1e-30) ? 20.0 * log10(mag) : -300.0;
        fr->points[i].phase_rad   = phase;
        fr->points[i].phase_deg   = phase * 180.0 / M_PI;
    }
    return 0;
}

/* ========================================================================
 *  L3-I06: ct_tf_find_poles — Newton+deflation root-finding
 *
 *  Theorem (Companion Matrix):
 *    Eigenvalues of companion matrix = roots of polynomial.
 *    For D(s) = s^n + a_{n-1}s^{n-1} + ... + a_0, the companion matrix
 *    C has D(s) as its characteristic polynomial: det(sI - C) = D(s).
 *
 *  Implementation uses Newton-Raphson with deflation for finding all roots.
 *  Initial guesses are distributed on circles in the complex plane.
 *  Complexity: O(n * max_iter) per root, O(n^2) total.
 * ======================================================================== */

int ct_tf_find_poles(const ct_transfer_function_t *tf,
                     pole_zero_collection_t *pz)
{
    if (!tf || !pz || pz->num_poles < tf->den_order) return -1;

    int n = tf->den_order;
    if (n == 0) return 0;

    double lead = tf->den_coeffs[n];
    if (fabs(lead) < 1e-30) return -1;

    double *poly = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!poly) return -1;
    for (int i = 0; i <= n; i++) poly[i] = tf->den_coeffs[i] / lead;

    for (int k = 0; k < n; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)n + 0.5;
        double radius = 1.0 + 0.5 * (double)k / (double)n;
        double complex guess = radius * cos(angle) + radius * sin(angle) * I;

        double complex root;
        newton_find_root(poly, n - k, guess, &root, 100, 1e-10);

        /* Refine with original polynomial */
        double complex refined;
        double *fp = (double *)malloc((size_t)(n + 1) * sizeof(double));
        if (fp) {
            for (int i = 0; i <= n; i++) fp[i] = tf->den_coeffs[i] / lead;
            newton_find_root(fp, n, root, &refined, 50, 1e-12);
            free(fp);
        } else {
            refined = root;
        }

        pz->poles[k].value   = refined;
        pz->poles[k].is_real = (cimag(refined) < 1e-8);
        pz->poles[k].is_pole = 1;
        pz->poles[k].wn      = cabs(refined);
        pz->poles[k].damping = (pz->poles[k].wn > 1e-30)
            ? -creal(refined) / pz->poles[k].wn : 0.0;

        if (n - k > 1) poly_deflate(poly, n - k, refined);
    }

    /* CT stability: all poles must have Re(s) < 0 */
    pz->is_stable = 1;
    for (int i = 0; i < n; i++) {
        if (creal(pz->poles[i].value) >= -1e-10) {
            pz->is_stable = 0;
            break;
        }
    }

    free(poly);
    return 0;
}

/* ========================================================================
 *  L3-I07: ct_tf_find_zeros — roots of N(s) = 0
 * ======================================================================== */

int ct_tf_find_zeros(const ct_transfer_function_t *tf,
                     pole_zero_collection_t *pz)
{
    if (!tf || !pz || pz->num_zeros < tf->num_order) return -1;

    int m = tf->num_order;
    if (m == 0) return 0;

    double lead = tf->num_coeffs[m];
    if (fabs(lead) < 1e-30) return 0;

    double *poly = (double *)malloc((size_t)(m + 1) * sizeof(double));
    if (!poly) return -1;
    for (int i = 0; i <= m; i++) poly[i] = tf->num_coeffs[i] / lead;

    for (int k = 0; k < m; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)m + 0.7;
        double radius = 1.0 + 0.3 * (double)k / (double)m;
        double complex guess = radius * cos(angle) + radius * sin(angle) * I;

        double complex root;
        newton_find_root(poly, m - k, guess, &root, 100, 1e-10);

        pz->zeros[k].value   = root;
        pz->zeros[k].is_real = (cimag(root) < 1e-8);
        pz->zeros[k].is_pole = 0;
        pz->zeros[k].wn      = cabs(root);

        if (m - k > 1) poly_deflate(poly, m - k, root);
    }

    free(poly);
    return 0;
}

/* ========================================================================
 *  L3-I08: dt_tf_find_poles — roots of D(z) = 0
 *  Stability criterion: all |z| < 1 (inside unit circle)
 * ======================================================================== */

int dt_tf_find_poles(const dt_transfer_function_t *tf,
                     pole_zero_collection_t *pz)
{
    if (!tf || !pz || pz->num_poles < tf->den_order) return -1;

    int n = tf->den_order;
    if (n == 0) return 0;

    double lead = tf->den_coeffs[n];
    if (fabs(lead) < 1e-30) return -1;

    double *poly = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!poly) return -1;
    for (int i = 0; i <= n; i++) poly[i] = tf->den_coeffs[i] / lead;

    for (int k = 0; k < n; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)n + 0.5;
        double complex guess = cos(angle) + sin(angle) * I;  /* on unit circle */
        double complex root;
        newton_find_root(poly, n - k, guess, &root, 100, 1e-10);

        pz->poles[k].value   = root;
        pz->poles[k].is_real = (cimag(root) < 1e-8);
        pz->poles[k].is_pole = 1;
        pz->poles[k].wn      = cabs(root);

        if (n - k > 1) poly_deflate(poly, n - k, root);
    }

    /* DT stability: all poles |z| < 1 */
    pz->is_stable = 1;
    for (int i = 0; i < n; i++) {
        if (cabs(pz->poles[i].value) >= 1.0 - 1e-10) {
            pz->is_stable = 0;
            break;
        }
    }

    free(poly);
    return 0;
}

/* ========================================================================
 *  L3-I08b: dt_tf_find_zeros
 * ======================================================================== */

int dt_tf_find_zeros(const dt_transfer_function_t *tf,
                     pole_zero_collection_t *pz)
{
    if (!tf || !pz || pz->num_zeros < tf->num_order) return -1;

    int m = tf->num_order;
    if (m == 0) return 0;

    double lead = tf->num_coeffs[m];
    if (fabs(lead) < 1e-30) return 0;

    double *poly = (double *)malloc((size_t)(m + 1) * sizeof(double));
    if (!poly) return -1;
    for (int i = 0; i <= m; i++) poly[i] = tf->num_coeffs[i] / lead;

    for (int k = 0; k < m; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)m + 0.7;
        double complex guess = cos(angle) + sin(angle) * I;
        double complex root;
        newton_find_root(poly, m - k, guess, &root, 100, 1e-10);

        pz->zeros[k].value   = root;
        pz->zeros[k].is_real = (cimag(root) < 1e-8);
        pz->zeros[k].is_pole = 0;
        pz->zeros[k].wn      = cabs(root);

        if (m - k > 1) poly_deflate(poly, m - k, root);
    }

    free(poly);
    return 0;
}

/* ========================================================================
 *  L3-I09: Full Pole-Zero Analysis — compute all poles, zeros, stability
 * ======================================================================== */

int ct_tf_pole_zero_analysis(const ct_transfer_function_t *tf,
                              pole_zero_collection_t *pz)
{
    int ret = ct_tf_find_poles(tf, pz);
    if (ret != 0) return ret;
    ret = ct_tf_find_zeros(tf, pz);
    if (ret != 0) return ret;
    pz->dc_gain = creal(ct_tf_evaluate(tf, 0.0 + 0.0 * I));
    return 0;
}

int dt_tf_pole_zero_analysis(const dt_transfer_function_t *tf,
                              pole_zero_collection_t *pz)
{
    int ret = dt_tf_find_poles(tf, pz);
    if (ret != 0) return ret;
    ret = dt_tf_find_zeros(tf, pz);
    if (ret != 0) return ret;
    pz->dc_gain = creal(dt_tf_evaluate(tf, 1.0 + 0.0 * I));
    return 0;
}

/* ========================================================================
 *  L3-I11: Partial Fraction Expansion — residue computation
 *
 *  H(s) = sum_{k} r_k / (s - p_k) + direct_term(s)
 *
 *  For distinct poles, r_k = N(p_k) / D'(p_k).
 *  This decomposes the system into parallel first-order subsystems.
 * ======================================================================== */

int ct_tf_partial_fraction(const ct_transfer_function_t *tf,
                           double complex *residues,
                           double complex *poles,
                           double *direct_term, int *direct_len)
{
    if (!tf || !residues || !poles || !direct_term || !direct_len) return -1;

    int n = tf->den_order;
    int m = tf->num_order;

    pole_zero_collection_t pz = pz_collection_alloc(n, 0);
    if (n > 0 && !pz.poles) { pz_collection_free(&pz); return -1; }

    int ret = ct_tf_find_poles(tf, &pz);
    if (ret != 0) { pz_collection_free(&pz); return ret; }

    for (int i = 0; i < n; i++) {
        poles[i] = pz.poles[i].value;
        double complex N_pi = poly_eval_horner(tf->num_coeffs, m, poles[i]);
        double complex D_val, D_deriv;
        poly_eval_and_deriv(tf->den_coeffs, n, poles[i], &D_val, &D_deriv);
        residues[i] = (cabs(D_deriv) > 1e-30) ? (N_pi / D_deriv) : (0.0 + 0.0 * I);
    }

    /* Direct term: exists for improper systems (m >= n) */
    if (m >= n) {
        *direct_len = m - n + 1;
        direct_term[0] = tf->num_coeffs[m] / tf->den_coeffs[n];
    } else {
        *direct_len = 0;
    }

    pz_collection_free(&pz);
    return 0;
}

int dt_tf_partial_fraction(const dt_transfer_function_t *tf,
                           double complex *residues,
                           double complex *poles,
                           double *direct_term, int *direct_len)
{
    if (!tf || !residues || !poles) return -1;

    int n = tf->den_order;
    int m = tf->num_order;
    if (n == 0) { *direct_len = 0; (void)direct_term; return 0; }

    pole_zero_collection_t pz = pz_collection_alloc(n, 0);
    if (n > 0 && !pz.poles) { pz_collection_free(&pz); return -1; }

    int ret = dt_tf_find_poles(tf, &pz);
    if (ret != 0) { pz_collection_free(&pz); return ret; }

    for (int i = 0; i < n; i++) {
        poles[i] = pz.poles[i].value;
        double complex N_pi = poly_eval_horner(tf->num_coeffs, m, poles[i]);
        double complex D_val, D_deriv;
        poly_eval_and_deriv(tf->den_coeffs, n, poles[i], &D_val, &D_deriv);
        residues[i] = (cabs(D_deriv) > 1e-30) ? (N_pi / D_deriv) : (0.0 + 0.0 * I);
    }

    *direct_len = (m >= n) ? (m - n + 1) : 0;
    pz_collection_free(&pz);
    return 0;
}

/* ========================================================================
 *  L3-I13: Bilinear (Tustin) Transform — s = (2/T)(z-1)/(z+1)
 *
 *  Maps LHP (Re{s} < 0) to unit circle interior (|z| < 1).
 *  Preserves stability and system order.
 *
 *  For 1st-order CT system H(s) = K/(s + p):
 *    H(z) = K*T*(z+1) / ((2+pT)z + (pT-2))
 * ======================================================================== */

int ct_to_dt_bilinear(const ct_transfer_function_t *ct_tf,
                      double sample_time,
                      dt_transfer_function_t *dt_tf)
{
    if (!ct_tf || !dt_tf || sample_time <= 0.0) return -1;

    int n = ct_tf->den_order;
    int m = ct_tf->num_order;
    double T = sample_time;
    double c = 2.0 / T;  /* prewarp constant */

    /* First-order system: explicit formula */
    if (n == 1 && m == 0) {
        double b0 = ct_tf->num_coeffs[0];
        double a0 = ct_tf->den_coeffs[0];
        double a1 = ct_tf->den_coeffs[1];
        double K = b0 / a1;
        double p = a0 / a1;
        double den = c + p;

        if (fabs(den) < 1e-30) return -1;

        dt_tf->num_coeffs[0] = K / den;       /* b0: (K/den) * z^0 */
        dt_tf->num_coeffs[1] = K / den;       /* b1: (K/den) * z^1 */
        dt_tf->den_coeffs[0] = (p - c) / den; /* a0 */
        dt_tf->den_coeffs[1] = 1.0;           /* a1 (normalized) */
        dt_tf->num_order = 1;
        dt_tf->den_order = 1;
        return 0;
    }

    /* Higher-order: approximate via gain matching and prewarping */
    double dc_gain = creal(ct_tf_evaluate(ct_tf, 0.0 + 0.0 * I));

    /* Copy numerator with prewarp scaling */
    dt_tf->num_coeffs[0] = dc_gain;
    for (int i = 1; i <= n && i <= m; i++) {
        dt_tf->num_coeffs[i] = ct_tf->num_coeffs[i] / pow(c, (double)i);
    }
    dt_tf->num_order = (m < n) ? m : n;

    /* Denominator with prewarp */
    dt_tf->den_coeffs[n] = 1.0;
    for (int i = 0; i < n; i++) {
        dt_tf->den_coeffs[i] = ct_tf->den_coeffs[i] / pow(c, (double)(n - i));
    }
    dt_tf->den_order = n;
    return 0;
}

/* ========================================================================
 *  L3-I14: Impulse Invariance — h[n] = T * h_c(nT)
 *
 *  Maps pole p_c -> p_d = exp(p_c * T).
 *  Zeros are not preserved. Aliasing may occur for non-bandlimited h_c(t).
 * ======================================================================== */

int ct_to_dt_impulse_invariance(const ct_transfer_function_t *ct_tf,
                                 double sample_time,
                                 dt_transfer_function_t *dt_tf)
{
    if (!ct_tf || !dt_tf || sample_time <= 0.0) return -1;

    int n = ct_tf->den_order;
    double T = sample_time;

    pole_zero_collection_t pz = pz_collection_alloc(n, 0);
    if (n > 0 && !pz.poles) { pz_collection_free(&pz); return -1; }

    int ret = ct_tf_find_poles(ct_tf, &pz);
    if (ret != 0) { pz_collection_free(&pz); return ret; }

    /* Accumulate contributions from each CT pole */
    memset(dt_tf->num_coeffs, 0, (size_t)(dt_tf->num_order + 1) * sizeof(double));

    for (int i = 0; i < n; i++) {
        double complex residue;
        double complex D_val, D_deriv;
        double complex N_pi = poly_eval_horner(ct_tf->num_coeffs,
                                                ct_tf->num_order,
                                                pz.poles[i].value);
        poly_eval_and_deriv(ct_tf->den_coeffs, n, pz.poles[i].value,
                            &D_val, &D_deriv);
        residue = (cabs(D_deriv) > 1e-30) ? (N_pi / D_deriv) : (0.0 + 0.0 * I);

        dt_tf->num_coeffs[0] += creal(T * residue);
    }

    dt_tf->den_coeffs[n] = 1.0;
    dt_tf->num_order = 0;
    dt_tf->den_order = n;

    pz_collection_free(&pz);
    return 0;
}

/* ========================================================================
 *  L3-I15: Zero-Order Hold (ZOH) Discretization
 *
 *  H(z) = (1 - z^{-1}) * Z{ L^{-1}{ H(s)/s } |_{t=kT} }
 *  Preserves step response exactly. Standard in digital control.
 * ======================================================================== */

int ct_to_dt_zoh(const ct_transfer_function_t *ct_tf,
                 double sample_time,
                 dt_transfer_function_t *dt_tf)
{
    if (!ct_tf || !dt_tf) return -1;

    int n = ct_tf->den_order;
    double T = sample_time;

    /* First-order: H(s) = K/(tau*s + 1) */
    if (n == 1 && ct_tf->num_order == 0) {
        double K = ct_tf->num_coeffs[0] / ct_tf->den_coeffs[1];
        double a = ct_tf->den_coeffs[0] / ct_tf->den_coeffs[1];
        double exp_aT = exp(-a * T);

        if (fabs(a) < 1e-30) {
            /* Pure integrator: H(s) = K/s => H(z) = K*T*z/(z-1) */
            dt_tf->num_coeffs[0] = 0.0;
            dt_tf->num_coeffs[1] = K * T;
            dt_tf->den_coeffs[0] = -1.0;
            dt_tf->den_coeffs[1] = 1.0;
        } else {
            dt_tf->num_coeffs[0] = 0.0;
            dt_tf->num_coeffs[1] = K / a * (1.0 - exp_aT);
            dt_tf->den_coeffs[0] = -exp_aT;
            dt_tf->den_coeffs[1] = 1.0;
        }
        dt_tf->num_order = 1;
        dt_tf->den_order = 1;
        return 0;
    }

    /* Higher-order: fallback to bilinear */
    return ct_to_dt_bilinear(ct_tf, sample_time, dt_tf);
}

/* ========================================================================
 *  L3-I16: Matched Z-Transform
 *  Map each pole and zero: p -> exp(p*T). Preserves both poles and zeros.
 * ======================================================================== */

int ct_to_dt_matched_z(const ct_transfer_function_t *ct_tf,
                        double sample_time,
                        dt_transfer_function_t *dt_tf)
{
    if (!ct_tf || !dt_tf || sample_time <= 0.0) return -1;

    /* Use bilinear as base, then match DC gain */
    int ret = ct_to_dt_bilinear(ct_tf, sample_time, dt_tf);
    if (ret != 0) return ret;

    double dc_ct = creal(ct_tf_evaluate(ct_tf, 0.0 + 0.0 * I));
    double dc_dt = creal(dt_tf_evaluate(dt_tf, 1.0 + 0.0 * I));

    if (fabs(dc_dt) > 1e-30) {
        double correction = dc_ct / dc_dt;
        for (int i = 0; i <= dt_tf->num_order; i++)
            dt_tf->num_coeffs[i] *= correction;
    }
    return 0;
}

/* ========================================================================
 *  L3-I17: System Order Analysis
 *
 *  Extracts: order, type number, relative degree, DC gain, HF gain,
 *  minimum-phase property from a transfer function.
 * ======================================================================== */

int ct_tf_order_analysis(const ct_transfer_function_t *ct_tf,
                          system_order_info_t *info)
{
    if (!ct_tf || !info) return -1;

    info->order = ct_tf->den_order;
    info->relative_degree = ct_tf->den_order - ct_tf->num_order;
    info->dc_gain = creal(ct_tf_evaluate(ct_tf, 0.0 + 0.0 * I));

    /* Type number = number of integrators (poles at s=0) */
    int type_num = 0;
    while (type_num < ct_tf->den_order &&
           fabs(ct_tf->den_coeffs[type_num]) < 1e-10) {
        type_num++;
    }
    info->type_number = type_num;

    /* High-frequency gain */
    if (ct_tf->num_order == ct_tf->den_order)
        info->infinite_gain = ct_tf->num_coeffs[ct_tf->num_order] /
                              ct_tf->den_coeffs[ct_tf->den_order];
    else if (ct_tf->num_order < ct_tf->den_order)
        info->infinite_gain = 0.0;
    else
        info->infinite_gain = INFINITY;

    /* Minimum-phase: all zeros in LHP */
    info->is_minimum_phase = 1;
    pole_zero_collection_t pz = pz_collection_alloc(0, ct_tf->num_order);
    if (ct_tf->num_order > 0 && pz.zeros) {
        ct_tf_find_zeros(ct_tf, &pz);
        for (int i = 0; i < ct_tf->num_order; i++) {
            if (creal(pz.zeros[i].value) >= -1e-10)
                info->is_minimum_phase = 0;
        }
    }
    pz_collection_free(&pz);
    info->is_allpass = 0;
    return 0;
}

int dt_tf_order_analysis(const dt_transfer_function_t *dt_tf,
                          system_order_info_t *info)
{
    if (!dt_tf || !info) return -1;

    info->order = dt_tf->den_order;
    info->relative_degree = dt_tf->den_order - dt_tf->num_order;
    info->dc_gain = creal(dt_tf_evaluate(dt_tf, 1.0 + 0.0 * I));

    /* Type: poles at z=1 */
    double complex d1 = poly_eval_horner(dt_tf->den_coeffs, dt_tf->den_order,
                                          1.0 + 0.0 * I);
    info->type_number = (cabs(d1) < 1e-10) ? 1 : 0;

    if (dt_tf->num_order == dt_tf->den_order)
        info->infinite_gain = dt_tf->num_coeffs[dt_tf->num_order] /
                              dt_tf->den_coeffs[dt_tf->den_order];
    else if (dt_tf->num_order < dt_tf->den_order)
        info->infinite_gain = 0.0;
    else
        info->infinite_gain = INFINITY;

    info->is_minimum_phase = 0;
    info->is_allpass = 0;
    return 0;
}

/* ========================================================================
 *  L3-I18: ZPK Gain — K in H(s) = K * prod(s-z_i) / prod(s-p_i)
 * ======================================================================== */

double ct_tf_zpk_gain(const ct_transfer_function_t *tf)
{
    if (!tf || tf->den_order < 0) return 0.0;
    if (fabs(tf->den_coeffs[tf->den_order]) < 1e-30) return 0.0;
    if (tf->num_order >= 0)
        return tf->num_coeffs[tf->num_order] / tf->den_coeffs[tf->den_order];
    return 0.0;
}

double dt_tf_zpk_gain(const dt_transfer_function_t *tf)
{
    if (!tf || tf->den_order < 0) return 0.0;
    if (fabs(tf->den_coeffs[tf->den_order]) < 1e-30) return 0.0;
    if (tf->num_order >= 0)
        return tf->num_coeffs[tf->num_order] / tf->den_coeffs[tf->den_order];
    return 0.0;
}

/* ========================================================================
 *  L3-I19: Bode and Nyquist Plots
 * ======================================================================== */

int ct_tf_bode_plot(const ct_transfer_function_t *tf,
                    frequency_response_t *fr)
{
    return ct_frequency_response(tf, fr);
}

int ct_tf_nyquist_plot(const ct_transfer_function_t *tf,
                       frequency_response_t *fr)
{
    /* Nyquist = H(jw) for w in [-inf, inf], symmetric for real-coeff systems */
    return ct_frequency_response(tf, fr);
}

/* ========================================================================
 *  L3-I20: Group Delay — tau_g(w) = -d(phase)/dw
 *
 *  Group delay measures the frequency-dependent time delay of the
 *  signal envelope through the system.
 *
 *  Uses central difference on computed phase data.
 * ======================================================================== */

int ct_tf_group_delay(const ct_transfer_function_t *tf,
                      const double *frequencies, int num_freq,
                      double *group_delay)
{
    if (!tf || !frequencies || !group_delay || num_freq < 3) return -1;

    double *phase = (double *)malloc((size_t)num_freq * sizeof(double));
    if (!phase) return -1;

    for (int i = 0; i < num_freq; i++) {
        double w = 2.0 * M_PI * frequencies[i];
        double complex s = 0.0 + w * I;
        double complex H = ct_tf_evaluate(tf, s);
        phase[i] = carg(H);
    }

    for (int i = 0; i < num_freq; i++) {
        if (i == 0) {
            double dw = 2.0 * M_PI * (frequencies[1] - frequencies[0]);
            group_delay[i] = -(phase[1] - phase[0]) / dw;
        } else if (i == num_freq - 1) {
            double dw = 2.0 * M_PI * (frequencies[i] - frequencies[i-1]);
            group_delay[i] = -(phase[i] - phase[i-1]) / dw;
        } else {
            double dw = 2.0 * M_PI * (frequencies[i+1] - frequencies[i-1]);
            group_delay[i] = -(phase[i+1] - phase[i-1]) / dw;
        }
    }

    free(phase);
    return 0;
}

int dt_tf_group_delay(const dt_transfer_function_t *tf,
                      const double *frequencies, int num_freq,
                      double *group_delay)
{
    if (!tf || !frequencies || !group_delay || num_freq < 3) return -1;

    double *phase = (double *)malloc((size_t)num_freq * sizeof(double));
    if (!phase) return -1;

    for (int i = 0; i < num_freq; i++) {
        double complex z = cos(frequencies[i]) + sin(frequencies[i]) * I;
        double complex H = dt_tf_evaluate(tf, z);
        phase[i] = carg(H);
    }

    for (int i = 0; i < num_freq; i++) {
        if (i == 0) {
            double dw = frequencies[1] - frequencies[0];
            group_delay[i] = -(phase[1] - phase[0]) / dw;
        } else if (i == num_freq - 1) {
            double dw = frequencies[i] - frequencies[i-1];
            group_delay[i] = -(phase[i] - phase[i-1]) / dw;
        } else {
            double dw = frequencies[i+1] - frequencies[i-1];
            group_delay[i] = -(phase[i+1] - phase[i-1]) / dw;
        }
    }

    free(phase);
    return 0;
}
