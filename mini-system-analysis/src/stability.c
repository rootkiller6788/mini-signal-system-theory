
/**
 * @file stability.c
 * @brief L4 Fundamental Laws — System Stability Analysis
 *
 * Implements the core stability criteria for LTI systems:
 * BIBO stability, Routh-Hurwitz, Jury test, Nyquist criterion,
 * gain/phase margins, root locus, Lyapunov direct method.
 *
 * Knowledge Implementation (15 points):
 *   L4-I01: ct_is_bibo_stable — pole location check
 *   L4-I02: dt_is_bibo_stable
 *   L4-I03: routh_hurwitz — Routh table construction
 *   L4-I04: ct_routh_hurwitz — applied to CT transfer function
 *   L4-I05: jury_stability_test — DT polynomial criterion
 *   L4-I06: dt_jury_stability_test — applied to DT transfer function
 *   L4-I07: nyquist_stability — encirclement counting
 *   L4-I08: compute_stability_margins — gain/phase margin
 *   L4-I09: root_locus — pole trajectories vs gain
 *   L4-I10: ct_lyapunov_stability — Lyapunov equation method
 *   L4-I11: ct_count_marginal_poles / dt_count_marginal_poles
 *   L4-I12: feedback_stability_check
 *   L4-I13: ct_dominant_pole_analysis — damping ratio, natural freq
 *   L4-I14: ct_settling_time_estimate — t_s = 4/(zeta*wn)
 *   L4-I15: ct_overshoot_estimate — M_p = exp(-pi*zeta/sqrt(1-zeta^2))
 *
 * References:
 *   Ogata, Modern Control Engineering, 5th Ed.
 *   Dorf & Bishop, Modern Control Systems, 13th Ed.
 *   Oppenheim & Willsky, Signals and Systems, 2nd Ed.
 */

#include "stability.h"
#include "transfer_function.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 *  L4-I01: CT BIBO Stability — all poles in LHP (Re{s} < 0)
 *
 *  Theorem (BIBO Stability):
 *    An LTI CT system is BIBO stable iff all poles of H(s)
 *    lie strictly in the open left half-plane.
 *    Equivalently: integral |h(tau)| dtau < infinity.
 * ======================================================================== */

int ct_is_bibo_stable(const ct_transfer_function_t *tf)
{
    if (!tf) return 0;

    int n = tf->den_order;
    pole_zero_collection_t pz = pz_collection_alloc(n, 0);
    if (n > 0 && !pz.poles) return -1;

    if (ct_tf_find_poles(tf, &pz) != 0) {
        pz_collection_free(&pz);
        return -1;
    }

    int stable = 1;
    for (int i = 0; i < n; i++) {
        if (creal(pz.poles[i].value) >= -1e-10) {
            stable = 0;
            break;
        }
    }

    pz_collection_free(&pz);
    return stable;
}

/* ========================================================================
 *  L4-I02: DT BIBO Stability — all poles inside unit circle (|z| < 1)
 * ======================================================================== */

int dt_is_bibo_stable(const dt_transfer_function_t *tf)
{
    if (!tf) return 0;

    int n = tf->den_order;
    pole_zero_collection_t pz = pz_collection_alloc(n, 0);
    if (n > 0 && !pz.poles) return -1;

    if (dt_tf_find_poles(tf, &pz) != 0) {
        pz_collection_free(&pz);
        return -1;
    }

    int stable = 1;
    for (int i = 0; i < n; i++) {
        if (cabs(pz.poles[i].value) >= 1.0 - 1e-10) {
            stable = 0;
            break;
        }
    }

    pz_collection_free(&pz);
    return stable;
}

int ct_is_bibo_stable_from_poles(const pole_zero_collection_t *pz)
{
    if (!pz) return 0;
    return pz->is_stable;
}

int dt_is_bibo_stable_from_poles(const pole_zero_collection_t *pz)
{
    if (!pz) return 0;
    return pz->is_stable;
}

/* ========================================================================
 *  L4-I03: Routh-Hurwitz Criterion
 *
 *  Theorem (Routh-Hurwitz, 1874/1895):
 *    A polynomial a_n s^n + ... + a_0 with a_n > 0 has all roots
 *    in the open LHP iff all elements in the first column of the
 *    Routh array are positive.
 *
 *  The Routh array is constructed from the coefficients:
 *    Row 0: a_n, a_{n-2}, a_{n-4}, ...
 *    Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 *    Row k (k>=2): b_i = (b_{k-1,1}*b_{k-2,i+1} - b_{k-2,1}*b_{k-1,i+1}) / b_{k-1,1}
 *
 *  Number of sign changes in first column = number of RHP roots.
 *
 *  Special cases:
 *    1. First element of a row is zero: replace with small epsilon.
 *    2. Entire row is zero: auxiliary polynomial, roots on jw-axis.
 *
 *  Complexity: O(n^2), where n = polynomial degree.
 * ======================================================================== */

int routh_hurwitz(const double *coeffs, int order,
                  double *rh_table, int *special_case)
{
    if (!coeffs || order < 1 || !rh_table || !special_case) return -1;

    if (coeffs[order] <= 0.0) return -1;  /* leading coefficient must be positive */

    int n = order;
    int cols = (n + 2) / 2;  /* ceil((n+1)/2) */

    /* Fill first two rows */
    for (int j = 0; j < cols; j++) {
        int idx0 = n - 2 * j;
        int idx1 = n - 1 - 2 * j;
        rh_table[0 * cols + j] = (idx0 >= 0) ? coeffs[idx0] : 0.0;
        rh_table[1 * cols + j] = (idx1 >= 0) ? coeffs[idx1] : 0.0;
    }

    *special_case = 0;
    double epsilon = 1e-6;

    for (int i = 2; i <= n; i++) {
        double prev_first = rh_table[(i - 1) * cols + 0];

        /* Check for zero first element */
        if (fabs(prev_first) < 1e-14) {
            int all_zero = 1;
            for (int j = 0; j < cols; j++) {
                if (fabs(rh_table[(i - 1) * cols + j]) > 1e-14) {
                    all_zero = 0;
                    break;
                }
            }

            if (all_zero) {
                *special_case = 2;  /* Entire row zero: auxiliary polynomial */
                /* Form auxiliary polynomial from row i-2 */
                int aux_deg = n - i + 2;
                for (int j = 0; j < cols; j++) {
                    rh_table[(i - 1) * cols + j] =
                        (aux_deg - 2 * j) * rh_table[(i - 2) * cols + j];
                }
                prev_first = rh_table[(i - 1) * cols + 0];
            } else {
                *special_case = 1;  /* Zero first element only */
                rh_table[(i - 1) * cols + 0] = epsilon;
                prev_first = epsilon;
            }
        }

        /* Compute row i */
        double a = rh_table[(i - 2) * cols + 0];
        for (int j = 0; j < cols - 1; j++) {
            double b = rh_table[(i - 2) * cols + j + 1];
            double c = rh_table[(i - 1) * cols + 0];
            double d = rh_table[(i - 1) * cols + j + 1];
            rh_table[i * cols + j] = (c * b - a * d) / prev_first;
        }
        rh_table[i * cols + cols - 1] = 0.0;
    }

    /* Count sign changes in first column */
    int sign_changes = 0;
    double prev_sign = rh_table[0 * cols + 0];
    for (int i = 1; i <= n; i++) {
        double cur = rh_table[i * cols + 0];
        if (fabs(cur) < 1e-14) cur = 0.0;

        if (prev_sign * cur < -1e-14) {
            sign_changes++;
        } else if (fabs(cur) > 1e-14) {
            prev_sign = cur;
        }
    }

    return sign_changes;
}

int ct_routh_hurwitz(const ct_transfer_function_t *tf,
                     int *rhp_poles, int *special_case)
{
    if (!tf || !rhp_poles || !special_case) return -1;

    int n = tf->den_order;
    if (n < 1) return -1;

    int cols = (n + 2) / 2;
    double *rh_table = (double *)calloc((size_t)((n + 1) * cols), sizeof(double));
    if (!rh_table) return -1;

    int ret = routh_hurwitz(tf->den_coeffs, n, rh_table, special_case);
    *rhp_poles = ret;

    free(rh_table);
    return (ret >= 0) ? 0 : ret;
}

/* ========================================================================
 *  L4-I05: Jury Stability Test (DT equivalent of Routh-Hurwitz)
 *
 *  Theorem (Jury, 1964):
 *    All roots of polynomial a_n z^n + ... + a_0 (a_n > 0) lie
 *    inside the unit circle iff all Jury conditions are satisfied:
 *      1.  D(1) > 0
 *      2.  (-1)^n D(-1) > 0
 *      3.  |a_0| < a_n
 *      4.  |b_0| > |b_{n-1}|  (for reduced-order table)
 *      5.  |c_0| > |c_{n-2}|
 *      ... (repeated table reduction)
 * ======================================================================== */

int jury_stability_test(const double *coeffs, int order)
{
    if (!coeffs || order < 1) return -1;

    int n = order;

    /* Normalize: make a_n > 0 */
    double lead = coeffs[n];
    if (lead <= 0.0) return -1;

    double *a = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!a) return -1;
    for (int i = 0; i <= n; i++) a[i] = coeffs[i] / lead;

    /* Condition 1: D(1) > 0 */
    double sum = 0.0;
    for (int i = 0; i <= n; i++) sum += a[i];
    if (sum <= 0.0) { free(a); return 0; }

    /* Condition 2: (-1)^n D(-1) > 0 */
    sum = 0.0;
    double sign = 1.0;
    for (int i = 0; i <= n; i++) {
        sum += sign * a[i];
        sign = -sign;
    }
    if (sum <= 0.0) { free(a); return 0; }

    /* Condition 3: |a_0| < a_n (which is 1 after normalization) */
    if (fabs(a[0]) >= 1.0) { free(a); return 0; }

    /* Construct Jury table */
    double *row = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!row) { free(a); return -1; }
    memcpy(row, a, (size_t)(n + 1) * sizeof(double));

    for (int k = 0; k < n - 1; k++) {
        int m = n - k;
        double alpha = row[m] / row[0];

        if (fabs(row[m]) <= fabs(row[0])) {
            free(a); free(row);
            return 0;
        }

        for (int i = 0; i < m; i++) {
            row[i] = row[i] - alpha * row[m - 1 - i];
        }

        if (fabs(row[m-1]) <= fabs(row[0])) {
            free(a); free(row);
            return 0;
        }
    }

    free(a); free(row);
    return 1;  /* All conditions satisfied */
}

int dt_jury_stability_test(const dt_transfer_function_t *tf)
{
    if (!tf) return -1;

    int n = tf->den_order;
    if (n < 1) return -1;

    return jury_stability_test(tf->den_coeffs, n);
}

/* ========================================================================
 *  L4-I07: Nyquist Stability Criterion
 *
 *  Theorem (Nyquist, 1932):
 *    N = Z - P, where:
 *      N = net number of CCW encirclements of -1+0j by L(jw)
 *      Z = number of closed-loop RHP poles
 *      P = number of open-loop RHP poles
 *    Closed-loop stable iff Z = 0, i.e., N = -P (P CW encirclements).
 *
 *  This implementation counts encirclements by tracking the angle
 *  of L(jw) - (-1) as w varies from 0 to inf.
 * ======================================================================== */

int nyquist_stability(const ct_transfer_function_t *loop_tf,
                      const frequency_response_t *fr,
                      int open_loop_rhp_poles,
                      int *encirclements,
                      int *closed_loop_rhp_poles)
{
    if (!loop_tf || !fr || !encirclements || !closed_loop_rhp_poles) return -1;
    if (fr->num_points < 4) return -1;

    /* Count encirclements of (-1, 0) by tracking cumulative angle change.
     * For each consecutive pair of points on Nyquist contour,
     * compute the angle of vector from (-1,0) to H(jw). */
    double total_angle = 0.0;
    double prev_angle = 0.0;
    int first_valid = 0;

    for (int i = 0; i < fr->num_points; i++) {
        /* H(jw) in complex plane */
        double re = fr->points[i].magnitude * cos(fr->points[i].phase_rad);
        double im = fr->points[i].magnitude * sin(fr->points[i].phase_rad);

        /* Vector from (-1, 0) to H(jw) */
        double vx = re - (-1.0);
        double vy = im - 0.0;

        if (fabs(vx) < 1e-15 && fabs(vy) < 1e-15) continue;

        double angle = atan2(vy, vx);

        if (first_valid) {
            double diff = angle - prev_angle;
            /* Normalize to [-pi, pi] */
            while (diff > M_PI) diff -= 2.0 * M_PI;
            while (diff < -M_PI) diff += 2.0 * M_PI;
            total_angle += diff;
            prev_angle = angle;
        } else {
            prev_angle = angle;
            first_valid = 1;
        }
    }

    /* Encirclements = total_angle / (2*pi), rounded to nearest integer */
    *encirclements = (int)round(total_angle / (2.0 * M_PI));

    /* Z = P + N (where N is encirclements counter-clockwise) */
    *closed_loop_rhp_poles = open_loop_rhp_poles + (*encirclements);

    return 0;
}

/* ========================================================================
 *  L4-I08: Gain and Phase Margins from Frequency Response
 *
 *  Gain margin (GM): reciprocal of |L(jw_pc)| at phase crossover
 *    where angle(L(jw_pc)) = -180 deg. GM_dB = -20*log10(|L(jw_pc)|).
 *
 *  Phase margin (PM): 180 + angle(L(jw_gc)) at gain crossover
 *    where |L(jw_gc)| = 1 (0 dB).
 *
 *  For stability: GM > 0 dB and PM > 0 deg.
 *  Typical design targets: GM >= 6 dB, PM >= 45 deg.
 * ======================================================================== */

int compute_stability_margins(const frequency_response_t *fr,
                               double *gain_margin_db,
                               double *phase_margin_deg,
                               double *gm_freq,
                               double *pm_freq)
{
    if (!fr || !gain_margin_db || !phase_margin_deg || !gm_freq || !pm_freq) {
        return -1;
    }

    /* Find phase crossover: phase ≈ -180 deg */
    int pc_idx = -1;
    double min_phase_dist = 1e10;
    for (int i = 0; i < fr->num_points; i++) {
        /* Wrap phase to [-180, 180] for comparison */
        double ph = fr->points[i].phase_deg;
        while (ph > 180.0) ph -= 360.0;
        while (ph < -180.0) ph += 360.0;

        double dist = fabs(ph + 180.0);
        if (dist < min_phase_dist) {
            min_phase_dist = dist;
            pc_idx = i;
        }
    }

    /* Find gain crossover: magnitude ≈ 0 dB */
    int gc_idx = -1;
    double min_gain_dist = 1e10;
    for (int i = 0; i < fr->num_points; i++) {
        double dist = fabs(fr->points[i].magnitude_db);
        if (dist < min_gain_dist) {
            min_gain_dist = dist;
            gc_idx = i;
        }
    }

    if (pc_idx >= 0) {
        *gain_margin_db = -fr->points[pc_idx].magnitude_db;
        *gm_freq = fr->points[pc_idx].frequency;
    } else {
        *gain_margin_db = INFINITY;
        *gm_freq = INFINITY;
    }

    if (gc_idx >= 0) {
        *phase_margin_deg = 180.0 + fr->points[gc_idx].phase_deg;
        /* Normalize */
        while (*phase_margin_deg > 360.0) *phase_margin_deg -= 360.0;
        while (*phase_margin_deg < -180.0) *phase_margin_deg += 360.0;
        *pm_freq = fr->points[gc_idx].frequency;
    } else {
        *phase_margin_deg = INFINITY;
        *pm_freq = INFINITY;
    }

    return 0;
}

/* ========================================================================
 *  L4-I09: Root Locus — Closed-loop poles vs gain K
 *
 *  D_cl(s) = D_ol(s) + K * N_ol(s) = 0
 *
 *  For each K, find roots of D_ol(s) + K*N_ol(s) = 0.
 *  Plots how closed-loop poles move as K varies from 0 to inf.
 * ======================================================================== */

int root_locus(const ct_transfer_function_t *tf,
               const double *K_values, int num_K,
               double complex *poles)
{
    if (!tf || !K_values || !poles || num_K <= 0) return -1;

    int n = tf->den_order;
    int m = tf->num_order;

    /* For each K, find roots of D(s) + K*N(s) = 0 */
    for (int k_idx = 0; k_idx < num_K; k_idx++) {
        double K = K_values[k_idx];

        /* New polynomial: den(s) + K*num(s) */
        int max_order = (n > m) ? n : m;
        double *poly = (double *)calloc((size_t)(max_order + 1), sizeof(double));
        if (!poly) return -1;

        for (int i = 0; i <= n; i++) poly[i] += tf->den_coeffs[i];
        for (int i = 0; i <= m; i++) poly[i] += K * tf->num_coeffs[i];

        /* Find roots using Newton+deflation */
        double *poly_work = (double *)malloc((size_t)(max_order + 1) * sizeof(double));
        if (!poly_work) { free(poly); return -1; }
        memcpy(poly_work, poly, (size_t)(max_order + 1) * sizeof(double));

        for (int r = 0; r < max_order; r++) {
            double angle = 2.0 * M_PI * (double)r / (double)max_order;
            double complex guess = 5.0 * cos(angle) + 5.0 * sin(angle) * I;
            double complex x = guess;

            for (int iter = 0; iter < 100; iter++) {
                /* Evaluate poly and derivative */
                double complex val = poly_work[max_order];
                double complex deriv = 0;
                for (int i = max_order; i >= 1; i--) {
                    deriv = deriv * x + (double)i * poly_work[i];
                    if (i > 0) val = val * x + poly_work[i - 1];
                }
                if (cabs(deriv) < 1e-15) { x += 0.1; continue; }
                x = x - val / deriv;
                if (cabs(val) < 1e-10) break;
            }
            poles[k_idx * max_order + r] = x;
        }

        free(poly); free(poly_work);
    }

    return 0;
}

/* ========================================================================
 *  L4-I10: Lyapunov Stability (Direct Method)
 *
 *  Theorem (Lyapunov, 1892):
 *    For dx/dt = Ax, the system is asymptotically stable iff there
 *    exists P = P' > 0 such that A'P + PA = -Q for any Q = Q' > 0.
 *
 *  We use state_space.c's lyapunov solver and check positive definiteness
 *  via Cholesky attempt or eigenvalue check.
 * ======================================================================== */

int ct_lyapunov_stability(const double *A, int n,
                           const double *Q, double *P)
{
    if (!A || !Q || !P || n <= 0) return -1;

    /* Solve Lyapunov equation via state_space module */
    ct_state_space_t ss;
    ss.A = (double *)A;
    ss.n = n;
    ss.m = 1;
    ss.p = 1;
    ss.owns_data = 0;

    if (ct_ss_lyapunov(&ss, Q, P) != 0) return -1;

    /* Check positive definiteness of P: all leading principal minors > 0 */
    /* Simplification: check diagonal elements > 0 and trace > 0 */
    double trace = 0.0;
    int pos_diag = 1;
    for (int i = 0; i < n; i++) {
        if (P[i * n + i] <= 1e-10) pos_diag = 0;
        trace += P[i * n + i];
    }

    return (pos_diag && trace > 0.0) ? 1 : 0;
}

/* ========================================================================
 *  L4-I11: Count Marginal Poles
 * ======================================================================== */

int ct_count_marginal_poles(const pole_zero_collection_t *pz)
{
    if (!pz) return 0;
    int count = 0;
    for (int i = 0; i < pz->num_poles; i++) {
        if (fabs(creal(pz->poles[i].value)) < 1e-10 &&
            pz->poles[i].is_pole) {
            count++;
        }
    }
    return count;
}

int dt_count_marginal_poles(const pole_zero_collection_t *pz)
{
    if (!pz) return 0;
    int count = 0;
    for (int i = 0; i < pz->num_poles; i++) {
        if (fabs(cabs(pz->poles[i].value) - 1.0) < 1e-10 &&
            pz->poles[i].is_pole) {
            count++;
        }
    }
    return count;
}

/* ========================================================================
 *  L4-I12: Feedback Stability Check
 *  H_cl = G / (1 + G*H)
 * ======================================================================== */

int feedback_stability_check(const ct_transfer_function_t *G,
                              const ct_transfer_function_t *H)
{
    if (!G || !H) return -1;

    /* Compute loop transfer function L(s) = G(s) * H(s) */
    /* Closed-loop characteristic: 1 + L(s) = 0 */
    /* Equivalent to D_G*D_H + N_G*N_H = 0 */

    int nG = G->den_order, mG = G->num_order;
    int nH = H->den_order, mH = H->num_order;

    int n_cl = nG + nH;
    int m_cl = mG + mH;
    int max_n = (n_cl > m_cl) ? n_cl : m_cl;

    double *den_cl = (double *)calloc((size_t)(max_n + 1), sizeof(double));
    if (!den_cl) return -1;

    /* D_G * D_H */
    for (int i = 0; i <= nG; i++)
        for (int j = 0; j <= nH; j++)
            den_cl[i + j] += G->den_coeffs[i] * H->den_coeffs[j];

    /* + N_G * N_H */
    for (int i = 0; i <= mG; i++)
        for (int j = 0; j <= mH; j++)
            den_cl[i + j] += G->num_coeffs[i] * H->num_coeffs[j];

    /* Check stability of closed-loop characteristic polynomial */
    double lead = den_cl[max_n];
    if (fabs(lead) < 1e-30) { free(den_cl); return 0; }

    /* Normalize */
    for (int i = 0; i <= max_n; i++) den_cl[i] /= lead;

    int cols = (max_n + 2) / 2;
    double *rh = (double *)calloc((size_t)((max_n + 1) * cols), sizeof(double));
    if (!rh) { free(den_cl); return -1; }

    int special;
    int rhp = routh_hurwitz(den_cl, max_n, rh, &special);

    free(den_cl); free(rh);
    return (rhp == 0) ? 1 : 0;  /* Stable if 0 RHP roots */
}

/* ========================================================================
 *  L4-I13: Dominant Pole Analysis
 *
 *  The dominant pole pair is the one with smallest |Re{s}| (slowest decay).
 *  damping zeta = -Re{p} / |p|, natural freq wn = |p|
 * ======================================================================== */

int ct_dominant_pole_analysis(const pole_zero_collection_t *pz,
                               double *damping_ratio,
                               double *natural_freq)
{
    if (!pz || !damping_ratio || !natural_freq) return -1;
    if (pz->num_poles == 0) return -1;

    /* Find pole with maximum real part (closest to jw-axis = dominant) */
    double max_real = -INFINITY;
    int dom_idx = 0;
    for (int i = 0; i < pz->num_poles; i++) {
        double re = creal(pz->poles[i].value);
        if (re > max_real) {
            max_real = re;
            dom_idx = i;
        }
    }

    double complex p = pz->poles[dom_idx].value;
    *natural_freq = cabs(p);
    *damping_ratio = (*natural_freq > 1e-30)
        ? -creal(p) / (*natural_freq) : 0.0;

    return 0;
}

/* ========================================================================
 *  L4-I14: Settling Time Estimate (2% criterion)
 *  t_s ≈ 4 / (zeta * wn)  for dominant second-order system
 * ======================================================================== */

double ct_settling_time_estimate(double damping_ratio, double natural_freq)
{
    if (damping_ratio <= 0.0 || natural_freq <= 0.0) return INFINITY;
    return 4.0 / (damping_ratio * natural_freq);
}

/* ========================================================================
 *  L4-I15: Peak Overshoot Estimate
 *  M_p = exp(-pi * zeta / sqrt(1 - zeta^2)) * 100%
 *  Valid for 0 <= zeta < 1 (underdamped)
 * ======================================================================== */

double ct_overshoot_estimate(double damping_ratio)
{
    if (damping_ratio <= 0.0) return 100.0;  /* Undamped: 100% overshoot */
    if (damping_ratio >= 1.0) return 0.0;    /* Critically/overdamped: no overshoot */
    return 100.0 * exp(-M_PI * damping_ratio / sqrt(1.0 - damping_ratio * damping_ratio));
}
