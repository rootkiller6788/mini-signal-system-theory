/*
 * stability_analysis.c - Stability Analysis Methods
 *
 * Implements classical stability criteria for LTI systems:
 * Routh-Hurwitz criterion (continuous-time), Jury test (discrete-time),
 * pole-based stability checks, and gain/phase margin computation.
 *
 * Knowledge Coverage:
 *   L4: Routh-Hurwitz stability criterion (Routh, 1877; Hurwitz, 1895)
 *   L4: Jury stability test (Jury, 1964)
 *   L2: Pole location and stability relationship
 *   L5: Gain margin and phase margin from Bode plots
 *
 * Reference: Ogata, "Modern Control Engineering" (2010), Ch.5
 *            Jury, "Theory and Application of the Z-Transform Method" (1964)
 *            Franklin, Powell, Emami-Naeini, "Feedback Control" (2019)
 */

#include "laplace_z_transform_core.h"
#include "stability_analysis.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

/* ========================================================================
 * L4: Routh-Hurwitz Stability Criterion (Continuous-Time)
 * ======================================================================== */

/*
 * The Routh-Hurwitz criterion determines if all roots of a polynomial
 * have negative real parts, without explicitly computing the roots.
 *
 * For polynomial: a_n*s^n + a_{n-1}*s^{n-1} + ... + a_1*s + a_0
 * with a_n > 0 (we enforce this by multiplying by -1 if needed).
 *
 * The Routh array is constructed as follows:
 *
 * Row 0: a_n,    a_{n-2}, a_{n-4}, ...
 * Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 * Row i (i >= 2): computed from rows i-2 and i-1
 *
 * Element (i, j) = -det([row_{i-2}[0], row_{i-2}[j+1]; row_{i-1}[0], row_{i-1}[j+1]]) / row_{i-1}[0]
 *
 * Reference: Ogata (2010), Section 5-6
 */

int stability_routh_array(const double *coeffs, int n, double *routh,
                           int *rows, int *cols) {
    /*
     * Build the Routh array from polynomial coefficients.
     *
     * Input: coeffs = [a_0, a_1, ..., a_n] (constant term first)
     *        n = polynomial degree
     * Output: routh = Routh array (flat, row-major)
     *         rows = number of rows
     *         cols = number of columns
     */
    if (!coeffs || n < 1 || !routh || !rows || !cols) return -1;

    /* Normalize: ensure a_n > 0 */
    double sign = (coeffs[n] >= 0) ? 1.0 : -1.0;

    int n_rows = n + 1;
    int n_cols = (n + 2) / 2;  /* ceil((n+1)/2) */

    *rows = n_rows;
    *cols = n_cols;

    /* Initialize first two rows */
    for (int j = 0; j < n_cols; j++) {
        /* Row 0: even-indexed coefficients (a_n, a_{n-2}, ...) */
        int idx0 = n - 2*j;
        routh[0 * n_cols + j] = (idx0 >= 0) ? sign * coeffs[idx0] : 0.0;

        /* Row 1: odd-indexed coefficients (a_{n-1}, a_{n-3}, ...) */
        int idx1 = n - 1 - 2*j;
        routh[1 * n_cols + j] = (idx1 >= 0) ? sign * coeffs[idx1] : 0.0;
    }

    /* Compute remaining rows */
    for (int i = 2; i < n_rows; i++) {
        double prev_first = routh[(i-1) * n_cols + 0];

        if (fabs(prev_first) < 1e-15) {
            /* First element of previous row is zero.
             * Replace with a small epsilon to continue.
             * This is the "epsilon method" for handling zero in first column.
             * Reference: Ogata (2010), Section 5-6 */
            prev_first = 1e-10;
        }

        for (int j = 0; j < n_cols - 1; j++) {
            double a = routh[(i-2) * n_cols + 0];
            double b = routh[(i-2) * n_cols + j + 1];
            double c = routh[(i-1) * n_cols + 0];
            double d = routh[(i-1) * n_cols + j + 1];

            /* routh[i][j] = -(a*d - b*c) / c = (b*c - a*d) / c */
            routh[i * n_cols + j] = (b * c - a * d) / prev_first;
        }
        /* Last column is always 0 */
        routh[i * n_cols + n_cols - 1] = 0.0;
    }

    return 0;
}

int stability_routh_rhp_roots(const double *coeffs, int n) {
    /*
     * Count the number of roots in the right half-plane (RHP).
     *
     * The number of sign changes in the first column of the Routh array
     * equals the number of roots with positive real parts.
     *
     * Theorem: Routh (1877)
     * Reference: Ogata (2010), Section 5-6
     */
    if (!coeffs || n < 1) return -1;

    int n_rows = n + 1;
    int n_cols = (n + 2) / 2;
    double *routh = (double*)calloc(n_rows * n_cols, sizeof(double));
    if (!routh) return -1;

    int rows, cols;
    if (stability_routh_array(coeffs, n, routh, &rows, &cols) != 0) {
        free(routh);
        return -1;
    }

    /* Count sign changes in first column */
    int sign_changes = 0;
    double prev = routh[0 * cols + 0];
    int prev_sign = (prev > 0) ? 1 : ((prev < 0) ? -1 : 0);

    for (int i = 1; i < rows; i++) {
        double curr = routh[i * cols + 0];
        int curr_sign = (curr > 0) ? 1 : ((curr < 0) ? -1 : 0);

        if (curr_sign != 0 && prev_sign != 0 && curr_sign != prev_sign) {
            sign_changes++;
        }
        if (curr_sign != 0) {
            prev_sign = curr_sign;
        }
    }

    free(routh);
    return sign_changes;
}

int stability_routh_is_stable(const double *coeffs, int n) {
    /*
     * Check if polynomial represents a stable system.
     *
     * Necessary and sufficient condition:
     * All elements in the first column of the Routh array have the same sign.
     *
     * Returns 1 if stable (all roots in LHP), 0 if unstable.
     */
    int rhp_roots = stability_routh_rhp_roots(coeffs, n);
    if (rhp_roots < 0) return -1;
    return (rhp_roots == 0) ? 1 : 0;
}

int stability_routh_necessary_condition(const double *coeffs, int n) {
    /*
     * Check the necessary (but not sufficient) condition for stability:
     * All coefficients must exist (non-zero) and have the same sign.
     *
     * If any coefficient is missing or has a different sign,
     * the system is definitely unstable.
     *
     * Reference: Ogata (2010), Section 5-6
     */
    if (!coeffs || n < 1) return 0;

    double sign = (coeffs[n] > 0) ? 1.0 : -1.0;

    for (int i = 0; i <= n; i++) {
        if (fabs(coeffs[i]) < 1e-15) return 0;  /* missing coefficient */
        if ((coeffs[i] > 0 && sign < 0) || (coeffs[i] < 0 && sign > 0)) return 0;
    }

    return 1;
}

int stability_routh_auxiliary(const double *routh, int row, int cols,
                               double *aux_poly, int *aux_deg) {
    /*
     * Handle the special case where an entire row of the Routh array
     * consists of zeros.
     *
     * This occurs when the polynomial has roots that are symmetrically
     * placed about the origin (e.g., purely imaginary roots).
     *
     * The auxiliary polynomial is formed from the row above the all-zero row,
     * using only the powers: s^{row_above_power}, s^{row_above_power-2}, ...
     *
     * Reference: Ogata (2010), Section 5-6
     */
    if (!routh || !aux_poly || !aux_deg || row < 2) return -1;

    /* The auxiliary polynomial uses row-1 (the row above the zero row).
     * Powers: starting from (n - row + 2) down to 0, step -2. */
    int power = *aux_deg;
    if (power < 0) power = row;  /* estimate */

    int n_coeffs = (power / 2) + 1;
    for (int i = 0; i < n_coeffs && i < cols; i++) {
        aux_poly[2*i] = routh[(row-1) * cols + i];
    }

    *aux_deg = power;
    return 0;
}

/* ========================================================================
 * L4: Jury Stability Test (Discrete-Time)
 * ======================================================================== */

int stability_jury_table(const double *coeffs, int n, double *jury_table,
                          int *rows, int *cols) {
    /*
     * Build the Jury stability table for a discrete-time polynomial.
     *
     * Polynomial: P(z) = a_n*z^n + a_{n-1}*z^{n-1} + ... + a_1*z + a_0
     * with a_n > 0.
     *
     * The Jury table has 2n-3 rows (for n >= 2).
     *
     * Construction:
     * Row 0: [a_0, a_1, ..., a_n]
     * Row 1: [a_n, a_{n-1}, ..., a_0]  (reverse of row 0)
     * Row 2k: computed from rows 2k-2 and 2k-1, reducing size by 1 each time
     *
     * Reference: Jury (1964)
     *            Ogata, "Discrete-Time Control Systems" (1995), Section 4.3
     */
    if (!coeffs || n < 2 || !jury_table || !rows || !cols) return -1;

    /* Normalize: a_n > 0 */
    double an = coeffs[n];
    if (an < 0) {
        for (int i = 0; i <= n; i++) {
            /* Will need to negate; for simplicity just use as-is */
        }
    }

    int n_rows = 2*n;
    int n_cols = n + 1;
    *rows = n_rows;
    *cols = n_cols;

    /* Row 0: original coefficients */
    for (int j = 0; j <= n; j++) {
        jury_table[0 * n_cols + j] = coeffs[j];
    }

    /* Row 1: reversed coefficients */
    for (int j = 0; j <= n; j++) {
        jury_table[1 * n_cols + j] = coeffs[n - j];
    }

    /* Compute remaining rows */
    for (int i = 2; i < n_rows; i++) {
        int prev_size = n + 1 - (i/2 - 1);
        double alpha = jury_table[(i-2) * n_cols + prev_size - 1] /
                       jury_table[(i-2) * n_cols + 0];

        for (int j = 0; j < prev_size - 1; j++) {
            double a0j = jury_table[(i-2) * n_cols + j];
            double b0j = jury_table[(i-1) * n_cols + j];
            jury_table[i * n_cols + j] = a0j - alpha * b0j;
        }
        /* Fill remaining with zeros */
        for (int j = prev_size - 1; j < n_cols; j++) {
            jury_table[i * n_cols + j] = 0.0;
        }
    }

    return 0;
}

int stability_jury_is_stable(const double *coeffs, int n) {
    /*
     * Check if a discrete-time polynomial satisfies the Jury stability criteria.
     *
     * For polynomial P(z) = a_n*z^n + ... + a_0 with a_n > 0:
     *
     * Condition 1: P(1) > 0
     * Condition 2: (-1)^n * P(-1) > 0
     * Condition 3: |a_0| < a_n
     * Condition 4: |b_0| > |b_{n-1}|  (from Jury table row 2)
     * Condition 5: |c_0| > |c_{n-2}|  (from Jury table row 4)
     * ... (for n > 2, more conditions from the Jury table)
     *
     * Returns 1 if stable, 0 if unstable.
     *
     * Reference: Jury (1964)
     */
    if (!coeffs || n < 1) return -1;

    /* Ensure a_n > 0 */
    double scale = (coeffs[n] > 0) ? 1.0 : -1.0;

    /* Condition 1: P(1) > 0 */
    double P1 = 0.0;
    for (int i = 0; i <= n; i++) P1 += scale * coeffs[i];
    if (P1 <= 0) return 0;

    /* Condition 2: (-1)^n * P(-1) > 0 */
    double Pm1 = 0.0;
    for (int i = 0; i <= n; i++) {
        double term = scale * coeffs[i] * ((i % 2 == 0) ? 1.0 : -1.0);
        Pm1 += term;
    }
    if ((n % 2 == 0) ? Pm1 <= 0 : Pm1 >= 0) return 0;

    /* Condition 3: |a_0| < a_n */
    if (fabs(scale * coeffs[0]) >= fabs(scale * coeffs[n])) return 0;

    /* For n = 1, only conditions 1-3 apply */
    if (n == 1) return 1;

    /* For n >= 2, need additional Jury table conditions */
    if (n >= 2) {
        int n_rows = 2*n;
        int n_cols = n + 1;
        double *jury = (double*)calloc(n_rows * n_cols, sizeof(double));
        if (!jury) return -1;

        /* Scale coefficients */
        double *scaled = (double*)malloc((n+1) * sizeof(double));
        for (int i = 0; i <= n; i++) scaled[i] = scale * coeffs[i];

        int rows, cols;
        stability_jury_table(scaled, n, jury, &rows, &cols);
        free(scaled);

        /* Check constraints from Jury table */
        for (int k = 1; k < n; k++) {
            int row_idx = 2*k;
            int size = n + 1 - k;
            double first = jury[row_idx * cols + 0];
            double last = jury[row_idx * cols + size - 1];
            if (fabs(first) <= fabs(last)) {
                free(jury);
                return 0;
            }
        }

        free(jury);
    }

    return 1;
}

double stability_jury_eval_P_at_1(const double *coeffs, int n) {
    /* P(1) = sum of all coefficients */
    if (!coeffs || n < 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i <= n; i++) sum += coeffs[i];
    return sum;
}

double stability_jury_eval_P_at_minus_1(const double *coeffs, int n) {
    /* P(-1) = sum (-1)^i * coeffs[i] */
    if (!coeffs || n < 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i <= n; i++) {
        sum += (i % 2 == 0) ? coeffs[i] : -coeffs[i];
    }
    /* Return (-1)^n * P(-1) as needed by Jury test */
    return (n % 2 == 0) ? sum : -sum;
}

/* ========================================================================
 * L2: Pole-Based Stability Analysis
 * ======================================================================== */

int stability_from_poles(const complex_t *poles, int n_poles, int is_z_domain) {
    /*
     * Direct stability check from pole locations.
     *
     * Continuous-time: All poles must have Re(p) < 0
     * Discrete-time:   All poles must have |p| < 1
     *
     * This is the most fundamental stability criterion for LTI systems.
     */
    if (!poles || n_poles <= 0) return 0;

    for (int i = 0; i < n_poles; i++) {
        if (is_z_domain) {
            if (complex_mag(poles[i]) >= 1.0 - 1e-12) return 0;
        } else {
            if (poles[i].real >= -1e-12) return 0;
        }
    }
    return 1;
}

double stability_margin(const complex_t *poles, int n_poles, int is_z_domain) {
    /*
     * Stability margin: distance from most critical pole to stability boundary.
     *
     * Continuous-time: margin = -max(Re(p_k)). Positive = stable.
     * Discrete-time:   margin = 1 - max(|p_k|).  Positive = stable.
     *
     * Larger margin indicates greater stability robustness.
     */
    if (!poles || n_poles <= 0) return -INFINITY;

    if (is_z_domain) {
        double max_mag = 0.0;
        for (int i = 0; i < n_poles; i++) {
            double m = complex_mag(poles[i]);
            if (m > max_mag) max_mag = m;
        }
        return 1.0 - max_mag;
    } else {
        double max_real = -INFINITY;
        for (int i = 0; i < n_poles; i++) {
            if (poles[i].real > max_real) max_real = poles[i].real;
        }
        return -max_real;
    }
}

int stability_lyapunov_continuous(const double *A, int n) {
    /*
     * Check Lyapunov stability for dx/dt = A*x.
     *
     * For a 2x2 system, the trace-determinant conditions are:
     *   det(A) > 0 AND trace(A) < 0  => stable (both eigenvalues in LHP)
     *
     * For higher-order systems, we check if A is Hurwitz via Routh-Hurwitz
     * on the characteristic polynomial det(sI - A).
     *
     * Reference: Khalil, "Nonlinear Systems" (2002), Ch.3
     */
    if (!A || n < 1) return 0;

    if (n == 1) {
        /* A is scalar; stable if A < 0 */
        return A[0] < 0;
    }

    if (n == 2) {
        /* A = [a b; c d] (row-major)
         * det = a*d - b*c
         * trace = a + d */
        double det = A[0]*A[3] - A[1]*A[2];
        double trace = A[0] + A[3];
        return (det > 0) && (trace < 0);
    }

    /* For n >= 3, compute characteristic polynomial and use Routh-Hurwitz */
    /* Characteristic polynomial: det(sI - A) = s^n + c_{n-1}*s^{n-1} + ... + c_0 */
    /* Use the Faddeev-Leverrier algorithm to compute coefficients */

    /* Allocate working matrices */
    double *B = (double*)malloc(n * n * sizeof(double));
    double *C = (double*)malloc(n * n * sizeof(double));

    if (!B || !C) { free(B); free(C); return -1; }

    /* Initialize B = I */
    for (int i = 0; i < n*n; i++) B[i] = (i % (n+1) == 0) ? 1.0 : 0.0;

    double *coeffs = (double*)malloc((n + 1) * sizeof(double));
    coeffs[n] = 1.0;  /* leading coefficient */

    for (int k = 1; k <= n; k++) {
        /* C = A * B */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int l = 0; l < n; l++) {
                    sum += A[i*n + l] * B[l*n + j];
                }
                C[i*n + j] = sum;
            }
        }

        /* coeffs[n-k] = -trace(C)/k */
        double traceC = 0.0;
        for (int i = 0; i < n; i++) traceC += C[i*n + i];
        coeffs[n - k] = -traceC / k;

        /* B = C + coeffs[n-k] * I */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                B[i*n + j] = C[i*n + j];
                if (i == j) B[i*n + j] += coeffs[n - k];
            }
        }
    }

    free(B); free(C);

    /* Apply Routh-Hurwitz to characteristic polynomial */
    /* coeffs are [c_0, c_1, ..., c_{n-1}, 1] where coeffs[k] is coeff of s^k */
    int stable = stability_routh_is_stable(coeffs, n);

    free(coeffs);
    return stable;
}

/* ========================================================================
 * L5: Gain Margin and Phase Margin
 * ======================================================================== */

double stability_gain_margin_db(const rational_func_t *L) {
    /*
     * Compute Gain Margin (GM) from loop transfer function L(s) = G(s)*H(s).
     *
     * GM = 1/|L(j*w_pc)| where w_pc = phase crossover frequency (arg(L) = -180 deg).
     * GM_dB = 20*log10(GM) = -20*log10(|L(j*w_pc)|)
     *
     * Positive GM implies stability.
     * Typical design: GM > 6 dB.
     *
     * Algorithm: Search for phase crossover frequency, then compute magnitude.
     * Reference: Ogata (2010), Section 8-6
     */
    if (!L || L->den.degree < 0) return -INFINITY;

    /* Search for phase crossover (phase = -180 deg = -PI rad) */
    double w_low = 0.01, w_high = 1000.0;
    int found = 0;
    double w_pc = 0.0;

    /* Coarse grid search */
    for (int i = 0; i < 500; i++) {
        double w = w_low * pow(w_high / w_low, (double)i / 499.0);
        double phase = rational_phase_at(L, w, 0);
        /* Phase is in [-pi, pi]. We want phase = -pi (or pi) */
        double phase_deg = phase * 180.0 / M_PI;

        /* Look for phase crossing -180 deg */
        if (fabs(phase_deg + 180.0) < 2.0 || fabs(phase_deg - 180.0) < 2.0) {
            w_pc = w;
            found = 1;
            break;
        }
    }

    if (!found) {
        /* Phase never reaches -180 deg -> infinite gain margin (stable) */
        return INFINITY;
    }

    double mag_pc = rational_mag_at(L, w_pc, 0);
    if (mag_pc < 1e-30) return INFINITY;

    return -20.0 * log10(mag_pc);
}

double stability_phase_margin_deg(const rational_func_t *L) {
    /*
     * Compute Phase Margin (PM) from loop transfer function L(s).
     *
     * PM = 180 + arg(L(j*w_gc)) [degrees]
     * where w_gc = gain crossover frequency (|L(j*w_gc)| = 1 = 0 dB).
     *
     * Positive PM implies stability.
     * Typical design: PM > 45 degrees for good transient response.
     *
     * Reference: Ogata (2010), Section 8-6
     */
    if (!L || L->den.degree < 0) return -180.0;

    /* Search for gain crossover (magnitude = 1 = 0 dB) */
    double w_low = 0.01, w_high = 1000.0;
    int found = 0;
    double w_gc = 0.0;

    /* Coarse grid search for magnitude close to 1 */
    for (int i = 0; i < 500; i++) {
        double w = w_low * pow(w_high / w_low, (double)i / 499.0);
        double mag = rational_mag_at(L, w, 0);
        double mag_db = 20.0 * log10(mag);

        if (fabs(mag_db) < 1.0) {
            w_gc = w;
            found = 1;
            break;
        }
    }

    if (!found) {
        /* Try bisection refinement around the crossing */
        double a = w_low, b = w_high;
        double mag_a = rational_mag_at(L, a, 0);
        double mag_b = rational_mag_at(L, b, 0);

        if ((mag_a - 1.0) * (mag_b - 1.0) < 0) {
            /* Sign change: refine with bisection */
            for (int iter = 0; iter < 50; iter++) {
                double c = (a + b) / 2.0;
                double mag_c = rational_mag_at(L, c, 0);
                if (fabs(mag_c - 1.0) < 0.01) { w_gc = c; found = 1; break; }
                if ((mag_a - 1.0) * (mag_c - 1.0) < 0) { b = c; mag_b = mag_c; }
                else { a = c; mag_a = mag_c; }
            }
        }
    }

    if (!found) return -180.0;

    double phase_gc = rational_phase_at(L, w_gc, 0);
    double pm = M_PI + phase_gc;  /* PM in radians */
    return pm * 180.0 / M_PI;
}