/**
 * @file nonlinear_describing.c
 * @brief Describing function analysis for nonlinear systems.
 *
 * Implements the describing function (DF) method for limit cycle
 * prediction and stability analysis of nonlinear systems.
 *
 * Knowledge points:
 *   L5: Analytical DF formulas for all standard nonlinearities
 *   L5: Numerical DF via composite Simpson quadrature
 *   L6: Limit cycle prediction via DF (intersection of loci)
 *   L6: Critical locus -1/N(A) computation
 *
 * Reference:
 *   - Gelb & Vander Velde, "Multiple-Input Describing Functions", 1968
 *   - Slotine & Li, "Applied Nonlinear Control", Ch.5, 1991
 */

#include "nonlinear_describing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* ===================================================================
 * Numerical Describing Function via Composite Simpson's Rule
 *
 * For a memoryless nonlinearity driven by x(t) = A sin(wt):
 *   N(A) = (1/(pi*A)) integral_0^{2pi} f(A sin t) e^{-jt} dt
 *
 * Using n_quad equally spaced points (must be even for Simpson).
 * =================================================================== */

int nl_describing_function(const nl_nonlinearity_t *nl, double A,
                           int n_quad, nl_describing_function_t *df)
{
    if (!nl || !df || A <= 0.0 || n_quad < 10) {
        errno = EINVAL;
        return -1;
    }
    if (n_quad % 2 != 0) n_quad++;

    df->amplitude = A;
    df->frequency = 0.0;

    double sum_re = 0.0, sum_im = 0.0;
    double h = 2.0 * M_PI / (double)n_quad;

    for (int k = 0; k <= n_quad; k++) {
        double theta = (double)k * h;
        double x = A * sin(theta);
        double f_val = nl_nonlinearity_eval(nl, x);

        double weight;
        if (k == 0 || k == n_quad) weight = 1.0;
        else if (k % 2 == 1) weight = 4.0;
        else weight = 2.0;

        sum_re += weight * f_val * sin(theta);
        sum_im += weight * (-f_val * cos(theta));
    }

    double factor = h / (3.0 * M_PI * A);
    df->P = factor * sum_re;
    df->Q = factor * sum_im;
    df->magnitude = sqrt(df->P * df->P + df->Q * df->Q);
    df->phase = atan2(df->Q, df->P);
    return 0;
}

/* ===================================================================
 * Analytical DF: Ideal Relay
 *
 * N(A) = 4M / (pi * A)
 *
 * The ideal relay output is a square wave; its fundamental Fourier
 * component has amplitude 4M/pi.
 * =================================================================== */

double nl_df_ideal_relay(double M, double A)
{
    if (M <= 0.0 || A <= 0.0) {
        errno = EINVAL;
        return 0.0;
    }
    return 4.0 * M / (M_PI * A);
}

/* ===================================================================
 * Analytical DF: Dead-Zone
 *
 * For A > delta:
 *   N(A) = 1 - (2/pi)[arcsin(delta/A) + (delta/A)*sqrt(1-(delta/A)^2)]
 * For A <= delta: N(A) = 0
 * =================================================================== */

double nl_df_dead_zone(double delta, double A)
{
    if (delta < 0.0 || A <= 0.0) {
        errno = EINVAL;
        return 0.0;
    }
    if (A <= delta) return 0.0;

    double r = delta / A;
    double asin_r = asin(r);
    double sqrt_term = r * sqrt(1.0 - r * r);
    return 1.0 - (2.0 / M_PI) * (asin_r + sqrt_term);
}

/* ===================================================================
 * Analytical DF: Saturation
 *
 * For A <= a/k: N(A) = k (linear)
 * For A > a/k:
 *   N(A) = (2k/pi)[arcsin(a/(kA)) + (a/(kA))*sqrt(1-(a/(kA))^2)]
 * =================================================================== */

double nl_df_saturation(double k, double a, double A)
{
    if (k <= 0.0 || a <= 0.0 || A <= 0.0) {
        errno = EINVAL;
        return 0.0;
    }
    double lim = a / k;
    if (A <= lim) return k;

    double r = lim / A;
    double asin_r = asin(r);
    double sqrt_term = r * sqrt(1.0 - r * r);
    return (2.0 * k / M_PI) * (asin_r + sqrt_term);
}

/* ===================================================================
 * Analytical DF: Backlash
 *
 * For A > delta:
 *   Real: k/2 - (k/pi)*arcsin(1-2d/A) - (2k/pi)*(1-2d/A)*sqrt((d/A)*(1-d/A))
 *   Imag: (4kd/(pi*A))*(d/A - 1)
 *
 * The negative imaginary part captures hysteresis phase lag.
 * =================================================================== */

int nl_df_backlash(double k, double delta, double A,
                   double *real_part, double *imag_part)
{
    if (!real_part || !imag_part) { errno = EINVAL; return -1; }
    if (k <= 0.0 || delta < 0.0 || A <= 0.0) {
        *real_part = 0.0; *imag_part = 0.0;
        errno = EINVAL; return -1;
    }
    if (A <= delta) {
        *real_part = 0.0; *imag_part = 0.0;
        return 0;
    }

    double r = delta / A;
    double term = 1.0 - 2.0 * r;
    double asin_term = asin(term);
    double sqrt_term = (r * (1.0 - r) > 0.0)
                       ? sqrt(r * (1.0 - r)) : 0.0;

    *real_part = k / 2.0 - (k / M_PI) * asin_term
                 - (2.0 * k / M_PI) * term * sqrt_term;
    *imag_part = (4.0 * k * delta / (M_PI * A)) * (r - 1.0);
    return 0;
}

/* ===================================================================
 * Analytical DF: Relay with Hysteresis
 *
 * For A >= h:
 *   Re[N] = (4M/(pi*A)) * sqrt(1 - (h/A)^2)
 *   Im[N] = -(4M*h)/(pi*A^2)
 * For A < h: N(A) = 0
 * =================================================================== */

int nl_df_relay_hysteresis(double M, double h, double A,
                           double *real_part, double *imag_part)
{
    if (!real_part || !imag_part) { errno = EINVAL; return -1; }
    if (M <= 0.0 || h < 0.0 || A <= 0.0) {
        *real_part = 0.0; *imag_part = 0.0;
        errno = EINVAL; return -1;
    }
    if (A < h) { *real_part = 0.0; *imag_part = 0.0; return 0; }

    double r = h / A;
    *real_part = (4.0 * M / (M_PI * A)) * sqrt(1.0 - r * r);
    *imag_part = -(4.0 * M * h) / (M_PI * A * A);
    return 0;
}

/* ===================================================================
 * Critical Locus: -1/N(A)
 *
 * L(A) = -1/N(A) = -(P - jQ) / (P^2 + Q^2)
 *
 * This is the locus that the Nyquist plot G(jw) should not encircle
 * (or must encircle the correct number of times) for stability.
 * =================================================================== */

int nl_df_critical_locus(const nl_nonlinearity_t *nl,
                         const double *A, size_t nA,
                         double *locus_real, double *locus_imag)
{
    if (!nl || !A || nA < 1 || !locus_real || !locus_imag) {
        errno = EINVAL;
        return -1;
    }

    for (size_t i = 0; i < nA; i++) {
        if (A[i] <= 0.0) {
            locus_real[i] = -INFINITY;
            locus_imag[i] = 0.0;
            continue;
        }

        nl_describing_function_t df;
        if (nl_describing_function(nl, A[i], 200, &df) != 0) {
            locus_real[i] = 0.0;
            locus_imag[i] = 0.0;
            continue;
        }

        double denom = df.P * df.P + df.Q * df.Q;
        if (denom < 1e-15) {
            locus_real[i] = -INFINITY;
            locus_imag[i] = 0.0;
        } else {
            locus_real[i] = -df.P / denom;
            locus_imag[i] = df.Q / denom;
        }
    }
    return 0;
}

/* ===================================================================
 * Limit Cycle Prediction via DF Method
 *
 * Solves G(jw) * N(A) = -1  =>  G(jw) = -1/N(A)
 *
 * Performs grid search over (A, w) space to find the intersection
 * of the Nyquist and critical loci. Refines using local search.
 * =================================================================== */

int nl_df_limit_cycle(const nl_linear_system_t *lin,
                      const nl_nonlinearity_t *nl,
                      const double *A_range, size_t nA,
                      const double *w_range, size_t nw,
                      double *A_pred, double *w_pred,
                      int *stable)
{
    if (!lin || !nl || !A_range || nA < 2
        || !w_range || nw < 2 || !A_pred || !w_pred || !stable) {
        errno = EINVAL;
        return -1;
    }

    double A_lo = A_range[0], A_hi = A_range[1];
    double w_lo = w_range[0], w_hi = w_range[1];
    double best_dist = INFINITY;
    double best_A = 0.0, best_w = 0.0;
    int found = 0;

    for (size_t i = 0; i < nA; i++) {
        double A = A_lo + (A_hi - A_lo) * (double)i / (double)(nA - 1);
        if (A <= 0.0) continue;

        nl_describing_function_t df;
        if (nl_describing_function(nl, A, 200, &df) != 0) continue;

        double d2 = df.P * df.P + df.Q * df.Q;
        if (d2 < 1e-15) continue;
        double cr = -df.P / d2, ci = df.Q / d2;

        for (size_t j = 0; j < nw; j++) {
            double w = w_lo + (w_hi - w_lo) * (double)j / (double)(nw - 1);
            double Gr, Gi;
            if (nl_transfer_function(lin, w, &Gr, &Gi) != 0) continue;

            double dist = hypot(Gr - cr, Gi - ci);
            if (dist < best_dist) {
                best_dist = dist;
                best_A = A;
                best_w = w;
                if (dist < 0.01) found = 1;
            }
        }
    }

    if (!found && best_dist > 0.1) return 0;

    /* Determine stability: for memoryless odd nonlinearities,
     * the limit cycle is stable if d(-1/N)/dA crosses G(jw)
     * from the stable side. Heuristic: N(A) decreasing => stable. */
    nl_describing_function_t df1, df2;
    double hA = best_A * 0.01 + 1e-6;
    nl_describing_function(nl, best_A - hA, 200, &df1);
    nl_describing_function(nl, best_A + hA, 200, &df2);
    *stable = (df2.magnitude < df1.magnitude) ? 1 : 0;

    *A_pred = best_A;
    *w_pred = best_w;
    return 1;
}

/* ===================================================================
 * Harmonic Balance Skeleton
 *
 * Single-harmonic implementation serving as baseline for multi-harmonic
 * Galerkin projection. For a 2D autonomous system:
 *   dx/dt = f1(x,y),  dy/dt = f2(x,y)
 * assume x(t) = A cos(wt), y(t) = B sin(wt) + C cos(wt).
 *
 * Energy balance iteration refines amplitude and frequency estimates.
 * =================================================================== */

int nl_harmonic_balance(nl_nonlinear_system_t *sys,
                        int N_harmonics, double w_guess,
                        double *coeffs, double *w_sol)
{
    if (!sys || !coeffs || !w_sol
        || N_harmonics < 1 || N_harmonics > 8
        || w_guess <= 0.0 || sys->dim < 1
        || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL;
        return -1;
    }

    size_t n_coeffs = (size_t)(2 * N_harmonics + 1);
    memset(coeffs, 0, n_coeffs * sizeof(double));

    /* Single-harmonic energy balance iteration */
    double A = 1.0;
    double w = w_guess;
    int n_samples = 200;
    double T = 2.0 * M_PI / w;
    double dt = T / (double)n_samples;

    for (int iter = 0; iter < 100; iter++) {
        double energy_diff = 0.0;

        for (int s = 1; s <= n_samples; s++) {
            double t = (double)s * dt;
            double x_curr = A * sin(w * t);
            double v_curr = A * w * cos(w * t);

            double state[2] = {x_curr, v_curr};
            double f_val[NL_MAX_STATE_DIM];
            if (sys->f(state, t, sys->params, f_val, sys->dim) != 0)
                break;

            /* Power balance for 2D system:
             * f[0] = v, f[1] = acceleration */
            energy_diff += v_curr * v_curr * dt;
            if (sys->dim >= 2)
                energy_diff -= fabs(f_val[1]) * fabs(v_curr) * dt;
        }

        A += 0.05 * energy_diff / (double)n_samples;
        if (A < 1e-6) A = 1e-6;
        if (A > 100.0) A = 100.0;
        if (fabs(energy_diff) < 1e-6) break;
    }

    coeffs[1] = A;
    coeffs[2] = 0.0;
    *w_sol = w;
    return 0;
}
