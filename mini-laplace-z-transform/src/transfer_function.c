/*
 * transfer_function.c - Transfer Function Analysis
 *
 * Implements time-domain response computation, frequency response analysis,
 * system interconnection laws, and transient response metrics for LTI systems
 * described by transfer functions H(s) and H(z).
 *
 * Knowledge Coverage:
 *   L1: Transfer function definition, impulse/step response
 *   L2: Frequency response, DC gain, cutoff frequency
 *   L4: Series, parallel, feedback interconnection laws
 *   L6: Canonical analysis problems (settling time, overshoot, damping)
 *
 * Reference: Ogata, "Modern Control Engineering" (2010), Ch.3, Ch.5
 *            Franklin, Powell, Emami-Naeini, "Feedback Control" (2019)
 *            Oppenheim & Schafer, "DSP" (2010), Ch.5
 */

#include "laplace_z_transform_core.h"
#include "transfer_function.h"
#include "partial_fraction.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* ========================================================================
 * L1: Transfer Function Construction
 * ======================================================================== */

rational_func_t tf_create(const double *num_coeffs, int num_len,
                           const double *den_coeffs, int den_len) {
    polynomial_t num = polynomial_from_coeffs(num_coeffs, num_len);
    polynomial_t den = polynomial_from_coeffs(den_coeffs, den_len);
    rational_func_t H; H.num = num; H.den = den;
    return H;
}

/* ========================================================================
 * L6: Time-Domain Response Computation
 * ======================================================================== */

double tf_impulse_response(const rational_func_t *H, double t) {
    /*
     * Impulse response h(t) = L^{-1}{H(s)}
     *
     * For a rational H(s) = sum r_k/(s-p_k),
     * h(t) = sum r_k * exp(p_k * t) * u(t)
     *
     * Algorithm: partial fraction expansion + inverse transform table lookup.
     * Reference: Oppenheim & Willsky (1997), Table 9.2
     */
    if (!H || t < 0.0) return 0.0;
    if (H->den.degree < 0) return 0.0;

    partial_fraction_t pf = pf_create(0);
    if (pf_decompose_laplace(H, &pf) != 0) {
        pf_free(&pf);
        return 0.0;
    }

    double result = 0.0;

    /* Direct term (improper part): represents impulses and derivatives of impulses */
    /* For simplicity, only evaluate at t>0 (impulse contribution is at t=0) */
    if (t > 0.0) {
        /* Direct polynomial contributes nothing for t > 0 in impulse response */
    }

    /* Sum contributions from partial fraction terms */
    for (int i = 0; i < pf.n_terms; i++) {
        pf_term_t *term = &pf.terms[i];
        complex_t p = term->pole;
        complex_t r = term->residue;
        int m = term->multiplicity;

        /* r/(s-p)^m <-> r * t^{m-1}/(m-1)! * exp(p*t) * u(t) */
        double factorial = 1.0;
        for (int k = 2; k < m; k++) factorial *= k;

        /* exp(p*t) = exp((sigma + jw)*t) = exp(sigma*t) * (cos(wt) + j*sin(wt)) */
        double exp_real = exp(p.real * t);
        double cos_wt = cos(p.imag * t);
        double sin_wt = sin(p.imag * t);

        double t_pow = (m == 1) ? 1.0 : pow(t, m - 1) / factorial;

        /* Re[r * exp(p*t)] * t^{m-1}/(m-1)! */
        double term_val = t_pow * exp_real *
            (r.real * cos_wt - r.imag * sin_wt);

        result += term_val;
    }

    pf_free(&pf);
    return result;
}

double tf_step_response(const rational_func_t *H, double t) {
    /*
     * Step response: y_step(t) = L^{-1}{H(s)/s}
     *
     * For H(s) = sum r_k/(s-p_k),
     * H(s)/s = r_0/s + sum r_k'/(s-p_k)
     *
     * y_step(t) = r_0*u(t) + sum r_k' * exp(p_k * t) * u(t)
     *
     * Algorithm: multiply H(s) by 1/s, then partial fraction expand.
     * Reference: Ogata (2010), Section 5-2
     */
    if (!H || t < 0.0) return 0.0;
    if (H->den.degree < 0) return 0.0;

    /* Build H(s)/s: multiply denominator by s */
    double s_coeffs[2] = {0.0, 1.0};
    polynomial_t s_poly = polynomial_from_coeffs(s_coeffs, 2);
    polynomial_t new_den = polynomial_mul(&s_poly, &H->den);
    polynomial_free(&s_poly);

    rational_func_t H_over_s;
    H_over_s.num = polynomial_copy(&H->num);
    H_over_s.den = new_den;

    /* Use impulse response of H(s)/s (which gives step response of H(s)) */
    double result = tf_impulse_response(&H_over_s, t);

    rational_free(&H_over_s);
    return result;
}

double tf_discrete_impulse_response(const rational_func_t *H, int n) {
    /*
     * Discrete impulse response: h[n] = Z^{-1}{H(z)}
     *
     * For H(z^{-1}) = sum r_k/(1 - p_k*z^{-1}),
     * h[n] = sum r_k * p_k^n * u[n]
     *
     * Reference: Oppenheim & Schafer (2010), Section 3.3.2
     */
    if (!H || n < 0) return 0.0;
    if (H->den.degree < 0) return 0.0;

    partial_fraction_t pf = pf_create(0);
    if (pf_decompose_z(H, &pf) != 0) {
        pf_free(&pf);
        return 0.0;
    }

    double result = 0.0;

    /* Direct term for n=0 (delta[n]) */
    if (n == 0 && pf.direct.degree >= 0) {
        result += pf.direct.coeff[0];
    }

    /* Sum contributions from partial fraction terms */
    for (int i = 0; i < pf.n_terms; i++) {
        pf_term_t *term = &pf.terms[i];
        complex_t p = term->pole;
        complex_t r = term->residue;
        int m = term->multiplicity;

        if (n == 0 && m == 1) {
            /* r/(1 - p*z^{-1}) -> r * delta[n] + r*p^n*u[n-1]
             * For n=0: h[0] = r */
            result += r.real;
            continue;
        }

        /* r/(1 - p*z^{-1})^m <-> r * C(n+m-1, m-1) * p^n * u[n]
         * where C(n+m-1, m-1) = (n+m-1)!/(n!*(m-1)!) */
        double binom = 1.0;
        if (m > 1) {
            for (int k = 1; k < m; k++) {
                binom *= (n + k) / (double)k;
            }
        }

        /* p^n: for complex p = |p|*exp(j*theta), p^n = |p|^n * exp(j*n*theta) */
        double p_mag = complex_mag(p);
        double p_phase = complex_phase(p);
        double pn_real = pow(p_mag, n) * cos(n * p_phase);
        double pn_imag = pow(p_mag, n) * sin(n * p_phase);

        double term_val = binom * (r.real * pn_real - r.imag * pn_imag);
        result += term_val;
    }

    pf_free(&pf);
    return result;
}

double tf_discrete_step_response(const rational_func_t *H, int n) {
    /*
     * Discrete step response: y_step[n] = Z^{-1}{H(z) * z/(z-1)}
     *
     * Multiply H(z) by z/(z-1) = 1/(1 - z^{-1}), then inverse Z-transform.
     */
    if (!H || n < 0) return 0.0;
    if (H->den.degree < 0) return 0.0;

    /* H(z) * 1/(1 - z^{-1}): multiply denominator by (1 - z^{-1}) */
    double step_den_coeffs[2] = {1.0, -1.0};
    polynomial_t step_den = polynomial_from_coeffs(step_den_coeffs, 2);
    polynomial_t new_den = polynomial_mul(&step_den, &H->den);
    polynomial_free(&step_den);

    rational_func_t H_step;
    H_step.num = polynomial_copy(&H->num);
    H_step.den = new_den;

    double result = tf_discrete_impulse_response(&H_step, n);

    rational_free(&H_step);
    return result;
}

/* ========================================================================
 * L2: Frequency Response Analysis
 * ======================================================================== */

void tf_freq_response_mag(const rational_func_t *H, const double *freqs,
                           double *mags, int n_freqs) {
    for (int i = 0; i < n_freqs; i++) {
        mags[i] = rational_mag_at(H, freqs[i], 0);
    }
}

void tf_freq_response_phase(const rational_func_t *H, const double *freqs,
                             double *phases, int n_freqs) {
    for (int i = 0; i < n_freqs; i++) {
        phases[i] = rational_phase_at(H, freqs[i], 0);
    }
}

double tf_bode_mag_db(const rational_func_t *H, double omega) {
    double mag = rational_mag_at(H, omega, 0);
    if (mag < 1e-30) return -600.0;  /* effectively -inf dB */
    return 20.0 * log10(mag);
}

double tf_bode_phase_deg(const rational_func_t *H, double omega) {
    return rational_phase_at(H, omega, 0) * 180.0 / M_PI;
}

double tf_cutoff_frequency(const rational_func_t *H, double f_low, double f_high) {
    /*
     * Find -3dB cutoff frequency using bisection search.
     *
     * At cutoff: |H(j*wc)| = |H(0)|/sqrt(2)
     * i.e., |H(j*wc)|^2 = |H(0)|^2 / 2
     *
     * Reference: Sedra & Smith, Microelectronic Circuits (2020), Section 1.6
     */
    double dc_gain = rational_mag_at(H, 0.0, 0);
    double target = dc_gain * dc_gain / 2.0;

    double a = f_low, b = f_high;
    double fa = rational_mag_at(H, a, 0);
    double fb = rational_mag_at(H, b, 0);

    /* Check if cutoff is within range */
    if ((fa*fa - target) * (fb*fb - target) > 0) {
        /* No sign change - cutoff may be outside range */
        return -1.0;
    }

    /* Bisection: O(log((b-a)/epsilon)) */
    for (int iter = 0; iter < 50; iter++) {
        double c = (a + b) / 2.0;
        double fc = rational_mag_at(H, c, 0);
        double val = fc * fc - target;

        if (fabs(b - a) < 1e-6) return c;

        double fa_val = fa * fa - target;
        if (fa_val * val < 0) {
            b = c; fb = fc;
        } else {
            a = c; fa = fc;
        }
    }
    return (a + b) / 2.0;
}

double tf_dc_gain(const rational_func_t *H, int is_z_domain) {
    if (is_z_domain) {
        /* DC gain: |H(e^{j0})| = |H(1)| */
        return rational_mag_at(H, 0.0, 1);
    } else {
        /* DC gain: |H(0)| */
        return rational_mag_at(H, 0.0, 0);
    }
}

/* ========================================================================
 * L4: System Interconnection Laws
 * ======================================================================== */

rational_func_t tf_series(const rational_func_t *H1, const rational_func_t *H2) {
    /* H(s) = H1(s) * H2(s)
     * Reference: Ogata (2010), Section 3-2 */
    return rational_mul(H1, H2);
}

rational_func_t tf_parallel(const rational_func_t *H1, const rational_func_t *H2) {
    /* H(s) = H1(s) + H2(s)
     * Reference: Ogata (2010), Section 3-2 */
    return rational_add(H1, H2);
}

rational_func_t tf_feedback(const rational_func_t *G, const rational_func_t *H) {
    /*
     * Negative feedback: H_cl(s) = G(s) / (1 + G(s)*H(s))
     *
     * This is the fundamental closed-loop transfer function.
     *
     * Derivation:
     *   Y(s) = G(s) * E(s)
     *   E(s) = X(s) - H(s) * Y(s)
     *   => Y(s) = G(s)*(X(s) - H(s)*Y(s))
     *   => Y(s)*(1 + G(s)*H(s)) = G(s)*X(s)
     *   => H_cl(s) = Y(s)/X(s) = G(s)/(1 + G(s)*H(s))
     *
     * Reference: Ogata (2010), Section 3-3
     */
    /* Compute 1 + G*H = (Dg*Dh + Ng*Nh)/(Dg*Dh) */
    polynomial_t DgDh = polynomial_mul(&G->den, &H->den);
    polynomial_t NgNh = polynomial_mul(&G->num, &H->num);
    polynomial_t open_loop_num = polynomial_add(&DgDh, &NgNh);
    polynomial_t open_loop_den = DgDh;

    /* H_cl = G / (1 + G*H) = (Ng/Dg) / (open_loop_num/open_loop_den)
     *      = Ng * open_loop_den / (Dg * open_loop_num) */
    polynomial_t cl_num = polynomial_mul(&G->num, &open_loop_den);
    polynomial_t cl_den = polynomial_mul(&G->den, &open_loop_num);

    polynomial_free(&NgNh);
    polynomial_free(&open_loop_num);
    polynomial_free(&open_loop_den);

    rational_func_t result; result.num = cl_num; result.den = cl_den;
    return result;
}

rational_func_t tf_feedback_positive(const rational_func_t *G, const rational_func_t *H) {
    /* Positive feedback: H_cl(s) = G(s) / (1 - G(s)*H(s)) */
    polynomial_t DgDh = polynomial_mul(&G->den, &H->den);
    polynomial_t NgNh = polynomial_mul(&G->num, &H->num);
    polynomial_t open_loop_num = polynomial_sub(&DgDh, &NgNh);
    polynomial_t open_loop_den = DgDh;

    polynomial_t cl_num = polynomial_mul(&G->num, &open_loop_den);
    polynomial_t cl_den = polynomial_mul(&G->den, &open_loop_num);

    polynomial_free(&NgNh);
    polynomial_free(&open_loop_num);
    polynomial_free(&open_loop_den);

    rational_func_t result; result.num = cl_num; result.den = cl_den;
    return result;
}

polynomial_t tf_char_eq_feedback(const rational_func_t *G, const rational_func_t *H) {
    /*
     * Characteristic equation for feedback: 1 + G(s)*H(s) = 0
     *
     * Returns the polynomial: Dg*Dh + Ng*Nh = 0
     *
     * The roots of this polynomial are the closed-loop poles.
     * Reference: Ogata (2010), Section 5-3
     */
    polynomial_t DgDh = polynomial_mul(&G->den, &H->den);
    polynomial_t NgNh = polynomial_mul(&G->num, &H->num);
    polynomial_t char_eq = polynomial_add(&DgDh, &NgNh);

    polynomial_free(&DgDh);
    polynomial_free(&NgNh);
    return char_eq;
}

/* ========================================================================
 * L2: Transient Response Metrics (Canonical Second-Order Analysis)
 * ======================================================================== */

int tf_dominant_pole(const rational_func_t *H, complex_t *dominant_pole) {
    /*
     * Find the dominant pole: the one closest to the jw-axis.
     * For continuous-time: pole with largest real part (least negative).
     * For discrete-time: pole with largest magnitude (closest to unit circle).
     *
     * Reference: Franklin et al. (2019), Section 3.3
     */
    if (!H || H->den.degree <= 0 || !dominant_pole) return -1;

    pole_zero_t pz = pz_from_rational(H);
    if (pz.n_poles <= 0) { pz_free(&pz); return -1; }

    /* Find pole with maximum real part */
    int idx = 0;
    double max_real = pz.poles[0].real;
    for (int i = 1; i < pz.n_poles; i++) {
        if (pz.poles[i].real > max_real) {
            max_real = pz.poles[i].real;
            idx = i;
        }
    }

    *dominant_pole = pz.poles[idx];
    pz_free(&pz);
    return 0;
}

double tf_settling_time(const rational_func_t *H) {
    /*
     * Estimate 2% settling time from dominant pole.
     *
     * For continuous-time: t_s = 4/|Re(p_dominant)|
     * (time to reach and stay within 2% of final value)
     *
     * For discrete-time: n_s = -4/ln(|p_dominant|)
     *
     * Reference: Ogata (2010), Section 5-4
     */
    complex_t dom_pole;
    if (tf_dominant_pole(H, &dom_pole) != 0) return -1.0;

    double sigma = -dom_pole.real; /* real part should be negative for stable systems */
    if (sigma <= 1e-12) return 1e300; /* marginally stable or unstable */

    return 4.0 / sigma;
}

double tf_overshoot_percent(const rational_func_t *H) {
    /*
     * Estimate percent overshoot from dominant complex pole pair.
     *
     * For a second-order system with damping ratio zeta:
     * Mp = exp(-pi*zeta/sqrt(1-zeta^2)) * 100%
     *
     * where zeta = -Re(p)/|p| (for the dominant pole).
     *
     * Reference: Ogata (2010), Section 5-4
     */
    complex_t dom_pole;
    if (tf_dominant_pole(H, &dom_pole) != 0) return -1.0;

    /* Only compute for complex conjugate poles (oscillatory response) */
    if (fabs(dom_pole.imag) < 1e-12) return 0.0; /* real pole, no overshoot */

    double sigma = -dom_pole.real;
    double wd = fabs(dom_pole.imag);
    double wn = sqrt(sigma*sigma + wd*wd);

    if (wn < 1e-12) return -1.0;

    double zeta = sigma / wn;

    if (zeta <= 0.0 || zeta >= 1.0) return 0.0; /* no overshoot for critically/overdamped */

    double mp = exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
    return mp * 100.0;
}

int tf_damping_params(const rational_func_t *H, double *zeta, double *wn) {
    /*
     * Extract damping ratio zeta and natural frequency wn
     * from the dominant pole pair of a transfer function.
     *
     * For dominant pole p = -sigma + j*wd:
     *   wn = sqrt(sigma^2 + wd^2)
     *   zeta = sigma/wn
     *
     * Reference: Ogata (2010), Section 5-3
     */
    if (!H || !zeta || !wn) return -1;

    complex_t dom_pole;
    if (tf_dominant_pole(H, &dom_pole) != 0) return -1;

    double sigma = -dom_pole.real;
    double wd = fabs(dom_pole.imag);

    *wn = sqrt(sigma*sigma + wd*wd);
    if (*wn < 1e-12) return -1;

    *zeta = sigma / (*wn);
    return 0;
}