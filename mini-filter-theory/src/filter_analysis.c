/**
 * filter_analysis.c — Filter Analysis and Performance Metrics
 *
 * L2: Frequency response computation
 * L3: Group delay, phase delay, pole-zero analysis
 * L6: SNR improvement, THD, step/impulse response
 *
 * Reference:
 *   Oppenheim & Willsky (1997) Signals and Systems
 *   Lyons, Understanding Digital Signal Processing (2010)
 *   Jackson, "On the Interaction of Roundoff Noise and
 *            Dynamic Range in Digital Filters" (1970)
 */

#include "filter_analysis.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * L2: Frequency Response at Single Point
 * ================================================================== */

/**
 * fir_freq_response_at — L2
 *
 * H(e^{jw}) = sum_{n=0}^{N-1} h[n] * e^{-j*w*n}
 *           = sum h[n]*cos(w*n) - j*sum h[n]*sin(w*n)
 *
 * Complexity: O(N) per frequency evaluation
 */
freq_resp_point_t fir_freq_response_at(const fir_filter_t *fir, double w) {
    freq_resp_point_t fp;
    memset(&fp, 0, sizeof(fp));
    if (!fir || !fir->coeff) return fp;

    fp.freq_rad = w;

    double real_sum = 0.0, imag_sum = 0.0;
    int n;
    for (n = 0; n < fir->length; n++) {
        real_sum += fir->coeff[n] * cos(w * n);
        imag_sum -= fir->coeff[n] * sin(w * n);
    }

    fp.real_part = real_sum;
    fp.imag_part = imag_sum;
    fp.magnitude = sqrt(real_sum * real_sum + imag_sum * imag_sum);
    fp.magnitude_db = 20.0 * log10(fp.magnitude > 1e-15 ? fp.magnitude : 1e-15);
    fp.phase_rad = atan2(imag_sum, real_sum);
    fp.phase_deg = fp.phase_rad * 180.0 / M_PI;
    fp.group_delay_s = 0.0;
    return fp;
}

/**
 * iir_freq_response_at — L2
 *
 * Evaluates cascade of biquads:
 *   H_k(e^{jw}) = (b0 + b1*e^{-jw} + b2*e^{-j2w})
 *                / (1 + a1*e^{-jw} + a2*e^{-j2w})
 *   H(e^{jw}) = gain * prod H_k(e^{jw})
 *
 * Complexity: O(num_sections) per frequency
 */
freq_resp_point_t iir_freq_response_at(const iir_filter_t *iir, double w) {
    freq_resp_point_t fp;
    memset(&fp, 0, sizeof(fp));
    if (!iir || !iir->sections) return fp;

    fp.freq_rad = w;

    double mag = iir->overall_gain;
    double phase = 0.0;

    int k;
    for (k = 0; k < iir->num_sections; k++) {
        biquad_section_t *bq = &iir->sections[k];
        double cos_w  = cos(w);
        double cos_2w = cos(2.0 * w);
        double sin_w  = sin(w);
        double sin_2w = sin(2.0 * w);

        /* Numerator: b0 + b1*e^{-jw} + b2*e^{-j2w} */
        double num_real = bq->b0 + bq->b1 * cos_w + bq->b2 * cos_2w;
        double num_imag = -(bq->b1 * sin_w + bq->b2 * sin_2w);

        /* Denominator: 1 + a1*e^{-jw} + a2*e^{-j2w} */
        double den_real = 1.0 + bq->a1 * cos_w + bq->a2 * cos_2w;
        double den_imag = -(bq->a1 * sin_w + bq->a2 * sin_2w);

        double den_mag = den_real * den_real + den_imag * den_imag;
        if (den_mag < 1e-30) {
            mag = INFINITY;
            break;
        }

        double num_mag = sqrt(num_real * num_real + num_imag * num_imag);
        mag *= num_mag / sqrt(den_mag);
        phase += atan2(num_imag, num_real) - atan2(den_imag, den_real);
    }

    /* Unwrap phase to be continuous */
    while (phase > M_PI)  phase -= 2.0 * M_PI;
    while (phase < -M_PI) phase += 2.0 * M_PI;

    fp.magnitude = mag;
    fp.magnitude_db = 20.0 * log10(mag > 1e-15 ? mag : 1e-15);
    fp.phase_rad = phase;
    fp.phase_deg = phase * 180.0 / M_PI;

    /* Re{H} and Im{H} */
    fp.real_part = mag * cos(phase);
    fp.imag_part = mag * sin(phase);

    return fp;
}

/* ==================================================================
 * L2: Complete Frequency Response Vector
 * ================================================================== */

freq_resp_t *fir_freq_response(const fir_filter_t *fir, double f_start,
                                double f_stop, int num_points, double fs) {
    if (!fir || num_points < 2 || fs <= 0.0) return NULL;
    if (f_start < 0.0 || f_stop <= f_start) return NULL;

    freq_resp_t *resp = (freq_resp_t *)malloc(sizeof(freq_resp_t));
    if (!resp) return NULL;

    resp->points = (freq_resp_point_t *)malloc(num_points * sizeof(freq_resp_point_t));
    if (!resp->points) { free(resp); return NULL; }

    resp->num_points = num_points;
    resp->f_start = f_start;
    resp->f_stop = f_stop;
    resp->f_resolution = (f_stop - f_start) / (num_points - 1);

    int i;
    for (i = 0; i < num_points; i++) {
        double f = f_start + i * resp->f_resolution;
        double w = 2.0 * M_PI * f / fs;

        resp->points[i] = fir_freq_response_at(fir, w);
        resp->points[i].freq_hz = f;

        /* Compute group delay via numerical differentiation */
        if (i > 0 && i < num_points - 1) {
            double dw = 2.0 * M_PI * resp->f_resolution / fs;
            double phi_plus  = resp->points[i + 1].phase_rad;
            double phi_minus = resp->points[i - 1].phase_rad;
            double dphi = phi_plus - phi_minus;
            /* Phase unwrapping correction */
            if (dphi > M_PI)  dphi -= 2.0 * M_PI;
            if (dphi < -M_PI) dphi += 2.0 * M_PI;
            resp->points[i].group_delay_s = -dphi / (2.0 * dw) / fs;
        }
    }

    return resp;
}

freq_resp_t *iir_freq_response(const iir_filter_t *iir, double f_start,
                                double f_stop, int num_points, double fs) {
    if (!iir || num_points < 2 || fs <= 0.0) return NULL;
    if (f_start < 0.0 || f_stop <= f_start) return NULL;

    freq_resp_t *resp = (freq_resp_t *)malloc(sizeof(freq_resp_t));
    if (!resp) return NULL;

    resp->points = (freq_resp_point_t *)malloc(num_points * sizeof(freq_resp_point_t));
    if (!resp->points) { free(resp); return NULL; }

    resp->num_points = num_points;
    resp->f_start = f_start;
    resp->f_stop = f_stop;
    resp->f_resolution = (f_stop - f_start) / (num_points - 1);

    int i;
    for (i = 0; i < num_points; i++) {
        double f = f_start + i * resp->f_resolution;
        double w = 2.0 * M_PI * f / fs;
        resp->points[i] = iir_freq_response_at(iir, w);
        resp->points[i].freq_hz = f;
    }

    return resp;
}

void freq_resp_free(freq_resp_t *resp) {
    if (resp) {
        free(resp->points);
        free(resp);
    }
}

freq_resp_point_t tf_continuous_response_at(const tf_continuous_t *tf,
                                              double omega_rad_s) {
    freq_resp_point_t fp;
    memset(&fp, 0, sizeof(fp));
    if (!tf) return fp;

    double complex H = tf_continuous_eval(tf, omega_rad_s * I);
    fp.freq_rad = omega_rad_s;
    fp.magnitude = cabs(H);
    fp.magnitude_db = 20.0 * log10(fp.magnitude > 1e-15 ? fp.magnitude : 1e-15);
    fp.phase_rad = carg(H);
    fp.phase_deg = fp.phase_rad * 180.0 / M_PI;
    fp.real_part = creal(H);
    fp.imag_part = cimag(H);
    return fp;
}

freq_resp_t *tf_continuous_freq_response(const tf_continuous_t *tf,
                                           double f_start, double f_stop,
                                           int num_points) {
    if (!tf || num_points < 2) return NULL;

    freq_resp_t *resp = (freq_resp_t *)malloc(sizeof(freq_resp_t));
    if (!resp) return NULL;

    resp->points = (freq_resp_point_t *)malloc(num_points * sizeof(freq_resp_point_t));
    if (!resp->points) { free(resp); return NULL; }

    resp->num_points = num_points;
    resp->f_start = f_start;
    resp->f_stop = f_stop;
    resp->f_resolution = (f_stop - f_start) / (num_points - 1);

    int i;
    for (i = 0; i < num_points; i++) {
        double f = f_start + i * resp->f_resolution;
        double w = 2.0 * M_PI * f;
        resp->points[i] = tf_continuous_response_at(tf, w);
        resp->points[i].freq_hz = f;
    }

    return resp;
}

/* ==================================================================
 * L2: Filter Performance Metrics
 * ================================================================== */

double passband_flatness(const freq_resp_t *resp,
                          double f_pass_start, double f_pass_stop) {
    if (!resp || !resp->points || resp->num_points < 2) return -1.0;

    double max_db = -INFINITY, min_db = INFINITY;
    int i;
    for (i = 0; i < resp->num_points; i++) {
        if (resp->points[i].freq_hz >= f_pass_start &&
            resp->points[i].freq_hz <= f_pass_stop) {
            double db = resp->points[i].magnitude_db;
            if (db > max_db) max_db = db;
            if (db < min_db) min_db = db;
        }
    }

    if (max_db == -INFINITY) return -1.0;
    return max_db - min_db;
}

double stopband_attenuation(const freq_resp_t *resp,
                             double f_stop_start, double f_stop_stop) {
    if (!resp || !resp->points || resp->num_points < 2) return -1.0;

    double passband_max = -INFINITY;
    int i;
    /* Find passband reference: max magnitude below f_stop_start */
    for (i = 0; i < resp->num_points; i++) {
        if (resp->points[i].freq_hz < f_stop_start) {
            if (resp->points[i].magnitude_db > passband_max) {
                passband_max = resp->points[i].magnitude_db;
            }
        }
    }

    /* Find minimum in stopband */
    double stopband_min = INFINITY;
    for (i = 0; i < resp->num_points; i++) {
        if (resp->points[i].freq_hz >= f_stop_start &&
            resp->points[i].freq_hz <= f_stop_stop) {
            if (resp->points[i].magnitude_db < stopband_min) {
                stopband_min = resp->points[i].magnitude_db;
            }
        }
    }

    if (passband_max == -INFINITY || stopband_min == INFINITY) return -1.0;
    return passband_max - stopband_min;
}

double find_3db_cutoff(const freq_resp_t *resp, double ref_db) {
    if (!resp || resp->num_points < 2) return -1.0;

    double target = ref_db - 3.0103; /* -3.01 dB */
    int i;
    for (i = 1; i < resp->num_points; i++) {
        double db0 = resp->points[i - 1].magnitude_db;
        double db1 = resp->points[i].magnitude_db;
        if ((db0 >= target && db1 <= target) ||
            (db0 <= target && db1 >= target)) {
            /* Linear interpolation */
            double f0 = resp->points[i - 1].freq_hz;
            double f1 = resp->points[i].freq_hz;
            return f0 + (target - db0) * (f1 - f0) / (db1 - db0);
        }
    }
    return -1.0;
}

/**
 * transition_bandwidth — L2 find transition band from frequency response
 *
 * The transition band is the frequency region between the passband edge
 * and stopband edge. This function locates these edges automatically
 * from the frequency response data.
 *
 * Algorithm:
 * 1. Find passband edge: where magnitude first drops by passband_ripple_db
 *    from the passband reference level.
 * 2. Find stopband edge: where magnitude first reaches stopband_atten_db
 *    attenuation below the passband reference.
 * 3. Transition bandwidth = stopband_edge - passband_edge
 *
 * This metric is critical for filter design validaton and comparison.
 *
 * @param resp              Frequency response data
 * @param passband_ripple_db Expected passband ripple (dB)
 * @param stopband_atten_db  Expected minimum stopband attenuation (dB)
 * @return                   Transition bandwidth in Hz, or -1 on error
 * Complexity: O(num_points)
 */
double transition_bandwidth(const freq_resp_t *resp,
                             double passband_ripple_db,
                             double stopband_atten_db) {
    if (!resp || !resp->points || resp->num_points < 4)
        return -1.0;
    if (passband_ripple_db < 0.0 || stopband_atten_db < 0.0)
        return -1.0;

    /* Find passband reference: max magnitude in the lower frequency region */
    double pass_ref = -INFINITY;
    int i, pass_ref_idx = 0;
    for (i = 0; i < resp->num_points / 3; i++) {
        if (resp->points[i].magnitude_db > pass_ref) {
            pass_ref = resp->points[i].magnitude_db;
            pass_ref_idx = i;
        }
    }
    if (pass_ref == -INFINITY) return -1.0;

    /* Find passband edge: first point ≤ pass_ref - passband_ripple_db
     * after the reference point */
    double pass_edge_thresh = pass_ref - passband_ripple_db;
    double pass_edge_hz = -1.0;
    for (i = pass_ref_idx + 1; i < resp->num_points; i++) {
        if (resp->points[i].magnitude_db <= pass_edge_thresh) {
            /* Linear interpolation for more accurate edge */
            if (i > 0 && resp->points[i-1].magnitude_db > pass_edge_thresh) {
                double f0 = resp->points[i-1].freq_hz;
                double f1 = resp->points[i].freq_hz;
                double m0 = resp->points[i-1].magnitude_db;
                double m1 = resp->points[i].magnitude_db;
                pass_edge_hz = f0 + (pass_edge_thresh - m0) * (f1 - f0) / (m1 - m0);
            } else {
                pass_edge_hz = resp->points[i].freq_hz;
            }
            break;
        }
    }
    if (pass_edge_hz < 0.0) return -1.0;

    /* Find stopband edge: first point ≤ pass_ref - stopband_atten_db
     * after the passband edge */
    double stop_edge_thresh = pass_ref - stopband_atten_db;
    double stop_edge_hz = -1.0;
    for (i = (int)(pass_edge_hz / resp->f_resolution) + 1;
         i < resp->num_points; i++) {
        if (resp->points[i].magnitude_db <= stop_edge_thresh) {
            if (i > 0 && resp->points[i-1].magnitude_db > stop_edge_thresh) {
                double f0 = resp->points[i-1].freq_hz;
                double f1 = resp->points[i].freq_hz;
                double m0 = resp->points[i-1].magnitude_db;
                double m1 = resp->points[i].magnitude_db;
                stop_edge_hz = f0 + (stop_edge_thresh - m0) * (f1 - f0) / (m1 - m0);
            } else {
                stop_edge_hz = resp->points[i].freq_hz;
            }
            break;
        }
    }
    if (stop_edge_hz < 0.0) return -1.0;

    return stop_edge_hz - pass_edge_hz;
}

double effective_noise_bandwidth(const freq_resp_t *resp) {
    if (!resp || resp->num_points < 2) return -1.0;

    /* ENB = integral(|H(f)|^2 df) / |H(f0)|^2 */
    double integral = 0.0;
    double max_mag_sq = 0.0;
    int i;
    double df = resp->f_resolution;

    for (i = 0; i < resp->num_points; i++) {
        double mag_sq = resp->points[i].magnitude * resp->points[i].magnitude;
        integral += mag_sq * df;
        if (mag_sq > max_mag_sq) max_mag_sq = mag_sq;
    }

    if (max_mag_sq < 1e-15) return -1.0;
    return integral / max_mag_sq;
}

/* ==================================================================
 * L3: Group Delay and Phase Analysis
 * ================================================================== */

double group_delay_at_fir(const fir_filter_t *fir, double w, double dw) {
    if (!fir || dw <= 0.0) return -1.0;

    freq_resp_point_t fp_plus  = fir_freq_response_at(fir, w + dw);
    freq_resp_point_t fp_minus = fir_freq_response_at(fir, w - dw);

    double dphi = fp_plus.phase_rad - fp_minus.phase_rad;
    if (dphi > M_PI)  dphi -= 2.0 * M_PI;
    if (dphi < -M_PI) dphi += 2.0 * M_PI;

    return -dphi / (2.0 * dw);
}

double group_delay_at_iir(const iir_filter_t *iir, double w, double dw) {
    if (!iir || dw <= 0.0) return -1.0;

    freq_resp_point_t fp_plus  = iir_freq_response_at(iir, w + dw);
    freq_resp_point_t fp_minus = iir_freq_response_at(iir, w - dw);

    double dphi = fp_plus.phase_rad - fp_minus.phase_rad;
    if (dphi > M_PI)  dphi -= 2.0 * M_PI;
    if (dphi < -M_PI) dphi += 2.0 * M_PI;

    return -dphi / (2.0 * dw);
}

double phase_delay_at_fir(const fir_filter_t *fir, double w) {
    if (!fir || fabs(w) < 1e-15) {
        return fir ? fir->group_delay : -1.0;
    }
    freq_resp_point_t fp = fir_freq_response_at(fir, w);
    return -fp.phase_rad / w;
}

double phase_delay_at_iir(const iir_filter_t *iir, double w) {
    if (!iir || fabs(w) < 1e-15) return 0.0;
    freq_resp_point_t fp = iir_freq_response_at(iir, w);
    return -fp.phase_rad / w;
}

double group_delay_ripple(const freq_resp_t *resp,
                           double f_start, double f_stop) {
    if (!resp || resp->num_points < 2) return -1.0;

    double max_gd = -INFINITY, min_gd = INFINITY;
    int i;
    for (i = 0; i < resp->num_points; i++) {
        if (resp->points[i].freq_hz >= f_start &&
            resp->points[i].freq_hz <= f_stop) {
            double gd = resp->points[i].group_delay_s;
            if (gd > max_gd) max_gd = gd;
            if (gd < min_gd) min_gd = gd;
        }
    }
    if (max_gd == -INFINITY) return -1.0;
    return max_gd - min_gd;
}

/* ==================================================================
 * L3: Pole-Zero Analysis
 * ================================================================== */

double pz_max_q(const pz_map_t *pz, double *qvals) {
    if (!pz) return -1.0;
    double max_q = 0.0;
    int i;
    for (i = 0; i < pz->num_pairs; i++) {
        double q = pz->pz_pairs[i].q_pole;
        if (qvals) qvals[i] = q;
        if (q > max_q) max_q = q;
    }
    return max_q;
}

int pz_is_minimum_phase(const pz_map_t *pz) {
    if (!pz) return -1;
    /* For continuous-time: all poles AND zeros in LHP (Re < 0) */
    int i;
    for (i = 0; i < pz->num_pairs; i++) {
        if (creal(pz->pz_pairs[i].zero) >= 0.0) return 0;
    }
    for (i = 0; i < pz->num_real_zeros; i++) {
        if (pz->real_zeros[i] >= 0.0) return 0;
    }
    return 1;
}

double pz_stability_margin(const pz_map_t *pz, int is_discrete) {
    if (!pz) return -INFINITY;
    double min_margin = INFINITY;
    int i;
    for (i = 0; i < pz->num_pairs; i++) {
        double margin;
        if (is_discrete) {
            margin = 1.0 - cabs(pz->pz_pairs[i].pole);
        } else {
            margin = -creal(pz->pz_pairs[i].pole);
        }
        if (margin < min_margin) min_margin = margin;
    }
    for (i = 0; i < pz->num_real_poles; i++) {
        double margin = is_discrete ?
            (1.0 - fabs(pz->real_poles[i])) : (-pz->real_poles[i]);
        if (margin < min_margin) min_margin = margin;
    }
    return min_margin;
}

/* ==================================================================
 * L4: Complete Stability Analysis
 * ================================================================== */

int poles_are_stable(const double complex *poles, int num_poles,
                      int is_discrete) {
    if (!poles || num_poles < 1) return -1;
    int i;
    for (i = 0; i < num_poles; i++) {
        if (is_discrete) {
            if (cabs(poles[i]) >= 1.0) return 0;
        } else {
            if (creal(poles[i]) >= 0.0) return 0;
        }
    }
    return 1;
}

/**
 * spectral_factor — L4 cepstral spectral factorization
 *
 * Given |H(e^{jw})|^2, finds a minimum-phase H(z) such that
 * |H(e^{jw})|^2 = H(e^{jw}) * H*(e^{jw}).
 *
 * Cepstral method:
 * 1. Take inverse Fourier transform of log(|H(w)|^2) → cepstrum c[n]
 * 2. Causal cepstrum: c_causal[0] = c[0]/2, c_causal[n] = c[n] for n>0
 * 3. Exponentiate: H(z) = exp(FFT(c_causal[n]))
 *
 * Simplified implementation using the autocorrelation method:
 * Given the magnitude-squared coefficients as a symmetric sequence,
 * find the minimum-phase spectral factor via Levinson-Durbin recursion.
 *
 * This is fundamental to:
 * - Optimal FIR filter design (FIRs that match a desired spectrum)
 * - Linear predictive coding (LPC) in speech processing
 * - Parametric spectral estimation
 *
 * Reference: Oppenheim & Schafer (2010), Ch. 12
 *             Sayed, Fundamentals of Adaptive Filtering (2003)
 * Complexity: O(half_len^2)
 */
int spectral_factor(const double *mag_sq_coeff, int half_len, double *H) {
    if (!mag_sq_coeff || !H || half_len < 1 || half_len > 500)
        return FILTER_ERR_NULL_PTR;

    /* Use the autocorrelation method:
     * Given r[0..p] (autocorrelation), find minimum-phase filter a[1..p]
     * via Levinson-Durbin, then H[0..p] are the filter coefficients.
     *
     * The magnitude-squared coefficients form a symmetric autocorrelation:
     * r[k] = mag_sq_coeff[k] for k = 0..half_len-1
     * r[-k] = r[k]
     */

    /* Levinson-Durbin recursion for linear prediction coefficients */
    int p = half_len - 1;
    if (p < 0) return FILTER_ERR_INVALID_ORDER;

    double *a = (double *)calloc((p + 1) * (p + 1), sizeof(double));
    double *E = (double *)calloc(p + 1, sizeof(double));
    if (!a || !E) { free(a); free(E); return FILTER_ERR_MEMORY; }

    /* Initialize: a[0][0] = 1, E[0] = r[0] */
    a[0 * (p + 1) + 0] = 1.0;
    E[0] = mag_sq_coeff[0];
    if (E[0] < 1e-15) { free(a); free(E); return FILTER_ERR_NUMERICAL; }

    int i, j;
    for (i = 0; i < p; i++) {
        /* Compute reflection coefficient k_{i+1} */
        double ki_num = mag_sq_coeff[i + 1];
        for (j = 0; j < i; j++) {
            ki_num -= a[i * (p + 1) + j] * mag_sq_coeff[i - j];
        }

        if (fabs(E[i]) < 1e-15) { free(a); free(E); return FILTER_ERR_NUMERICAL; }

        double ki = -ki_num / E[i];

        /* Update a_{i+1}[j] */
        int idx_new = (i + 1) * (p + 1);
        a[idx_new + 0] = 1.0;
        for (j = 1; j <= i; j++) {
            a[idx_new + j] = a[i * (p + 1) + j]
                           + ki * a[i * (p + 1) + (i - j + 1)];
        }
        a[idx_new + i + 1] = ki;

        /* Update prediction error */
        E[i + 1] = E[i] * (1.0 - ki * ki);
    }

    /* The minimum-phase spectral factor H is the prediction error filter */
    /* H[0] = 1, H[1..p] = a_p[1..p] (with sign flip for LPC convention) */
    H[0] = 1.0;
    for (j = 1; j <= p; j++) {
        H[j] = a[p * (p + 1) + j];
    }

    free(a); free(E);
    return FILTER_OK;
}

/* ==================================================================
 * L6: Filter Performance Comparison
 * ================================================================== */

double freq_resp_compare(const freq_resp_t *resp1, const freq_resp_t *resp2) {
    if (!resp1 || !resp2) return -1.0;
    int n = (resp1->num_points < resp2->num_points) ?
            resp1->num_points : resp2->num_points;
    if (n < 2) return -1.0;

    double sum_sq = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        double diff = resp1->points[i].magnitude_db
                    - resp2->points[i].magnitude_db;
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq / n);
}

double snr_improvement(const freq_resp_t *resp,
                        double signal_band_low, double signal_band_high) {
    if (!resp || resp->num_points < 2) return 0.0;
    if (signal_band_low >= signal_band_high) return 0.0;

    /* White noise: noise power = integral |H(f)|^2 over all frequencies
     * Signal power = integral |H(f)|^2 over signal band
     * SNR improvement = signal_power/total_power / (band_ratio)
     */
    double signal_power = 0.0, total_power = 0.0;
    int sig_count = 0, total_count = 0;
    int i;

    for (i = 0; i < resp->num_points; i++) {
        double mag_sq = resp->points[i].magnitude * resp->points[i].magnitude;
        total_power += mag_sq;
        total_count++;
        if (resp->points[i].freq_hz >= signal_band_low &&
            resp->points[i].freq_hz <= signal_band_high) {
            signal_power += mag_sq;
            sig_count++;
        }
    }

    if (total_power < 1e-15 || sig_count == 0) return 0.0;
    double band_ratio = (double)sig_count / total_count;
    double snr = (signal_power / total_power) / band_ratio;
    return 10.0 * log10(snr);
}

double filter_thd(const freq_resp_t *resp, double fundamental_hz,
                   int num_harmonics) {
    if (!resp || num_harmonics < 2) return -1.0;

    double fund_mag = 1.0;
    int i;
    /* Find fundamental magnitude */
    for (i = 0; i < resp->num_points; i++) {
        if (fabs(resp->points[i].freq_hz - fundamental_hz) < resp->f_resolution) {
            fund_mag = resp->points[i].magnitude;
            break;
        }
    }
    if (fund_mag < 1e-15) return -1.0;

    /* Sum harmonic powers */
    double harmonic_sum_sq = 0.0;
    int h;
    for (h = 2; h <= num_harmonics; h++) {
        double h_freq = fundamental_hz * h;
        double h_mag = 0.0;
        for (i = 0; i < resp->num_points; i++) {
            if (fabs(resp->points[i].freq_hz - h_freq) < resp->f_resolution) {
                h_mag = resp->points[i].magnitude;
                break;
            }
        }
        harmonic_sum_sq += h_mag * h_mag;
    }

    double thd = sqrt(harmonic_sum_sq) / fund_mag;
    return 20.0 * log10(thd);
}

/* ==================================================================
 * L6: Step and Impulse Response
 * ================================================================== */

void fir_step_response(const double *h, int N, double *step, int N_out) {
    if (!h || !step || N < 1 || N_out < 1) return;
    int n;
    double accum = 0.0;
    for (n = 0; n < N_out; n++) {
        int k;
        for (k = 0; k <= n && k < N; k++) {
            accum += h[k];
        }
        step[n] = accum;
        accum = 0.0;
        /* Re-accumulate correctly */
        step[n] = 0.0;
        for (k = 0; k <= n && k < N; k++) {
            step[n] += h[k];
        }
    }
}

void iir_impulse_response(const iir_filter_t *iir, int N_out, double *h) {
    if (!iir || !h || N_out < 1) return;

    filter_state_t *state = (filter_state_t *)malloc(sizeof(filter_state_t));
    if (!state) return;
    memset(state, 0, sizeof(filter_state_t));

    int n;
    for (n = 0; n < N_out; n++) {
        double x = (n == 0) ? 1.0 : 0.0;
        h[n] = iir_process_sample_df2(iir, state, x);
    }

    filter_state_free(state);
}

double step_overshoot(const double *step, int N) {
    if (!step || N < 2) return -1.0;
    double final_val = step[N - 1];
    if (fabs(final_val) < 1e-15) return -1.0;

    double max_val = step[0];
    int i;
    for (i = 1; i < N; i++) {
        if (step[i] > max_val) max_val = step[i];
    }
    return (max_val - final_val) / final_val * 100.0;
}

int settling_time(const double *step, int N, double tolerance) {
    if (!step || N < 2) return -1;
    double final_val = step[N - 1];
    double upper = final_val * (1.0 + tolerance);
    double lower = final_val * (1.0 - tolerance);

    int i;
    for (i = N - 1; i >= 0; i--) {
        if (step[i] > upper || step[i] < lower) {
            return (i < N - 1) ? (i + 1) : -1;
        }
    }
    return 0;
}

int rise_time(const double *step, int N) {
    if (!step || N < 2) return -1;
    double final_val = step[N - 1];
    if (fabs(final_val) < 1e-15) return -1;

    double lo_thresh = 0.10 * final_val;
    double hi_thresh = 0.90 * final_val;

    int t_lo = -1, t_hi = -1;
    int i;
    for (i = 0; i < N; i++) {
        if (t_lo < 0 && step[i] >= lo_thresh) t_lo = i;
        if (t_hi < 0 && step[i] >= hi_thresh) t_hi = i;
    }
    if (t_lo < 0 || t_hi < 0) return -1;
    return t_hi - t_lo;
}

/* ==================================================================
 * L6: Energy and Power Metrics
 * ================================================================== */

double filter_energy(const double *coeff, int length) {
    if (!coeff || length < 1) return -1.0;
    double energy = 0.0;
    int i;
    for (i = 0; i < length; i++) energy += coeff[i] * coeff[i];
    return energy;
}

double filter_l1_norm(const double *coeff, int length) {
    if (!coeff || length < 1) return -1.0;
    double l1 = 0.0;
    int i;
    for (i = 0; i < length; i++) l1 += fabs(coeff[i]);
    return l1;
}

int biquad_l2_scale(biquad_section_t *sections, int num_sections) {
    if (!sections || num_sections < 1) return FILTER_ERR_NULL_PTR;

    int k;
    for (k = 0; k < num_sections; k++) {
        /* L2 norm of impulse response from k-th section to output */
        /* Simplified: compute L2 norm of numerator */
        double l2 = sqrt(sections[k].b0 * sections[k].b0
                       + sections[k].b1 * sections[k].b1
                       + sections[k].b2 * sections[k].b2);
        if (l2 > 1e-15) {
            sections[k].b0 /= l2;
            sections[k].b1 /= l2;
            sections[k].b2 /= l2;
            if (k + 1 < num_sections) {
                sections[k + 1].b0 *= l2;
                sections[k + 1].b1 *= l2;
                sections[k + 1].b2 *= l2;
            }
        }
    }
    return FILTER_OK;
}
