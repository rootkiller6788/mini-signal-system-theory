/**
 * @file frequency_response.c
 * @brief L2/L5 Frequency Response Analysis and Filter Characterization
 *
 * Knowledge Implementation (14 points):
 *   L2-I01: compute_magnitude_response - |H(jw)| and dB
 *   L2-I02: compute_phase_response - angle(H(jw)) 
 *   L2-I03: compute_frequency_response_full - combined mag+phase
 *   L2-I04: find_3db_bandwidth - -3dB bandwidth detection
 *   L2-I05: find_3db_cutoff_frequency
 *   L2-I06: find_half_power_points - f_low and f_high
 *   L2-I07: find_resonance_frequency - peak in |H(jw)|
 *   L2-I08: compute_q_factor - quality factor from bandwidth
 *   L2-I09: compute_q_factor_from_poles - Q from pole damping
 *   L2-I10: compute_rolloff_db_per_decade - asymptotic slope
 *   L2-I11: compute_rolloff_db_per_octave
 *   L2-I12: is_lowpass / is_highpass / is_bandpass / is_bandstop
 *   L2-I13: compute_passband_ripple_db / compute_stopband_attenuation_db
 *   L2-I14: compute_dc_gain_ct / compute_dc_gain_dt / compute_hf_gain_ct
 *   L2-I15: bilinear_prewarp_frequency
 *   L2-I16: frequency_domain_multiply
 *
 * References:
 *   Sedra & Smith, Microelectronic Circuits, 8th Ed.
 *   Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd Ed.
 */

#include "frequency_response.h"
#include "transfer_function.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 *  L2-I01: Compute magnitude response |H(jw)|
 * ======================================================================== */

int compute_magnitude_response(const ct_transfer_function_t *tf,
                                const double *frequencies_hz, int num_freq,
                                double *magnitude_linear,
                                double *magnitude_db)
{
    if (!tf || !frequencies_hz || num_freq <= 0) return -1;

    for (int i = 0; i < num_freq; i++) {
        double w = 2.0 * M_PI * frequencies_hz[i];
        double complex s = 0.0 + w * I;
        double complex H = ct_tf_evaluate(tf, s);
        double mag = cabs(H);

        if (magnitude_linear) magnitude_linear[i] = mag;
        if (magnitude_db) magnitude_db[i] = (mag > 1e-30) ? 20.0 * log10(mag) : -300.0;
    }
    return 0;
}

/* ========================================================================
 *  L2-I02: Compute phase response angle(H(jw))
 * ======================================================================== */

int compute_phase_response(const ct_transfer_function_t *tf,
                            const double *frequencies_hz, int num_freq,
                            double *phase_rad, double *phase_deg)
{
    if (!tf || !frequencies_hz || num_freq <= 0) return -1;

    for (int i = 0; i < num_freq; i++) {
        double w = 2.0 * M_PI * frequencies_hz[i];
        double complex s = 0.0 + w * I;
        double complex H = ct_tf_evaluate(tf, s);
        double ph = carg(H);

        if (phase_rad) phase_rad[i] = ph;
        if (phase_deg) phase_deg[i] = ph * 180.0 / M_PI;
    }
    return 0;
}

/* ========================================================================
 *  L2-I03: Full frequency response (magnitude + phase)
 * ======================================================================== */

int compute_frequency_response_full(const ct_transfer_function_t *tf,
                                     const double *frequencies_hz, int num_freq,
                                     double *mag_linear, double *mag_db,
                                     double *phase_rad, double *phase_deg)
{
    if (!tf || !frequencies_hz || num_freq <= 0) return -1;

    for (int i = 0; i < num_freq; i++) {
        double w = 2.0 * M_PI * frequencies_hz[i];
        double complex s = 0.0 + w * I;
        double complex H = ct_tf_evaluate(tf, s);
        double mag = cabs(H);
        double ph  = carg(H);

        if (mag_linear) mag_linear[i] = mag;
        if (mag_db)     mag_db[i]     = (mag > 1e-30) ? 20.0 * log10(mag) : -300.0;
        if (phase_rad)  phase_rad[i]  = ph;
        if (phase_deg)  phase_deg[i]  = ph * 180.0 / M_PI;
    }
    return 0;
}

/* ========================================================================
 *  L2-I04: Find -3 dB bandwidth
 *
 *  Bandwidth is the frequency range where |H(jw)| >= max|H| / sqrt(2)
 *  i.e., no more than 3 dB below the maximum.
 * ======================================================================== */

double find_3db_bandwidth(const frequency_response_t *fr)
{
    if (!fr || !fr->points || fr->num_points < 2) return -1.0;

    /* Find max magnitude */
    double max_mag_db = -INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].magnitude_db > max_mag_db)
            max_mag_db = fr->points[i].magnitude_db;
    }

    double threshold_db = max_mag_db - 3.0;

    /* Find first crossing from below */
    double f_low = -1.0, f_high = -1.0;
    for (int i = 1; i < fr->num_points; i++) {
        double prev = fr->points[i-1].magnitude_db;
        double cur  = fr->points[i].magnitude_db;

        if (f_low < 0 && prev < threshold_db && cur >= threshold_db) {
            f_low = fr->points[i].frequency;
        }
        if (f_low > 0 && prev >= threshold_db && cur < threshold_db) {
            f_high = fr->points[i].frequency;
            break;
        }
    }

    if (f_low > 0 && f_high > 0)
        return f_high - f_low;
    else if (f_high > 0)
        return f_high;  /* Lowpass: bandwidth = cutoff frequency */
    return -1.0;
}

/* ========================================================================
 *  L2-I05: Find -3 dB cutoff frequency
 * ======================================================================== */

double find_3db_cutoff_frequency(const frequency_response_t *fr)
{
    if (!fr || !fr->points || fr->num_points < 2) return -1.0;

    double max_mag_db = -INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].magnitude_db > max_mag_db)
            max_mag_db = fr->points[i].magnitude_db;
    }

    double threshold_db = max_mag_db - 3.0;

    /* For lowpass: find first crossing from above */
    for (int i = 1; i < fr->num_points; i++) {
        if (fr->points[i-1].magnitude_db >= threshold_db &&
            fr->points[i].magnitude_db < threshold_db) {
            /* Linear interpolation */
            double f1 = fr->points[i-1].frequency;
            double f2 = fr->points[i].frequency;
            double m1 = fr->points[i-1].magnitude_db;
            double m2 = fr->points[i].magnitude_db;
            if (fabs(m2 - m1) > 1e-15) {
                return f1 + (f2 - f1) * (threshold_db - m1) / (m2 - m1);
            }
            return f2;
        }
    }
    return fr->points[fr->num_points - 1].frequency;  /* Last resort */
}

/* ========================================================================
 *  L2-I06: Find half-power (-3 dB) points
 * ======================================================================== */

int find_half_power_points(const frequency_response_t *fr,
                            double *f_low, double *f_high)
{
    if (!fr || !f_low || !f_high) return -1;

    double max_mag_db = -INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].magnitude_db > max_mag_db)
            max_mag_db = fr->points[i].magnitude_db;
    }

    double threshold = max_mag_db - 3.0;
    *f_low = -1.0; *f_high = -1.0;

    for (int i = 1; i < fr->num_points; i++) {
        double prev = fr->points[i-1].magnitude_db;
        double cur  = fr->points[i].magnitude_db;

        if (*f_low < 0 && prev < threshold && cur >= threshold)
            *f_low = fr->points[i].frequency;
        if (*f_low > 0 && prev >= threshold && cur < threshold) {
            *f_high = fr->points[i].frequency;
            return 0;
        }
    }
    return (*f_low > 0) ? 0 : -1;
}

/* ========================================================================
 *  L2-I07: Find resonance frequency (peak in magnitude)
 * ======================================================================== */

double find_resonance_frequency(const frequency_response_t *fr)
{
    if (!fr || !fr->points || fr->num_points < 2) return -1.0;

    double max_mag = -INFINITY;
    int max_idx = 0;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].magnitude > max_mag) {
            max_mag = fr->points[i].magnitude;
            max_idx = i;
        }
    }
    return fr->points[max_idx].frequency;
}

/* ========================================================================
 *  L2-I08: Q Factor from -3 dB bandwidth
 *
 *  Q = f_resonance / BW  (for bandpass filters)
 *  Q = 1 / (2 * zeta)    (for second-order lowpass)
 * ======================================================================== */

double compute_q_factor(const frequency_response_t *fr)
{
    if (!fr) return -1.0;

    double f_res = find_resonance_frequency(fr);
    if (f_res <= 0.0) return -1.0;

    double bw = find_3db_bandwidth(fr);
    if (bw <= 0.0) return -1.0;

    return f_res / bw;
}

/* ========================================================================
 *  L2-I09: Q Factor from dominant pole pair damping
 * ======================================================================== */

double compute_q_factor_from_poles(const pole_zero_collection_t *pz)
{
    if (!pz || pz->num_poles < 2) return -1.0;

    /* For second-order: Q = 1/(2*zeta) */
    /* Find complex conjugate pair */
    double min_damping = INFINITY;
    for (int i = 0; i < pz->num_poles; i++) {
        if (pz->poles[i].damping < min_damping && pz->poles[i].damping > 1e-10) {
            min_damping = pz->poles[i].damping;
        }
    }

    if (min_damping >= INFINITY) return -1.0;
    if (min_damping < 1e-10) return INFINITY;

    return 1.0 / (2.0 * min_damping);
}

/* ========================================================================
 *  L2-I10: Roll-off rate in dB/decade
 * ======================================================================== */

double compute_rolloff_db_per_decade(const frequency_response_t *fr,
                                      double f_start, double f_end)
{
    if (!fr || !fr->points || fr->num_points < 2) return -1.0;

    int i_start = -1, i_end = -1;
    for (int i = 0; i < fr->num_points; i++) {
        if (i_start < 0 && fr->points[i].frequency >= f_start) i_start = i;
        if (fr->points[i].frequency <= f_end) i_end = i;
    }

    if (i_start < 0 || i_end < 0 || i_start >= i_end) return 0.0;

    double mag_start = fr->points[i_start].magnitude_db;
    double mag_end   = fr->points[i_end].magnitude_db;
    double decades   = log10(f_end / f_start);

    if (decades < 1e-15) return 0.0;
    return (mag_end - mag_start) / decades;
}

/* ========================================================================
 *  L2-I11: Roll-off rate in dB/octave
 * ======================================================================== */

double compute_rolloff_db_per_octave(const frequency_response_t *fr,
                                      double f_start, double f_end)
{
    double db_per_decade = compute_rolloff_db_per_decade(fr, f_start, f_end);
    /* 1 decade = log2(10) ≈ 3.322 octaves */
    return db_per_decade / 3.321928094887362;
}

/* ========================================================================
 *  L2-I12: Filter type classification
 * ======================================================================== */

int is_lowpass(const frequency_response_t *fr)
{
    if (!fr || fr->num_points < 3) return 0;
    double mag_first = fr->points[0].magnitude_db;
    double mag_last  = fr->points[fr->num_points - 1].magnitude_db;
    /* Lowpass: high at DC, rolls off at HF */
    return (mag_first > mag_last + 10.0) ? 1 : 0;
}

int is_highpass(const frequency_response_t *fr)
{
    if (!fr || fr->num_points < 3) return 0;
    double mag_first = fr->points[0].magnitude_db;
    double mag_last  = fr->points[fr->num_points - 1].magnitude_db;
    return (mag_last > mag_first + 10.0) ? 1 : 0;
}

int is_bandpass(const frequency_response_t *fr)
{
    if (!fr || fr->num_points < 3) return 0;
    /* Bandpass: peaks in the middle, low at both ends */
    double mag_first = fr->points[0].magnitude_db;
    double mag_last  = fr->points[fr->num_points - 1].magnitude_db;
    double mag_mid   = fr->points[fr->num_points / 2].magnitude_db;
    return (mag_mid > mag_first + 3.0 && mag_mid > mag_last + 3.0) ? 1 : 0;
}

int is_bandstop(const frequency_response_t *fr)
{
    if (!fr || fr->num_points < 3) return 0;
    double mag_first = fr->points[0].magnitude_db;
    double mag_last  = fr->points[fr->num_points - 1].magnitude_db;
    double mag_mid   = fr->points[fr->num_points / 2].magnitude_db;
    return (mag_mid < mag_first - 3.0 && mag_mid < mag_last - 3.0) ? 1 : 0;
}

/* ========================================================================
 *  L2-I13: Passband ripple and stopband attenuation
 * ======================================================================== */

double compute_passband_ripple_db(const frequency_response_t *fr,
                                   double f_low, double f_high)
{
    if (!fr || fr->num_points < 2) return -1.0;

    double max_db = -INFINITY, min_db = INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].frequency >= f_low &&
            fr->points[i].frequency <= f_high) {
            if (fr->points[i].magnitude_db > max_db)
                max_db = fr->points[i].magnitude_db;
            if (fr->points[i].magnitude_db < min_db)
                min_db = fr->points[i].magnitude_db;
        }
    }
    return (max_db > -INFINITY) ? (max_db - min_db) : -1.0;
}

double compute_stopband_attenuation_db(const frequency_response_t *fr,
                                        double f_low, double f_high)
{
    if (!fr || fr->num_points < 2) return -1.0;

    /* Find max in passband (below f_low for lowpass) */
    double max_passband = -INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].frequency < f_low) {
            if (fr->points[i].magnitude_db > max_passband)
                max_passband = fr->points[i].magnitude_db;
        }
    }

    double min_stopband = INFINITY;
    for (int i = 0; i < fr->num_points; i++) {
        if (fr->points[i].frequency >= f_low &&
            fr->points[i].frequency <= f_high) {
            if (fr->points[i].magnitude_db < min_stopband)
                min_stopband = fr->points[i].magnitude_db;
        }
    }

    return (max_passband > -INFINITY && min_stopband < INFINITY)
        ? (max_passband - min_stopband) : -1.0;
}

/* ========================================================================
 *  L2-I14: DC and High-Frequency Gain
 * ======================================================================== */

double compute_dc_gain_ct(const ct_transfer_function_t *tf)
{
    if (!tf) return 0.0;
    return creal(ct_tf_evaluate(tf, 0.0 + 0.0 * I));
}

double compute_dc_gain_dt(const dt_transfer_function_t *tf)
{
    if (!tf) return 0.0;
    return creal(dt_tf_evaluate(tf, 1.0 + 0.0 * I));
}

double compute_hf_gain_ct(const ct_transfer_function_t *tf)
{
    if (!tf) return 0.0;
    /* lim_{s->inf} H(s): depends on relative degree */
    if (tf->num_order == tf->den_order)
        return tf->num_coeffs[tf->num_order] / tf->den_coeffs[tf->den_order];
    if (tf->num_order > tf->den_order) return INFINITY;
    return 0.0;
}

/* ========================================================================
 *  L2-I15: Bilinear Prewarping
 *
 *  Digital frequency w_d = 2*arctan(w_a * T / 2) / T
 *  Analog frequency  w_a = (2/T) * tan(w_d * T / 2)
 *
 *  Use prewarping to match a specific analog frequency to its
 *  digital counterpart when designing digital filters via bilinear.
 * ======================================================================== */

double bilinear_prewarp_frequency(double analog_freq_hz,
                                   double sample_rate_hz)
{
    if (sample_rate_hz <= 0.0) return analog_freq_hz;

    double T = 1.0 / sample_rate_hz;
    double w_a = 2.0 * M_PI * analog_freq_hz;
    double w_d = 2.0 * atan(w_a * T / 2.0) / T;
    return w_d / (2.0 * M_PI);
}

/* ========================================================================
 *  L2-I16: Frequency-Domain Multiplication (for convolution)
 *
 *  Y[k] = X[k] * H[k] for each frequency bin k.
 *  Used in FFT-based filtering.
 * ======================================================================== */

int frequency_domain_multiply(const double complex *X, const double complex *H,
                               int N, double complex *Y)
{
    if (!X || !H || !Y || N <= 0) return -1;

    for (int k = 0; k < N; k++) {
        Y[k] = X[k] * H[k];
    }
    return 0;
}
