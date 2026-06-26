/**
 * @file frequency_response.h
 * @brief L2/L5 - Frequency Response Analysis and Filter Characterization
 *
 * Knowledge Coverage:
 *   L2-F01: Magnitude response |H(jw)| computation
 *   L2-F02: Phase response angle(H(jw)) computation
 *   L2-F03: -3 dB bandwidth detection
 *   L2-F04: Resonance frequency and Q factor
 *   L2-F05: Half-power (-3 dB) points
 *   L5-F06: DC gain and high-frequency gain
 *   L5-F07: Roll-off rate (dB/decade)
 *   L5-F08: Cutoff frequency for Butterworth filters
 *   L5-F09: Passband/stopband ripple analysis
 *   L5-F10: Frequency warping (bilinear prewarp)
 *
 * Course Mapping:
 *   MIT 6.003 Ch.6; Stanford EE102A Ch.5;
 *   Berkeley EE123 Ch.5; Tsinghua Signals Ch.3
 */

#ifndef FREQUENCY_RESPONSE_H
#define FREQUENCY_RESPONSE_H

#include "system_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Frequency Response Computation ---- */

int compute_magnitude_response(const ct_transfer_function_t *tf,
                                const double *frequencies_hz, int num_freq,
                                double *magnitude_linear,
                                double *magnitude_db);

int compute_phase_response(const ct_transfer_function_t *tf,
                            const double *frequencies_hz, int num_freq,
                            double *phase_rad, double *phase_deg);

int compute_frequency_response_full(const ct_transfer_function_t *tf,
                                     const double *frequencies_hz, int num_freq,
                                     double *mag_linear, double *mag_db,
                                     double *phase_rad, double *phase_deg);

/* ---- Bandwidth Analysis ---- */

double find_3db_bandwidth(const frequency_response_t *fr);

double find_3db_cutoff_frequency(const frequency_response_t *fr);

int find_half_power_points(const frequency_response_t *fr,
                            double *f_low, double *f_high);

/* ---- Resonance Analysis ---- */

double find_resonance_frequency(const frequency_response_t *fr);

double compute_q_factor(const frequency_response_t *fr);

double compute_q_factor_from_poles(const pole_zero_collection_t *pz);

/* ---- Roll-off Analysis ---- */

double compute_rolloff_db_per_decade(const frequency_response_t *fr,
                                      double f_start, double f_end);

double compute_rolloff_db_per_octave(const frequency_response_t *fr,
                                      double f_start, double f_end);

/* ---- Filter Characterization ---- */

int is_lowpass(const frequency_response_t *fr);

int is_highpass(const frequency_response_t *fr);

int is_bandpass(const frequency_response_t *fr);

int is_bandstop(const frequency_response_t *fr);

double compute_passband_ripple_db(const frequency_response_t *fr,
                                   double f_low, double f_high);

double compute_stopband_attenuation_db(const frequency_response_t *fr,
                                        double f_low, double f_high);

/* ---- DC and HF Gain ---- */

double compute_dc_gain_ct(const ct_transfer_function_t *tf);

double compute_dc_gain_dt(const dt_transfer_function_t *tf);

double compute_hf_gain_ct(const ct_transfer_function_t *tf);

/* ---- Bilinear Prewarping ---- */

double bilinear_prewarp_frequency(double analog_freq_hz,
                                   double sample_rate_hz);

/* ---- Frequency-Domain Convolution ---- */

int frequency_domain_multiply(const double complex *X, const double complex *H,
                               int N, double complex *Y);

#ifdef __cplusplus
}
#endif

#endif /* FREQUENCY_RESPONSE_H */
