/**
 * filter_analysis.h — Frequency Response, Stability, and Performance Analysis
 *
 * L2: Frequency response computation, Bode/Nyquist/Nichols plots
 * L3: Group delay, phase delay, pole-zero analysis
 * L4: Stability margins, Routh-Hurwitz, Jury test
 * L6: Filter performance metrics: SNR improvement, THD, passband flatness
 *
 * Courses: MIT 6.003, Stanford EE102A, Berkeley EE123, ETH 227-0427
 * Reference: Oppenheim & Willsky (1997) Signals and Systems
 *             Lyons, Understanding Digital Signal Processing (2010)
 */

#ifndef FILTER_ANALYSIS_H
#define FILTER_ANALYSIS_H

#include "filter_defs.h"
#include "filter_transfer.h"
#include "filter_digital.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L2: Frequency Response Computation
 * ================================================================== */

/** Compute the frequency response of an FIR filter at a single frequency.
 *  L2: H(e^{jw}) = sum_{n=0}^{N-1} h[n] * e^{-j*w*n}
 *  Uses direct evaluation (not FFT-based) for accuracy at arbitrary w.
 *  Complexity: O(N) per frequency point */
freq_resp_point_t fir_freq_response_at(const fir_filter_t *fir, double w);

/** Compute the frequency response of an IIR filter at a single frequency.
 *  L2: Evaluates each biquad: H_k(e^{jw}) = (b0 + b1*e^{-jw} + b2*e^{-j2w})
 *                                            / (1 + a1*e^{-jw} + a2*e^{-j2w})
 *  Then multiplies magnitudes and adds phases.
 *  Complexity: O(num_sections) per frequency point */
freq_resp_point_t iir_freq_response_at(const iir_filter_t *iir, double w);

/** Compute the complete frequency response vector for an FIR filter.
 *  L2/L3: Evaluates H(e^{jw}) at num_points logarithmically or
 *  linearly spaced frequencies between f_start and f_stop.
 *  Also computes group delay via numerical differentiation.
 *  Complexity: O(N * num_points) */
freq_resp_t *fir_freq_response(const fir_filter_t *fir, double f_start,
                                double f_stop, int num_points, double fs);

/** Compute the complete frequency response vector for an IIR filter. */
freq_resp_t *iir_freq_response(const iir_filter_t *iir, double f_start,
                                double f_stop, int num_points, double fs);

/** Free a frequency response structure. */
void freq_resp_free(freq_resp_t *resp);

/** Compute frequency response of a continuous-time transfer function.
 *  L2: H(jw) = N(jw)/D(jw) where w is real frequency in rad/s.
 *  Complexity: O(max(num_len, den_len)) per point */
freq_resp_point_t tf_continuous_response_at(const tf_continuous_t *tf,
                                              double omega_rad_s);

/** Compute complete frequency response of continuous-time TF. */
freq_resp_t *tf_continuous_freq_response(const tf_continuous_t *tf,
                                           double f_start, double f_stop,
                                           int num_points);

/* ==================================================================
 * L2: Filter Performance Metrics
 * ================================================================== */

/** Compute the passband flatness (maximum deviation in dB).
 *  L2: Measures how much the magnitude varies within the passband.
 *  A perfectly flat filter has flatness = 0 dB.
 *  @param resp  Frequency response
 *  @param f_pass_start  Passband start frequency (Hz)
 *  @param f_pass_stop   Passband stop frequency (Hz)
 *  @return      Maximum deviation from average passband gain in dB */
double passband_flatness(const freq_resp_t *resp,
                          double f_pass_start, double f_pass_stop);

/** Compute the minimum stopband attenuation.
 *  L2: Minimum rejection in the stopband region (most critical point).
 *  @return Minimum attenuation in dB (positive number) */
double stopband_attenuation(const freq_resp_t *resp,
                             double f_stop_start, double f_stop_stop);

/** Find the -3 dB cutoff frequency of a filter.
 *  L2: Binary search for frequency where magnitude drops to
 *  max_magnitude/sqrt(2) -> -3.01 dB from passband reference.
 *  @param resp   Frequency response data
 *  @param ref_db Reference level (0 for DC gain, or passband average)
 *  @return       -3dB frequency in Hz, or -1 if not found */
double find_3db_cutoff(const freq_resp_t *resp, double ref_db);

/** Compute the transition bandwidth (stopband edge - passband edge).
 *  L2: Determined from the magnitude response.
 *  @return Transition bandwidth in Hz, or -1 on error */
double transition_bandwidth(const freq_resp_t *resp,
                             double passband_ripple_db,
                             double stopband_atten_db);

/** Compute the effective noise bandwidth (ENB) of a filter.
 *  L2: ENB = integral(|H(f)|^2 df) / |H(f0)|^2
 *  Measures the equivalent rectangular bandwidth for noise power calculations.
 *  @return ENB in Hz */
double effective_noise_bandwidth(const freq_resp_t *resp);

/* ==================================================================
 * L3: Group Delay and Phase Analysis
 * ================================================================== */

/** Compute group delay at a single frequency via numerical differentiation.
 *  L3: tau_g(w) = -d(phi(w))/dw ≈ -(phi(w+dw) - phi(w-dw))/(2*dw)
 *  Group delay measures the delay of the envelope of a signal.
 *  Constant group delay = linear phase = no waveform distortion.
 *  Complexity: O(N) per point (2 frequency evaluations) */
double group_delay_at_fir(const fir_filter_t *fir, double w, double dw);

/** Compute group delay for IIR filter at a single frequency. */
double group_delay_at_iir(const iir_filter_t *iir, double w, double dw);

/** Compute phase delay at a single frequency.
 *  L3: tau_p(w) = -phi(w)/w
 *  Phase delay measures the delay of the carrier (individual sinusoid).
 *  For linear-phase FIR: phase_delay = group_delay = (N-1)/2 samples. */
double phase_delay_at_fir(const fir_filter_t *fir, double w);

/** Compute phase delay for IIR filter at a single frequency. */
double phase_delay_at_iir(const iir_filter_t *iir, double w);

/** Compute group delay ripple (peak-to-peak variation in passband).
 *  L3: Ideally zero for linear-phase filters.
 *  Non-zero group delay ripple causes dispersion — different
 *  frequencies experience different delays.
 *  @return Peak-to-peak group delay variation in seconds */
double group_delay_ripple(const freq_resp_t *resp,
                           double f_start, double f_stop);

/* ==================================================================
 * L3: Pole-Zero Analysis
 * ================================================================== */

/** Compute pole Q-factors from pole-zero map.
 *  L3: Q = |pole|/(2*|Re(pole)|) for continuous-time
 *  Q = |pole|/(2*|Re(ln(pole))|) for discrete-time (approximate)
 *  Higher Q means sharper resonance, longer ringing, higher sensitivity.
 *  @param pz    Pole-zero map
 *  @param qvals Output array of Q values (length >= num_pairs)
 *  @return      Maximum Q value (most critical section) */
double pz_max_q(const pz_map_t *pz, double *qvals);

/** Check if a continuous-time pole-zero map is minimum-phase.
 *  L3: A system is minimum-phase if all poles AND zeros are in the
 *  left half-plane. Minimum-phase systems have the smallest possible
 *  group delay for a given magnitude response.
 *  @return 1 if minimum-phase, 0 otherwise */
int pz_is_minimum_phase(const pz_map_t *pz);

/** Compute the pole-zero plot distance from stability boundary.
 *  L3: For continuous-time: minimum |Re(pole)| = stability margin
 *  For discrete-time: minimum (1 - |pole|) = stability margin
 *  Larger margins indicate more robust stability.
 *  @return Minimum distance (positive = stable, negative = unstable) */
double pz_stability_margin(const pz_map_t *pz, int is_discrete);

/* ==================================================================
 * L4: Complete Stability Analysis
 * ================================================================== */

/** Check if pole locations indicate BIBO stability.
 *  L4: BIBO (Bounded-Input Bounded-Output) stability:
 *    Continuous-time: all poles must have Re(p) < 0
 *    Discrete-time:   all poles must have |p| < 1
 *  This is a necessary and sufficient condition for LTI systems.
 *  @param poles     Array of pole locations
 *  @param num_poles Number of poles
 *  @param is_discrete 1 for discrete-time, 0 for continuous-time
 *  @return 1 if stable, 0 if unstable */
int poles_are_stable(const double complex *poles, int num_poles,
                      int is_discrete);

/** Compute the spectral factorization of a magnitude-squared function.
 *  L4: Given |H(e^{jw})|^2, find a stable, causal H(z) such that
 *  |H(e^{jw})|^2 = H(e^{jw}) * H*(e^{jw}).
 *  This is the basis of many optimal filter design methods.
 *  @param mag_sq_coeff  Coefficients of the magnitude-squared polynomial
 *  @param half_len      Half-length (positive-lag) coefficients count
 *  @param H             Output stable spectral factor (minimum-phase)
 *  @return              FILTER_OK or error code */
int spectral_factor(const double *mag_sq_coeff, int half_len,
                     double *H);

/* ==================================================================
 * L6: Filter Performance Comparison
 * ================================================================== */

/** Compare two frequency responses and compute the difference.
 *  L6: Useful for comparing designed filter vs ideal specification.
 *  @return RMS error between the two magnitude responses in dB */
double freq_resp_compare(const freq_resp_t *resp1,
                          const freq_resp_t *resp2);

/** Compute the SNR improvement provided by a filter.
 *  L6: SNR_out/SNR_in = integral(|H(f)|^2 * S_signal(f) df)
 *                      / integral(|H(f)|^2 * S_noise(f) df)
 *  Assumes white noise (constant noise PSD).
 *  @param resp         Filter frequency response
 *  @param signal_band_low   Lower edge of signal band (Hz)
 *  @param signal_band_high  Upper edge of signal band (Hz)
 *  @return SNR improvement in dB */
double snr_improvement(const freq_resp_t *resp,
                        double signal_band_low, double signal_band_high);

/** Compute total harmonic distortion after filtering.
 *  L6: Measures how much a filter attenuates harmonics relative to
 *  the fundamental. Important for audio and power electronics.
 *  @param resp           Filter frequency response
 *  @param fundamental_hz Fundamental frequency
 *  @param num_harmonics  Number of harmonics to consider
 *  @return THD in dBc (dB relative to carrier after filtering) */
double filter_thd(const freq_resp_t *resp, double fundamental_hz,
                   int num_harmonics);

/* ==================================================================
 * L6: Step and Impulse Response
 * ================================================================== */

/** Compute the step response of an FIR filter.
 *  L6: s[n] = sum_{k=0}^{n} h[k] for n = 0, 1, ..., N-1
 *  (Assumes initial rest: y[-1] = 0)
 *  Step response reveals overshoot, settling time, and rise time.
 *  @param h      FIR coefficients
 *  @param N      Filter length
 *  @param step   Output step response (length >= N)
 *  @param N_out  Number of output samples to compute */
void fir_step_response(const double *h, int N, double *step, int N_out);

/** Compute the impulse response of an IIR filter (truncated).
 *  L6: Computes h[n] by filtering a unit impulse delta[n].
 *  The response is truncated after N_out samples.
 *  @param iir    IIR filter
 *  @param N_out  Number of output samples
 *  @param h      Output impulse response (length >= N_out) */
void iir_impulse_response(const iir_filter_t *iir, int N_out, double *h);

/** Compute overshoot of step response.
 *  L2: Overshoot (%) = (max(step) - step[inf]) / step[inf] * 100
 *  @return Overshoot as a percentage */
double step_overshoot(const double *step, int N);

/** Compute settling time of step response.
 *  L2: Time to settle within tolerance band of final value.
 *  @param step       Step response
 *  @param N          Length of step response
 *  @param tolerance  Tolerance band (e.g., 0.01 for 1%)
 *  @return           Settling time in samples, or -1 if not settled */
int settling_time(const double *step, int N, double tolerance);

/** Compute rise time (10% to 90%) of step response.
 *  L2: Standard metric for speed of response.
 *  @return Rise time in samples, or -1 if cannot compute */
int rise_time(const double *step, int N);

/* ==================================================================
 * L6: Energy and Power Metrics
 * ================================================================== */

/** Compute the energy of filter coefficients (sum of squares).
 *  L3: E = sum h[n]^2 (Parseval: also equals integral of |H(f)|^2 df)
 *  Relevant for scaling to prevent overflow in fixed-point. */
double filter_energy(const double *coeff, int length);

/** Compute L1-norm of impulse response.
 *  L4: L1 = sum |h[n]| — bounds the peak output magnitude.
 *  |y[n]| <= L1 * max|x[n]| (BIBO gain bound) */
double filter_l1_norm(const double *coeff, int length);

/** Compute L2-norm scaling factor for cascade biquads.
 *  L2: Scale each section such that overflow probability is minimized.
 *  Uses L2-norm (Parseval) scaling: gain_k = 1 / L2_k.
 *  Reference: Jackson, "On the Interaction of Roundoff Noise and
 *             Dynamic Range in Digital Filters" (1970) */
int biquad_l2_scale(biquad_section_t *sections, int num_sections);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_ANALYSIS_H */
