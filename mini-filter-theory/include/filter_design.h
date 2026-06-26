/**
 * filter_design.h — Complete Digital Filter Design Methods
 *
 * L5: FIR design — window, frequency sampling, Parks-McClellan
 * L5: IIR design — bilinear transform, impulse invariance, least-squares
 * L6: Canonical problems — design filters meeting given specifications
 *
 * Courses: MIT 6.003, Berkeley EE123, Illinois ECE 310, CMU 18-396
 * Reference: Oppenheim & Schafer (2010), Ch. 7
 *             Parks & Burrus, Digital Filter Design (1987)
 *             Antoniou, Digital Filters: Analysis, Design, Applications (1993)
 */

#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H

#include "filter_defs.h"
#include "filter_transfer.h"
#include "filter_analog.h"
#include "filter_digital.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L5: Window Function Generation
 * ================================================================== */

/** Generate a window function.
 *  L5: All window functions are real and symmetric about their center.
 *
 *  Rectangular:  w[n] = 1.0
 *  Hann:         w[n] = 0.5*(1 - cos(2*pi*n/(N-1)))
 *  Hamming:      w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1))
 *  Blackman:     w[n] = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1))
 *  Blackman-Harris (4-term):
 *    w[n] = 0.35875 - 0.48829*cos(2*pi*n/(N-1))
 *         + 0.14128*cos(4*pi*n/(N-1)) - 0.01168*cos(6*pi*n/(N-1))
 *  Kaiser:       w[n] = I0(beta*sqrt(1-(2n/(N-1)-1)^2)) / I0(beta)
 *    where I0 is the modified Bessel function of the first kind.
 *  Flat-top:     w[n] = a0 - a1*cos(2*pi*n/(N-1)) + a2*cos(4*pi*n/(N-1))
 *                       - a3*cos(6*pi*n/(N-1)) + a4*cos(8*pi*n/(N-1))
 *  Bartlett:     w[n] = 1 - |2n/(N-1) - 1|
 *  Nuttall (4-term continuous):
 *    w[n] = 0.355768 - 0.487396*cos(2*pi*n/(N-1))
 *         + 0.144232*cos(4*pi*n/(N-1)) - 0.012604*cos(6*pi*n/(N-1))
 *  Tukey:        w[n] = 1 for |n-M| <= alpha*M/2,
 *                0.5*(1+cos(pi*(|n-M|-alpha*M/2)/((1-alpha)*M/2))) elsewhere
 *
 *  Reference: Harris, Proc. IEEE, 66(1):51-83, 1978
 *             Nuttall, NBS Report 8738, 1980
 *  @param win  Output array (length >= params.length)
 *  @param params  Window specification
 *  @return  FILTER_OK or error code
 *  Complexity: O(N) */
int window_generate(double *win, const window_params_t *params);

/** Generate an ideal lowpass impulse response (sinc function).
 *  L5: h_ideal[n] = wc/pi * sinc(wc*(n - (N-1)/2) / pi)
 *  where sinc(x) = sin(pi*x)/(pi*x), and wc = cutoff in rad/sample.
 *  @param h        Output coefficient array (length N)
 *  @param N        Filter length
 *  @param wc       Normalized cutoff frequency (rad/sample, 0 to pi)
 *  Complexity: O(N) */
void ideal_lowpass_impulse(double *h, int N, double wc);

/** Generate an ideal bandpass impulse response.
 *  L5: h_ideal[n] = (w2*sinc(w2*(n-M)) - w1*sinc(w1*(n-M))) / pi
 *  @param h   Output array (length N)
 *  @param N   Filter length
 *  @param w1  Lower cutoff (rad/sample)
 *  @param w2  Upper cutoff (rad/sample) */
void ideal_bandpass_impulse(double *h, int N, double w1, double w2);

/** Generate an ideal highpass impulse response.
 *  L5: h_ideal[n] = delta[n-M] - wc/pi*sinc(wc*(n-M)/pi)
 *  where delta is the Kronecker delta function. */
void ideal_highpass_impulse(double *h, int N, double wc);

/* ==================================================================
 * L5: FIR Design — Window Method (Complete)
 * ================================================================== */

/** Design an FIR filter by windowing the ideal impulse response.
 *  L5/L6: This is the most widely used FIR design method.
 *  Steps: 1) compute h_ideal[n], 2) generate window w[n],
 *         3) h[n] = h_ideal[n] * w[n].
 *
 *  @param spec  Filter specification (fc1, fc2 define passband edges)
 *  @param N     Filter length (number of taps)
 *  @param win_type  Window function
 *  @param win_alpha Window parameter (beta for Kaiser, alpha for Tukey)
 *  @return      Designed FIR filter, or NULL on error */
fir_filter_t *fir_design_by_window(const filter_spec_t *spec, int N,
                                    window_type_t win_type, double win_alpha);

/** Estimate required FIR filter length using Kaiser's formula.
 *  L5: N ≈ (-20*log10(sqrt(dp*ds)) - 13) / (14.6*(ws-wp)/(2*pi))
 *  where dp = passband ripple (linear), ds = stopband ripple (linear),
 *  ws-wp = transition bandwidth (rad/sample).
 *  Reference: Kaiser (1974), simplified by Bellanger
 *  @param transition_bw  Transition bandwidth (rad/sample)
 *  @param atten_db       Minimum stopband attenuation (dB)
 *  @return               Estimated filter length */
int fir_estimate_length(double transition_bw, double atten_db);

/** Design an FIR filter using the frequency sampling method.
 *  L5: Specify desired frequency response at N equally spaced points,
 *  then inverse DFT to get impulse response.
 *  Allows arbitrary magnitude specification, not just LP/HP/BP/BS.
 *  @param desired_mag  Desired magnitude at N/2+1 frequencies [0..fs/2]
 *  @param N            Filter length
 *  @param is_linear_phase  If 1, construct Type I/II linear phase filter
 *  @return             Designed FIR filter */
fir_filter_t *fir_design_freq_sampling(const double *desired_mag, int N,
                                        int is_linear_phase);

/* ==================================================================
 * L5: IIR Design via Analog Prototypes
 * ================================================================== */

/** Design an IIR filter using the bilinear transform method.
 *  L5/L6: Complete digital filter design pipeline:
 *  1. Pre-warp digital cutoff frequencies to analog domain
 *  2. Design analog prototype with specified approximation
 *  3. Apply frequency transformation (LP/HP/BP/BS)
 *  4. Convert to digital via bilinear transform: s = 2*fs*(z-1)/(z+1)
 *  5. Factor into cascade of biquad sections for implementation
 *
 *  @param spec  Filter specification (sample_rate must be > 0)
 *  @return      IIR filter as cascade of biquads, or NULL on error */
iir_filter_t *iir_design_bilinear(const filter_spec_t *spec);

/** Design an IIR filter using the impulse invariance method.
 *  L5: Sample the analog impulse response: h[n] = T * h_a(n*T).
 *  Only suitable for lowpass and bandpass filters (aliasing in other types).
 *  Preserves the impulse response shape exactly at sample instants.
 *  Does NOT use pre-warping — the analog-to-digital frequency mapping
 *  is many-to-one (aliasing).
 *  @param spec  Filter specification
 *  @return      IIR filter as cascade of biquads, or NULL on error */
iir_filter_t *iir_design_impulse_invariance(const filter_spec_t *spec);

/** Convert an IIR filter between implementation structures.
 *  L5: Transform between Direct Form I, II, Transposed, etc.
 *  Note: Cascade to parallel conversion requires partial fraction expansion.
 *  @param iir    Input IIR filter
 *  @param target Target structure
 *  @return       New IIR filter with target structure, or NULL */
iir_filter_t *iir_convert_structure(const iir_filter_t *iir,
                                     digital_structure_t target);

/* ==================================================================
 * L5: IIR Filter Order Estimation
 * ================================================================== */

/** Estimate required IIR filter order for Butterworth.
 *  L5: N >= log10((10^{As/10}-1)/(10^{Ap/10}-1)) / (2*log10(ws/wp)) */
int butterworth_estimate_order(double passband_ripple_db,
                                double stopband_atten_db,
                                double wp, double ws);

/** Estimate required IIR filter order for Bessel.
 *  L5: Approximate formula based on group delay flatness requirement.
 *  Bessel requires the highest order for a given magnitude spec. */
int bessel_estimate_order(double passband_ripple_db,
                           double stopband_atten_db,
                           double wp, double ws);

/** Search for optimal Kaiser window beta parameter.
 *  L5: beta = 0 for As < 21 dB
 *       beta = 0.5842*(As - 21)^0.4 + 0.07886*(As - 21) for 21 <= As <= 50
 *       beta = 0.1102*(As - 8.7) for As > 50
 *  @param atten_db Desired stopband attenuation (dB)
 *  @return         Optimal beta */
double kaiser_beta_from_attenuation(double atten_db);

/** Compute Kaiser window main lobe width as a function of beta.
 *  L5: Main lobe width ≈ (beta*sqrt(pi^2 - ...))/...
 *  Useful for predicting transition bandwidth.
 *  @param beta Kaiser window parameter
 *  @param N    Window length
 *  @return     Approximate main lobe width (rad/sample) */
double kaiser_mainlobe_width(double beta, int N);

/* ==================================================================
 * L5: Least-Squares FIR Design
 * ================================================================== */

/** Design an FIR filter using the least-squares method.
 *  L5: Minimize integral of squared error between desired and actual
 *  frequency response over specified frequency bands.
 *  This is a quadratic optimization problem with a closed-form solution.
 *  Unlike equiripple design, LS minimizes total error energy rather
 *  than maximum error.
 *
 *  @param bands   Frequency band edges array [w0,w1,w2,w3,...] (rad/sample)
 *  @param desired Desired response in each band [g0,g1,...]
 *  @param weights Weight for each band [w0,w1,...]
 *  @param num_bands Number of frequency bands
 *  @param N       Filter length
 *  @return        Least-squares optimal FIR filter */
fir_filter_t *fir_design_least_squares(const double *bands,
                                        const double *desired,
                                        const double *weights,
                                        int num_bands, int N);

/* ==================================================================
 * L5: Quantization Effects Analysis
 * ================================================================== */

/** Quantize filter coefficients to fixed-point representation.
 *  L5: Models the effect of finite wordlength on filter response.
 *  coef_quantized = round(coef * 2^frac_bits) / 2^frac_bits
 *  This changes pole/zero locations and may cause instability in IIR filters.
 *  @param coef        Input coefficient array (modified in-place)
 *  @param num_coefs   Number of coefficients
 *  @param frac_bits   Number of fractional bits
 *  Complexity: O(num_coefs) */
void coef_quantize(double *coef, int num_coefs, int frac_bits);

/** Analyze coefficient sensitivity of an SOS filter.
 *  L5: Compute how pole/zero locations change with coefficient perturbation.
 *  High-Q poles are most sensitive to coefficient errors.
 *  This is why high-order filters MUST use cascade form, not direct form.
 *  @param sos    Second-order section matrix
 *  @param sens   Output sensitivity matrix (same dimensions as SOS)
 *  @return       Maximum pole sensitivity (dB), or -1 on error */
double sos_sensitivity(const sos_matrix_t *sos, double *sens);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_DESIGN_H */
