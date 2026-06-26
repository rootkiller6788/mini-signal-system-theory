/**
 * filter_digital.h — Digital Filter Structures and Realization
 *
 * L1: FIR and IIR digital filter data structures
 * L2: Direct Form I/II, Transposed, Cascade, Parallel realizations
 * L5: Digital filter implementation — streaming, block processing
 * L6: Canonical: implement and run a digital filter in real time
 *
 * Courses: MIT 6.003, Berkeley EE123, Michigan EECS 351, TU Munich
 * Reference: Oppenheim & Schafer (2010) Discrete-Time Signal Processing
 *             Proakis & Manolakis (2007) Digital Signal Processing
 */

#ifndef FILTER_DIGITAL_H
#define FILTER_DIGITAL_H

#include "filter_defs.h"
#include "filter_transfer.h"
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L1/L2: FIR Filter Implementation
 * ================================================================== */

/** Allocate an FIR filter with given coefficients.
 *  L1: Creates a deep copy of the coefficient array.
 *  Automatically determines linear-phase type from coefficient symmetry.
 *  @param coeff   Coefficient array h[0..length-1]
 *  @param length  Number of taps
 *  @return        Allocated FIR filter, or NULL on error
 *  Complexity: O(N) */
fir_filter_t *fir_filter_alloc(const double *coeff, int length);

/** Free an FIR filter and its internal arrays. */
void fir_filter_free(fir_filter_t *fir);

/** Initialize filter state for streaming FIR processing.
 *  L2: Creates a circular delay line of length N, zero-initialized.
 *  Complexity: O(N) */
filter_state_t *fir_state_alloc(int length);

/** Free filter state. */
void filter_state_free(filter_state_t *state);

/** Reset filter state to zero (clear delay lines).
 *  L2: Prepares filter for a new signal segment without transients. */
void filter_state_reset(filter_state_t *state);

/** Process a single sample through an FIR filter (streaming).
 *  L2: y[n] = sum_{k=0}^{N-1} h[k] * x[n-k]
 *  Uses circular buffer for efficient O(N) per sample.
 *  Complexity: O(N) per sample, O(1) memory update */
double fir_process_sample(const fir_filter_t *fir, filter_state_t *state,
                           double x);

/** Process a block of samples through an FIR filter.
 *  L2: Vectorized block processing.
 *  Complexity: O(N * block_len) */
void fir_process_block(const fir_filter_t *fir, filter_state_t *state,
                        const double *x, double *y, int block_len);

/** Design an FIR filter using the window method.
 *  L5: h[n] = h_ideal[n] * w[n] where h_ideal is the inverse DTFT
 *  of the desired frequency response, and w[n] is a window function.
 *
 *  For lowpass: h_ideal[n] = sin(wc*(n-M/2)) / (pi*(n-M/2))
 *  where M = N-1 and wc is the normalized cutoff frequency.
 *  (At n = M/2: h_ideal = wc/pi)
 *
 *  Reference: Kaiser (1974), "Nonrecursive Digital Filter Design
 *             Using the I_0-sinh Window Function"
 *  @param cutoff      Normalized cutoff frequency (0 to 0.5, fraction of fs)
 *  @param num_taps    Number of taps (filter length N)
 *  @param win_type    Window function to apply
 *  @param win_alpha   Window parameter (Kaiser beta, Tukey alpha)
 *  @return            Designed FIR filter, or NULL on error */
fir_filter_t *fir_design_window(double cutoff, int num_taps,
                                 window_type_t win_type, double win_alpha);

/** Design an FIR differentiator.
 *  L5: h_ideal[n] = cos(pi*(n-M/2))/(n-M/2) - sin(pi*(n-M/2))/(pi*(n-M/2)^2)
 *  for n != M/2, and h_ideal[M/2] = 0.
 *  Type III or IV linear phase.
 *  @param num_taps  Filter length
 *  @param win_type  Window function
 *  @param win_alpha Window parameter
 *  @return          FIR differentiator */
fir_filter_t *fir_design_differentiator(int num_taps, window_type_t win_type,
                                         double win_alpha);

/** Design an FIR Hilbert transformer.
 *  L5: h_ideal[n] = (2/pi)*sin^2(pi*(n-M/2)/2)/(n-M/2) for n-M/2 odd,
 *  0 for n-M/2 even. Type III or IV.
 *  Produces 90-degree phase shift across the band.
 *  @param num_taps  Filter length (must be odd for Type III)
 *  @param win_type  Window function
 *  @return          FIR Hilbert transformer */
fir_filter_t *fir_design_hilbert(int num_taps, window_type_t win_type,
                                  double win_alpha);

/* ==================================================================
 * L1/L2: IIR Filter Implementation
 * ================================================================== */

/** Allocate an IIR filter as cascade of biquad sections.
 *  L1: Deep copies the biquad array.
 *  Complexity: O(num_sections) */
iir_filter_t *iir_filter_alloc(const biquad_section_t *sections,
                                int num_sections, double overall_gain);

/** Free an IIR filter and its internal arrays. */
void iir_filter_free(iir_filter_t *iir);

/** Process a single sample through an IIR filter (DF-II biquad cascade).
 *  L2: For each biquad section k:
 *    w[n] = x[n] - a1k*w[n-1] - a2k*w[n-2]
 *    y[n] = b0k*w[n] + b1k*w[n-1] + b2k*w[n-2]
 *  Uses Direct Form II which minimizes state memory (2 states per biquad).
 *  Complexity: O(num_sections) per sample */
double iir_process_sample_df2(const iir_filter_t *iir, filter_state_t *state,
                                double x);

/** Process a single sample through an IIR filter (DF-I biquad cascade).
 *  L2: For each biquad:
 *    y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
 *  DF-I requires 4 states per biquad but is less sensitive to
 *  limit cycles in fixed-point implementations.
 *  Complexity: O(num_sections) per sample */
double iir_process_sample_df1(const iir_filter_t *iir, filter_state_t *state,
                                double x);

/** Process a block of samples through an IIR filter.
 *  L2: Block processing — cannot be vectorized due to feedback.
 *  Complexity: O(num_sections * block_len) */
void iir_process_block(const iir_filter_t *iir, filter_state_t *state,
                        const double *x, double *y, int block_len);

/* ==================================================================
 * L5: Cascaded Integrator-Comb (CIC) Filter
 * ================================================================== */

/** CIC filter — multiplierless, efficient decimation/interpolation.
 *  L5: H(z) = ((1 - z^{-RM})/(1 - z^{-1}))^N
 *  where R = rate change, M = differential delay (usually 1 or 2),
 *  N = number of stages.
 *  Extremely efficient for high-rate sample rate conversion.
 *  No multipliers needed — only adders and registers.
 *  Reference: Hogenauer, "An Economical Class of Digital Filters for
 *             Decimation and Interpolation" (1981)
 *  @param rate_change  R (decimation or interpolation factor)
 *  @param num_stages   N (number of integrator-comb pairs)
 *  @param diff_delay   M (differential delay, typically 1)
 *  @return             Combined frequency response magnitude at each freq */
double cic_magnitude(double freq_norm, int rate_change, int num_stages,
                      int diff_delay);

/** Compute CIC compensation FIR filter coefficients.
 *  L5: Compensates for the sinc^N passband droop of CIC filters.
 *  The compensation filter has response ~ 1/sinc^N(f) in the passband.
 *  Usually implemented as a short FIR after CIC decimation.
 *  @param rate_change  CIC rate change factor R
 *  @param num_stages   CIC number of stages N
 *  @param diff_delay   CIC differential delay M
 *  @param comp_taps    Number of compensation FIR taps
 *  @return             FIR compensation filter */
fir_filter_t *cic_compensation_filter(int rate_change, int num_stages,
                                       int diff_delay, int comp_taps);

/* ==================================================================
 * L2: Half-Band Filter
 * ================================================================== */

/** Design a half-band FIR filter.
 *  L2: Half-band filters have fc = fs/4 and every other coefficient
 *  (except the center tap) is zero. This reduces computational cost by ~50%.
 *  h[2k] = 0 for all k != (N-1)/2.
 *  Useful for efficient decimation/interpolation by 2.
 *  The passband and stopband are symmetric about fs/4.
 *  @param num_taps  Number of taps (must be odd, = 4K+3 for some K)
 *  @param win_type  Window function
 *  @return          Half-band FIR filter */
fir_filter_t *fir_design_halfband(int num_taps, window_type_t win_type);

/** Check if an FIR filter has half-band symmetry.
 *  Verifies h[2k] ≈ 0 for k != (N-1)/2.
 *  @param tolerance  Maximum allowed deviation from zero for even-indexed taps
 *  @return           1 if half-band, 0 otherwise */
int fir_is_halfband(const fir_filter_t *fir, double tolerance);

/* ==================================================================
 * L2: Linear Phase Analysis
 * ================================================================== */

/** Determine the linear-phase type of FIR coefficients.
 *  L2: Checks symmetry and anti-symmetry:
 *    even symmetry + odd length  -> Type I
 *    even symmetry + even length -> Type II
 *    odd symmetry + odd length   -> Type III
 *    odd symmetry + even length  -> Type IV
 *  @param coeff   Coefficient array
 *  @param length  Number of taps
 *  @param tol     Tolerance for symmetry check
 *  @return        FIR linear phase type, or -1 if not linear phase */
int fir_determine_linear_phase_type(const double *coeff, int length,
                                     double tol);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_DIGITAL_H */
