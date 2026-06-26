/**
 * @file    fourier_window.h
 * @brief   Window Functions for Spectral Analysis — L2 Concepts + L5 Algorithms
 *
 * @details Comprehensive library of window functions used in spectral analysis
 *          to mitigate spectral leakage caused by finite-length signal truncation.
 *          Each window is characterized by its mainlobe width, sidelobe level,
 *          coherent gain, equivalent noise bandwidth (ENBW), and scalloping loss.
 *
 *          The truncation of a signal x[n] to N samples is multiplication by
 *          a rectangular window w_R[n] = 1 for 0 ≤ n < N.  The DFT of w_R is the
 *          Dirichlet kernel (periodic sinc), whose sidelobes decay as 1/f — only
 *          -13.3 dB for the first sidelobe.  Better windows offer faster sidelobe
 *          decay at the cost of wider mainlobe (reduced frequency resolution).
 *
 * Knowledge Mapping:
 *   L2 - Core Concepts:
 *     - Spectral leakage phenomenon (finite observation window)
 *     - Window function design trade-off: mainlobe width ↔ sidelobe level
 *     - Coherent gain, ENBW, processing loss, scalloping loss
 *     - Overlap processing (Welch's method)
 *   L5 - Algorithms:
 *     - Hann (Hanning) window:     w[n] = 0.5·(1 - cos(2π·n/(N-1)))
 *     - Hamming window:           w[n] = 0.54 - 0.46·cos(2π·n/(N-1))
 *     - Blackman window:          w[n] = 0.42 - 0.5·cos(...) + 0.08·cos(...)
 *     - Nuttall window (4-term):  minimum sidelobe for given order
 *     - Blackman-Harris window (4-term, 7-term)
 *     - Flat-top window:          amplitude-accurate (±0.01 dB)
 *     - Kaiser window:            adjustable β parameter, optimal for given sidelobe
 *     - Bartlett (triangular) window
 *     - Gaussian window
 *     - Dolph-Chebyshev window:   equiripple sidelobes, min mainlobe for given level
 *     - Tukey (cosine-tapered) window
 *     - Lanczos window:           sinc-based
 *   L4 - Theorems:
 *     - Fourier transform of rectangular window = Dirichlet kernel
 *     - Uncertainty principle:  Δt · Δf ≥ 1/(4π)  (Gaussian minimizes product)
 *
 * Reference: Harris, F.J., "On the Use of Windows for Harmonic Analysis
 *            with the Discrete Fourier Transform", Proc. IEEE 66(1), 1978.
 *            Nuttall, A.H., "Some Windows with Very Good Sidelobe Behavior",
 *            IEEE Trans. ASSP 29(1), 1981.
 */

#ifndef FOURIER_WINDOW_H
#define FOURIER_WINDOW_H

#include "fourier_core.h"

/* ─── L5: Window Generation Functions ───────────────────────────────────── */

/**
 * @brief  Generate a rectangular window of length N.
 *
 *         w[n] = 1.0  for n = 0..N-1
 *
 *         Properties: mainlobe 3dB width = 0.89 bins, highest sidelobe = -13.3 dB,
 *         sidelobe decay rate = 6 dB/octave, coherent gain = 1.0, ENBW = 1.0 bins.
 *
 * @param  N  Window length
 * @return    Window config with coefficient array populated
 */
window_config_t  window_rectangular(size_t N);

/**
 * @brief  Generate a Hann (Hanning) window:  w[n] = 0.5·(1 - cos(2π·n/(N-1)))
 *
 *         Sidelobe decay: 18 dB/octave.  Highest sidelobe: -31.5 dB.
 *         Coherent gain: 0.5.  ENBW: 1.5 bins.  Good general-purpose window.
 *
 * @param  N  Window length (N ≥ 2)
 * @return    Window config
 */
window_config_t  window_hann(size_t N);

/**
 * @brief  Generate a Hamming window:  w[n] = 0.54 - 0.46·cos(2π·n/(N-1))
 *
 *         Optimized to place first sidelobe null at the same position as
 *         the Hann window but with lower first sidelobe (-42.6 dB).
 *         Sidelobe decay: 6 dB/octave (slow beyond first sidelobe).
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_hamming(size_t N);

/**
 * @brief  Generate a Blackman window (3-term):
 *         w[n] = 0.42 - 0.5·cos(2π·n/(N-1)) + 0.08·cos(4π·n/(N-1))
 *
 *         Highest sidelobe: -58.1 dB.  Sidelobe decay: 18 dB/octave.
 *         ENBW: 1.73 bins.
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_blackman(size_t N);

/**
 * @brief  Generate a Blackman-Harris 4-term window (minimum sidelobe):
 *         w[n] = 0.35875 - 0.48829·cos(2π·n/N) + 0.14128·cos(4π·n/N) - 0.01168·cos(6π·n/N)
 *
 *         Highest sidelobe: -92 dB.  ENBW: 2.00 bins.
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_blackman_harris4(size_t N);

/**
 * @brief  Generate a Nuttall 4-term window (Blackman-Nuttall):
 *         w[n] = 0.3635819 - 0.4891775·cos(2π·n/N) + 0.1365995·cos(4π·n/N) - 0.0106411·cos(6π·n/N)
 *
 *         Highest sidelobe: -98.2 dB.  ENBW: 2.02 bins.
 *         Better sidelobe than Blackman-Harris 4-term.
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_nuttall4(size_t N);

/**
 * @brief  Generate a flat-top window for precision amplitude measurement.
 *         Coefficients from the SFT3M flat-top (D'Antona, Ferrero):
 *         a0 = 0.26526, a1 = 0.50000, a2 = 0.23474
 *
 *         Passband ripple: < 0.01 dB (exceptional amplitude accuracy).
 *         Highest sidelobe: -44.7 dB.  ENBW: 3.77 bins (very wide).
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_flattop(size_t N);

/**
 * @brief  Generate a Kaiser window with adjustable shape parameter β:
 *         w[n] = I₀(β·√(1 - (2n/(N-1) - 1)²)) / I₀(β)
 *
 *         where I₀ is the modified Bessel function of the first kind.
 *
 *         β control:
 *           β = 0    → rectangular
 *           β ≈ 5.44 → similar to Hamming (sidelobe ≈ -42 dB)
 *           β ≈ 8.17 → similar to Blackman (sidelobe ≈ -60 dB)
 *           β ≈ 11.0 → sidelobe ≈ -80 dB
 *           β ≈ 14.0 → sidelobe ≈ -100 dB
 *
 *         Kaiser window is optimal: for a given sidelobe level, it minimizes
 *         mainlobe width (and vice versa), based on prolate spheroidal wave
 *         functions.
 *
 * @param  N     Window length
 * @param  beta  Shape parameter (≥ 0)
 * @return       Window config
 */
window_config_t  window_kaiser(size_t N, double beta);

/**
 * @brief  Generate a Bartlett (triangular) window:
 *         w[n] = 2n/(N-1) for 0 ≤ n ≤ (N-1)/2
 *              = 2 - 2n/(N-1) for (N-1)/2 < n ≤ N-1
 *
 *         Equivalent to convolution of two rectangular windows of length N/2.
 *         Highest sidelobe: -26.5 dB.  Sidelobe decay: 12 dB/octave.
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_bartlett(size_t N);

/**
 * @brief  Generate a Gaussian window:
 *         w[n] = exp(-½ · (α·(n - (N-1)/2) / ((N-1)/2))²)
 *
 *         The Gaussian minimizes the time-frequency uncertainty product.
 *         α controls the width: larger α → narrower window → wider mainlobe.
 *
 * @param  N     Window length
 * @param  alpha Width parameter (typically 2.5 ≤ α ≤ 3.5)
 * @return       Window config
 */
window_config_t  window_gaussian(size_t N, double alpha);

/**
 * @brief  Generate a Dolph-Chebyshev window: equiripple sidelobes at a
 *         specified level, achieving the minimum mainlobe width for a
 *         given sidelobe attenuation.
 *
 *         Constructed from Chebyshev polynomials in the frequency domain.
 *
 * @param  N               Window length
 * @param  sidelobe_db     Desired sidelobe level (dB below mainlobe, > 0)
 * @return                 Window config
 */
window_config_t  window_dolph_chebyshev(size_t N, double sidelobe_db);

/**
 * @brief  Generate a Tukey (cosine-tapered) window:
 *         Parameter α ∈ [0, 1] controls the taper fraction.
 *         - α = 0 → rectangular window
 *         - α = 1 → Hann window
 *
 *         w[n] = 1                                          for |n| ≤ (1-α)·N/2
 *              = 0.5·(1 + cos(π·(n - (1-α)·N/2) / (α·N/2)))  otherwise
 *
 * @param  N      Window length
 * @param  alpha  Taper fraction [0, 1]
 * @return        Window config
 */
window_config_t  window_tukey(size_t N, double alpha);

/**
 * @brief  Generate a Lanczos window:
 *         w[n] = sinc(2n/(N-1) - 1) = sin(π·(2n/(N-1) - 1)) / (π·(2n/(N-1) - 1))
 *
 *         Used in image resampling and interpolation.  Mainlobe: sinc-shaped.
 *
 * @param  N  Window length
 * @return    Window config
 */
window_config_t  window_lanczos(size_t N);

/* ─── L2: Window Analysis Functions ──────────────────────────────────────── */

/**
 * @brief  Compute window metrics: coherent gain, ENBW, processing loss,
 *         scalloping loss, mainlobe width, and sidelobe level.
 *
 *         These metrics are computed analytically where possible and
 *         via FFT-based oversampled frequency response for numeric metrics.
 *
 * @param  coeffs  Window coefficients
 * @param  N       Window length
 * @param  config  Output: window metrics filled in
 */
void  window_compute_metrics(double *coeffs, size_t N, window_config_t *config);

/**
 * @brief  Apply window to a signal: y[n] = w[n] · x[n]
 *
 * @param  x       Input signal (length N)
 * @param  window  Window coefficients (length N)
 * @param  N       Signal/window length
 * @param  y       Output: windowed signal (length N, can alias x)
 */
void  window_apply(const double *x, const double *window, size_t N, double *y);

/**
 * @brief  Compute window frequency response |W(ω)| at specified frequencies.
 *
 *         W(ω) = Σ_{n=0}^{N-1} w[n]·e^{-j·ω·n}
 *
 *         This is the DTFT of the window function, evaluated numerically.
 *
 * @param  window       Window coefficients (length N)
 * @param  N            Window length
 * @param  freq_points  Number of evaluation points
 * @param  W_out        Output: complex frequency response (length freq_points)
 * @param  omega_grid   Output: ω values (length freq_points, caller-allocated)
 */
void  window_frequency_response(const double *window, size_t N,
                                 size_t freq_points,
                                 fourier_complex_t *W_out, double *omega_grid);

/**
 * @brief  Estimate Kaiser β parameter from desired sidelobe attenuation.
 *
 *         Empirical formula (Kaiser & Schafer, 1980):
 *           β = 0.1102·(A - 8.7)           for A > 50 dB
 *           β = 0.5842·(A - 21)^0.4 + 0.07886·(A - 21)  for 21 ≤ A ≤ 50 dB
 *           β = 0                            for A < 21 dB
 *
 *         where A is the desired sidelobe attenuation in dB.
 *
 * @param  attenuation_db  Desired sidelobe attenuation in dB
 * @return                 Estimated β value
 */
double  window_kaiser_beta_from_attenuation(double attenuation_db);

/**
 * @brief  Estimate required window length N from desired mainlobe width.
 *
 *         Approximate relationship:  N ≈ K / BW_normalized
 *         where K is a window-dependent constant:
 *           Rectangular: K=0.89, Hann: K=1.44, Hamming: K=1.30,
 *           Blackman: K=1.68, Kaiser(β=8): K≈1.9
 *
 * @param  window_type        Window type
 * @param  mainlobe_width_hz  Desired 3dB mainlobe width (Hz)
 * @param  fs                 Sampling rate (Hz)
 * @return                    Recommended minimum window length
 */
size_t  window_recommend_length(window_type_t window_type,
                                 double mainlobe_width_hz, double fs);

/**
 * @brief  Free memory allocated for window_config_t coefficients.
 */
void  window_config_free(window_config_t *wc);

#endif /* FOURIER_WINDOW_H */
