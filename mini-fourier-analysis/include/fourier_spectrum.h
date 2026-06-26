/**
 * @file    fourier_spectrum.h
 * @brief   Spectral Estimation and Analysis — L6 Canonical Problems
 *
 * @details Spectral estimation — the art and science of estimating the
 *          power spectral density (PSD) of a signal from finite observations.
 *          Covers classical non-parametric methods (periodogram, Welch,
 *          Bartlett, Blackman-Tukey) and parametric methods (Yule-Walker
 *          autoregressive, Burg's maximum entropy method, MUSIC).
 *
 * Knowledge Mapping:
 *   L6 - Canonical Problems:
 *     - Periodogram estimation and its statistical properties
 *     - Welch's averaged modified periodogram (bias-variance trade-off)
 *     - Bartlett's method (averaging non-overlapping periodograms)
 *     - Blackman-Tukey correlogram method (windowing autocorrelation)
 *     - Parametric AR spectral estimation (Yule-Walker equations)
 *     - Burg's maximum entropy method (MEM) for AR parameter estimation
 *     - MUSIC (MUltiple SIgnal Classification) pseudo-spectrum
 *     - Cross-spectral density and coherence estimation
 *     - Cepstrum analysis (homomorphic deconvolution)
 *   L5 - Algorithms:
 *     - Levinson-Durbin recursion for Toeplitz systems
 *     - Welch segmentation with configurable overlap
 *     - Eigen-decomposition for MUSIC algorithm
 *   L4 - Theorems:
 *     - Wiener-Khinchin theorem: PSD = Fourier transform of autocorrelation
 *     - Periodogram is asymptotically unbiased but inconsistent
 *     - Cramér spectral representation of stationary processes
 *   L2 - Concepts:
 *     - Stationarity, ergodicity in spectral estimation
 *     - Bias-variance trade-off in nonparametric methods
 *     - Resolution limit of classical estimators
 *
 * Reference: Hayes, M.H., "Statistical Digital Signal Processing and Modeling"
 *            (1996), Ch.8.  Stoica & Moses, "Spectral Analysis of Signals"
 *            (2005).  Kay, S.M., "Modern Spectral Estimation" (1988).
 */

#ifndef FOURIER_SPECTRUM_H
#define FOURIER_SPECTRUM_H

#include "fourier_core.h"
#include "fourier_window.h"

/* ─── L6: Periodogram ───────────────────────────────────────────────────── */

/**
 * @brief  Compute the raw periodogram of a signal.
 *
 *         P̂[k] = (1/N)·|X[k]|²   where X = DFT{x[n]}
 *
 *         Statistical properties for Gaussian white noise:
 *           E[P̂[k]] → P(f) as N → ∞ (asymptotically unbiased)
 *           Var[P̂[k]] ≈ P²(f) (inconsistent — variance does not decrease with N)
 *
 * @param  x    Input signal
 * @param  N    Signal length
 * @param  fs   Sampling rate (Hz)
 * @return      PSD result with periodogram estimate
 */
psd_result_t  periodogram_raw(const double *x, size_t N, double fs);

/**
 * @brief  Compute the modified periodogram (windowed periodogram):
 *
 *         P̂_w[k] = (1/N·U)·|Σ_{n=0}^{N-1} w[n]·x[n]·e^{-j·2π·k·n/N}|²
 *
 *         where U = (1/N)·Σ w²[n] is the window normalization factor
 *         to ensure an unbiased estimate.
 *
 * @param  x       Input signal
 * @param  N       Signal length
 * @param  fs      Sampling rate
 * @param  win     Window configuration
 * @return         PSD result
 */
psd_result_t  periodogram_modified(const double *x, size_t N, double fs,
                                    const window_config_t *win);

/* ─── L6: Welch's Method ────────────────────────────────────────────────── */

/**
 * @brief  Welch's averaged periodogram method for PSD estimation.
 *
 *         Divides signal into K overlapping segments of length L, applies
 *         window to each segment, computes modified periodogram of each,
 *         averages the K periodograms:
 *
 *           P̂_welch(f) = (1/K)·Σ_{i=0}^{K-1} P̂_i(f)
 *
 *         Variance reduction: Var[P̂_welch] ≈ (1/K)·P²(f)
 *         (approx. 1/K times the single-periodogram variance for 0% overlap;
 *          overlap improves variance further but segments become correlated.)
 *
 *         Typical choices: L = N/4 to N/8, overlap = 50%, Hann window.
 *
 * @param  x          Input signal
 * @param  N          Signal length
 * @param  fs         Sampling rate
 * @param  seg_len    Segment length L (FFT size, power of 2 recommended)
 * @param  overlap    Overlap fraction [0, 1), e.g., 0.5 for 50% overlap
 * @param  win_type   Window type for each segment
 * @return            PSD result with Welch estimate
 */
psd_result_t  welch_psd(const double *x, size_t N, double fs,
                         size_t seg_len, double overlap,
                         window_type_t win_type);

/* ─── L6: Bartlett's Method ─────────────────────────────────────────────── */

/**
 * @brief  Bartlett's method: average non-overlapping periodograms.
 *
 *         Special case of Welch with 0% overlap.  Simpler but lower
 *         variance reduction for the same total data length.
 *
 * @param  x        Input signal
 * @param  N        Signal length
 * @param  fs       Sampling rate
 * @param  seg_len  Segment length
 * @return          PSD result
 */
psd_result_t  bartlett_psd(const double *x, size_t N, double fs,
                            size_t seg_len);

/* ─── L6: Blackman-Tukey Correlogram Method ─────────────────────────────── */

/**
 * @brief  Blackman-Tukey spectral estimator:
 *
 *         P̂_BT(f) = Σ_{τ=-(M-1)}^{M-1} r̂[τ]·w[τ]·e^{-j·2π·f·τ}
 *
 *         where r̂[τ] is the biased sample autocorrelation estimate,
 *         w[τ] is a lag window of length 2M-1 (M ≪ N).
 *
 *         Windowing the autocorrelation reduces variance at the cost
 *         of bias (smoothing).  Window must be positive-definite.
 *
 * @param  x        Input signal
 * @param  N        Signal length
 * @param  fs       Sampling rate
 * @param  max_lag  Maximum lag M for autocorrelation
 * @param  win_type Lag window type (Bartlett recommended for BT)
 * @return          PSD result
 */
psd_result_t  blackman_tukey_psd(const double *x, size_t N, double fs,
                                  size_t max_lag, window_type_t win_type);

/* ─── L5-L6: Parametric Spectral Estimation ─────────────────────────────── */

/**
 * @brief  Autoregressive (AR) model coefficients via Yule-Walker equations.
 *
 *         Models signal as:  x[n] = -Σ_{k=1}^{p} a_k·x[n-k] + e[n]
 *         where e[n] is white noise with variance σ².
 *
 *         Yule-Walker equations (Toeplitz system):
 *           R·a = -r    where R_{ij} = r[|i-j|], r_k = E[x[n]·x[n+k]]
 *
 *         Solved via Levinson-Durbin recursion in O(p²) time.
 *
 *         PSD from AR(p) model:
 *           P_AR(f) = σ² / |1 + Σ_{k=1}^{p} a_k·e^{-j·2π·f·k}|²
 *
 * @param  x           Input signal
 * @param  N           Signal length
 * @param  order       AR model order p
 * @param  a           Output: AR coefficients a[1..p] (length p, caller-allocated)
 * @param  noise_var   Output: estimated white noise variance σ²
 *
 * Complexity: O(p² + N·p) time
 */
void  yule_walker_ar(const double *x, size_t N, size_t order,
                     double *a, double *noise_var);

/**
 * @brief  Solve Yule-Walker equations via Levinson-Durbin recursion.
 *
 *         Classic O(p²) algorithm for Toeplitz systems.  Also yields
 *         reflection coefficients (PARCOR) and prediction error powers
 *         at each intermediate order.
 *
 * @param  r            Autocorrelation sequence r[0..p] (length p+1)
 * @param  p            AR order
 * @param  a            Output: AR coefficients a[1..p]
 * @param  reflection   Output: reflection (PARCOR) coefficients (length p)
 * @param  error_power  Output: prediction error power at each order 0..p
 * @return              0 on success, 1 if singular
 */
int  levinson_durbin(const double *r, size_t p,
                     double *a, double *reflection, double *error_power);

/**
 * @brief  Burg's Maximum Entropy Method for AR spectral estimation.
 *
 *         Estimates AR parameters by minimizing forward and backward
 *         prediction errors subject to the Levinson recursion constraint.
 *         Guarantees a stable filter (poles inside unit circle).
 *
 *         Advantage over Yule-Walker: better frequency resolution for
 *         short data records and closely-spaced sinusoids.
 *
 * @param  x           Input signal
 * @param  N           Signal length
 * @param  order       AR order p
 * @param  a           Output: AR coefficients
 * @param  noise_var   Output: noise variance
 */
void  burg_mem(const double *x, size_t N, size_t order,
               double *a, double *noise_var);

/**
 * @brief  Compute AR spectrum from model coefficients.
 *
 *         P_AR(f_k) = σ² / |1 + Σ_{n=1}^{p} a_n·e^{-j·2π·k·n/N}|²
 *
 * @param  a            AR coefficients a[1..p]
 * @param  order        AR order p
 * @param  noise_var    Noise variance σ²
 * @param  num_points   Number of frequency points
 * @param  psd_out      Output PSD values (length num_points)
 * @param  freq_out     Output frequency grid (length num_points)
 */
void  ar_spectrum(const double *a, size_t order, double noise_var,
                  size_t num_points, double *psd_out, double *freq_out);

/* ─── L8: MUSIC Pseudo-Spectrum ─────────────────────────────────────────── */

/**
 * @brief  MUSIC (MUltiple SIgnal Classification) pseudo-spectrum.
 *
 *         Decomposes the signal into signal subspace (dimension d) and
 *         noise subspace.  The MUSIC pseudo-spectrum is:
 *
 *           P_MUSIC(f) = 1 / Σ_{i=d+1}^{M} |e^H(f)·v_i|²
 *
 *         where v_i are the noise eigenvectors of the autocorrelation
 *         matrix and e(f) = [1, e^{jω}, ..., e^{jω(M-1)}]^T is the
 *         frequency steering vector.
 *
 *         MUSIC resolves closely-spaced sinusoids far beyond the Fourier
 *         resolution limit (Rayleigh criterion), making it a super-resolution
 *         method.
 *
 * @param  x             Input signal
 * @param  N             Signal length
 * @param  signal_dim    Estimated number of sinusoids (signal subspace dim)
 * @param  M             Autocorrelation matrix size (M > signal_dim)
 * @param  num_points    Frequency grid points
 * @param  pseudo_psd    Output: MUSIC pseudo-spectrum
 * @param  freq_out      Output: frequency grid
 *
 * Complexity: O(M²·N + M³) time, O(M²) space
 *
 * Reference: Schmidt, R.O., "Multiple Emitter Location and Signal Parameter
 *            Estimation", IEEE Trans. AP 34(3), 1986.
 */
void  music_pseudospectrum(const double *x, size_t N,
                            size_t signal_dim, size_t M,
                            size_t num_points,
                            double *pseudo_psd, double *freq_out);

/* ─── L6: Cross-Spectral Analysis ───────────────────────────────────────── */

/**
 * @brief  Compute cross-spectral density between two signals.
 *
 *         P_xy[k] = DFT{x}·conj(DFT{y}) / N
 *
 *         The cross-spectrum is complex: magnitude indicates shared power,
 *         phase indicates phase relationship.
 *
 * @param  x     First signal
 * @param  y     Second signal
 * @param  N     Signal length
 * @param  fs    Sampling rate
 * @param  csd   Output: cross spectral density (complex, length N/2+1)
 * @param  freq  Output: frequency grid (length N/2+1)
 */
void  cross_spectral_density(const double *x, const double *y,
                              size_t N, double fs,
                              fourier_complex_t *csd, double *freq);

/**
 * @brief  Compute magnitude squared coherence (MSC):
 *
 *         C_xy(f) = |P_xy(f)|² / (P_xx(f)·P_yy(f))
 *
 *         MSC ∈ [0, 1] measures the linear relationship between x and y
 *         at each frequency.  MSC = 1 indicates perfect linear coupling.
 *
 * @param  x     First signal
 * @param  y     Second signal
 * @param  N     Signal length
 * @param  fs    Sampling rate
 * @param  msc   Output: MSC values (length N/2+1)
 * @param  freq  Output: frequency grid
 */
void  magnitude_squared_coherence(const double *x, const double *y,
                                   size_t N, double fs,
                                   double *msc, double *freq);

/* ─── L6: Cepstrum Analysis ─────────────────────────────────────────────── */

/**
 * @brief  Compute the real cepstrum of a signal:
 *
 *         c[n] = IDFT{ log|DFT{x[n]}| }
 *
 *         The cepstrum is a homomorphic transform that converts convolution
 *         in the time domain into addition in the cepstral domain:
 *           x[n] = s[n] * h[n]  →  c_x[n] = c_s[n] + c_h[n]
 *
 *         Applications: pitch detection, formant tracking, echo detection,
 *         deconvolution (separating source and filter in speech).
 *
 * @param  x       Input signal
 * @param  N       Signal length
 * @param  cepstr  Output: real cepstrum (length N, caller-allocated)
 *
 * Reference: Oppenheim & Schafer, "From Frequency to Quefrency: A History
 *            of the Cepstrum", IEEE Signal Processing Mag. 21(5), 2004.
 */
void  real_cepstrum(const double *x, size_t N, double *cepstr);

/**
 * @brief  Liftering: windowing in the cepstral domain to separate
 *         spectral envelope (low quefrency) from fine structure (high quefrency).
 *
 *         c_filtered[n] = c[n] · l[n]   where l[n] is the lifter window.
 *
 * @param  cepstrum    Input cepstrum
 * @param  N           Cepstrum length
 * @param  cutoff      Quefrency cutoff (sample index)
 * @param  pass_low    1 = keep low quefrency (envelope), 0 = keep high (pitch)
 * @param  filtered    Output: liftered cepstrum
 */
void  cepstral_liftering(const double *cepstrum, size_t N, size_t cutoff,
                          int pass_low, double *filtered);

/* ─── L2: Bandwidth Analysis ────────────────────────────────────────────── */

/**
 * @brief  Compute bandwidth metrics (3dB, ENB, occupied) from a spectrum.
 *
 *         For a given positive spectrum S[k]:
 *           - 3 dB bandwidth: frequencies where S[k] ≥ max(S)/2
 *           - ENBW: Σ S[k] · Δf / max(S)
 *           - Occupied bandwidth: narrowest band containing fraction of total power
 *
 * @param  spectrum   Power/amplitude spectrum (positive values)
 * @param  N          Number of bins
 * @param  fs         Sampling rate
 * @param  power_frac Fraction for occupied bandwidth (0.99 for 99%)
 * @return            Bandwidth metrics
 */
bandwidth_metrics_t  spectrum_bandwidth(const double *spectrum, size_t N,
                                         double fs, double power_frac);

/**
 * @brief  Free memory in psd_result_t.
 */
void  psd_result_free(psd_result_t *psd);

/**
 * @brief  Free memory in esd_result_t.
 */
void  esd_result_free(esd_result_t *esd);

#endif /* FOURIER_SPECTRUM_H */
