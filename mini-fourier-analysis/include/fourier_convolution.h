/**
 * @file    fourier_convolution.h
 * @brief   Convolution, Correlation, and Integral Transforms — L3 Math + L4 Laws
 *
 * @details Convolution and correlation are the central operations of linear
 *          time-invariant (LTI) system theory.  Their relationships with the
 *          Fourier transform form the theoretical backbone of every signal
 *          processing application.
 *
 *          This module provides:
 *            - Direct convolution (time domain, O(N·M))
 *            - Fast convolution via FFT (overlap-add, overlap-save)
 *            - Circular convolution (DFT property verification)
 *            - Cross-correlation and autocorrelation
 *            - Analytic signal and Hilbert transform
 *            - Fractional Fourier transform (FrFT)
 *
 * Knowledge Mapping:
 *   L3 - Mathematical Structures:
 *     - Linear convolution: y[n] = Σ_k x[k]·h[n-k]
 *     - Circular convolution: y[n] = Σ_k x[k]·h[(n-k) mod N]
 *     - Cross-correlation: r_xy[n] = Σ_k x[k]·y[k+n]
 *     - Autocorrelation: r_xx[n] = Σ_k x[k]·x[k+n]
 *     - Hilbert transform (90° phase shifter)
 *     - Analytic signal: z(t) = x(t) + j·Ĥ[x(t)]
 *   L4 - Theorems:
 *     - Convolution theorem:  DFT{x * h} = DFT{x} · DFT{h}
 *     - Cross-correlation theorem: DFT{r_xy} = DFT{x}·conj(DFT{y})
 *     - Wiener-Khinchin theorem: PSD = DFT{autocorrelation}
 *     - Parseval's theorem: Σ |x[n]|² = (1/N)·Σ |X[k]|²
 *     - Plancherel's theorem (L² isometry of Fourier transform)
 *     - Modulation theorem (frequency shifting)
 *     - Hilbert transform:  H(ω) = -j·sgn(ω)
 *   L5 - Algorithms:
 *     - Direct convolution (O(N·M))
 *     - FFT-based fast convolution (O(N·log N))
 *     - Hilbert transform via FFT
 *     - Fractional Fourier transform (FrFT) via eigendecomposition
 *   L6 - Problems:
 *     - Matched filtering (pulse compression, radar)
 *     - Deconvolution (Wiener filter)
 *     - Envelope detection via Hilbert transform
 *
 * Reference: Oppenheim & Willsky (1997), Ch.2-4.  Bracewell (2000),
 *            "The Fourier Transform and Its Applications", Ch.3, 6, 12.
 */

#ifndef FOURIER_CONVOLUTION_H
#define FOURIER_CONVOLUTION_H

#include "fourier_core.h"

/* ─── L3: Linear Convolution ────────────────────────────────────────────── */

/**
 * @brief  Compute linear convolution y = x * h in the time domain.
 *
 *         y[n] = Σ_{k=0}^{M-1} h[k]·x[n-k]   for n = 0..(N+M-2)
 *
 *         This is the fundamental operation of LTI filtering.
 *
 * @param  x       Input signal (length N)
 * @param  N       Length of x
 * @param  h       Impulse response (length M)
 * @param  M       Length of h
 * @param  y       Output (length N+M-1, caller-allocated)
 *
 * Complexity: O(N·M) time, O(N+M) space
 */
void  convolution_linear(const double *x, size_t N,
                         const double *h, size_t M,
                         double *y);

/**
 * @brief  Compute circular convolution y = x ⊛ h of length L.
 *
 *         y[n] = Σ_{k=0}^{L-1} x[k]·h[(n-k) mod L]
 *
 *         The DFT diagonalizes circular convolution:
 *           X[k] = DFT{x}, H[k] = DFT{h}
 *           Y[k] = X[k]·H[k]
 *           y = IDFT{Y}
 *
 *         Linear convolution can be obtained from circular convolution
 *         by zero-padding to length ≥ N+M-1.
 *
 * @param  x    First sequence (length L)
 * @param  h    Second sequence (length L)
 * @param  L    Length of both sequences
 * @param  y    Output (length L, caller-allocated)
 */
void  convolution_circular(const double *x, const double *h,
                            size_t L, double *y);

/**
 * @brief  Compute circular convolution via DFT (verify convolution theorem).
 *
 *         Uses:  y = IDFT{ DFT{x} · DFT{h} }
 *
 * @param  x    First sequence
 * @param  h    Second sequence
 * @param  L    Length (both sequences)
 * @param  y    Output
 */
void  convolution_circular_fft(const double *x, const double *h,
                                size_t L, double *y);

/* ─── L3: Correlation Functions ─────────────────────────────────────────── */

/**
 * @brief  Compute cross-correlation r_xy[τ]:
 *
 *         r_xy[τ] = Σ_{n=0}^{N-1-τ} x[n]·y[n+τ]   for τ ≥ 0
 *                   r_yx[-τ]                         for τ < 0
 *
 * @param  x       First signal (length N)
 * @param  y       Second signal (length N)
 * @param  N       Signal length
 * @param  r       Output: cross-correlation (length 2N-1, caller-allocated)
 *                  r[N-1+τ] = r_xy[τ] for τ = -(N-1)..(N-1)
 */
void  cross_correlation(const double *x, const double *y,
                        size_t N, double *r);

/**
 * @brief  Compute autocorrelation r_xx[τ] (special case of cross-correlation
 *         with x = y).
 *
 *         Properties:
 *           r_xx[0] = Σ x²[n] = energy of signal
 *           r_xx[τ] = r_xx[-τ] (even symmetry)
 *           |r_xx[τ]| ≤ r_xx[0] (maximum at zero lag)
 *
 * @param  x    Input signal (length N)
 * @param  N    Signal length
 * @param  r    Output: autocorrelation (length 2N-1)
 */
void  autocorrelation(const double *x, size_t N, double *r);

/**
 * @brief  Compute autocorrelation via FFT (Wiener-Khinchin theorem).
 *
 *         r_xx[τ] = IDFT{ |DFT{x}|² }
 *
 *         Uses zero-padding to avoid circular wraparound.
 *
 * @param  x    Input signal
 * @param  N    Signal length
 * @param  r    Output: autocorrelation (length 2N-1)
 */
void  autocorrelation_fft(const double *x, size_t N, double *r);

/**
 * @brief  Compute normalized cross-correlation (correlation coefficient):
 *
 *         ρ_xy[τ] = r_xy[τ] / √(r_xx[0]·r_yy[0])
 *
 *         ρ_xy ∈ [-1, 1], with 1 indicating perfect positive correlation.
 *
 * @param  x    First signal
 * @param  y    Second signal
 * @param  N    Signal length
 * @param  rho  Output: normalized correlation (length 2N-1)
 */
void  normalized_cross_correlation(const double *x, const double *y,
                                    size_t N, double *rho);

/* ─── L4: Theorems — Function Verification ──────────────────────────────── */

/**
 * @brief  Verify the convolution theorem:
 *         DFT{x * y} = DFT{x} · DFT{y}
 *
 *         Returns the absolute maximum difference between the left
 *         and right sides (should be near machine epsilon for exact DFT).
 *
 * @param  x    First sequence
 * @param  y    Second sequence
 * @param  N    Sequence length
 * @return      Maximum absolute error in the theorem
 */
double  verify_convolution_theorem(const double *x, const double *y, size_t N);

/**
 * @brief  Verify Parseval's theorem (energy conservation):
 *         Σ_{n=0}^{N-1} |x[n]|² = (1/N)·Σ_{k=0}^{N-1} |X[k]|²
 *
 * @param  x    Time-domain signal (length N)
 * @param  N    Signal length
 * @return      Absolute difference between left and right sides
 */
double  verify_parseval_theorem(const double *x, size_t N);

/**
 * @brief  Verify the Wiener-Khinchin theorem:
 *         PSD via periodogram == DFT of biased autocorrelation estimate
 *
 * @param  x    Signal
 * @param  N    Signal length
 * @return      Maximum absolute difference
 */
double  verify_wiener_khinchin(const double *x, size_t N);

/* ─── L6: Deconvolution ─────────────────────────────────────────────────── */

/**
 * @brief  Wiener deconvolution (optimal in MSE sense).
 *
 *         Estimates original signal s[n] from degraded measurement
 *         y[n] = s[n] * h[n] + noise[n] in the frequency domain:
 *
 *           Ŝ(f) = Y(f)·H*(f) / (|H(f)|² + NSR(f))
 *
 *         where NSR(f) = P_noise(f) / P_signal(f) is the noise-to-signal
 *         ratio.  For NSR → 0, this approaches inverse filtering.
 *
 * @param  y       Degraded measurement
 * @param  y_len   Measurement length
 * @param  h       Impulse response (point spread function / blur kernel)
 * @param  h_len   Kernel length
 * @param  nsr     Noise-to-signal ratio (constant approximation)
 * @param  s_out   Output: estimated signal (length y_len, caller-allocated)
 */
void  wiener_deconvolution(const double *y, size_t y_len,
                            const double *h, size_t h_len,
                            double nsr, double *s_out);

/* ─── L6: Hilbert Transform and Analytic Signal ─────────────────────────── */

/**
 * @brief  Compute the Hilbert transform of a real signal via FFT.
 *
 *         Frequency response:  H(ω) = -j·sgn(ω)
 *           for ω > 0: H(ω) = -j  (90° phase lag)
 *           for ω = 0: H(ω) = 0
 *           for ω < 0: H(ω) = +j  (90° phase lead)
 *
 *         Algorithm:  X = FFT{x},  Ĥ[f] = -j·sgn(f)·X[f],  Ĥ = IFFT{Ĥ}
 *
 * @param  x     Input real signal (length N)
 * @param  N     Signal length
 * @param  ht    Output: Hilbert-transformed signal (length N)
 */
void  hilbert_transform(const double *x, size_t N, double *ht);

/**
 * @brief  Construct the analytic signal:
 *         z[n] = x[n] + j·Ĥ[x[n]]
 *
 *         The analytic signal is complex-valued with no negative frequency
 *         components.  Its envelope |z[n]| and instantaneous phase ∠z[n]
 *         are well-defined.
 *
 * @param  x          Input real signal (length N)
 * @param  N          Signal length
 * @param  z          Output: analytic signal (complex, length N)
 * @param  envelope   Output: instantaneous envelope |z[n]| (length N)
 * @param  phase      Output: instantaneous phase (length N, optionally NULL)
 * @param  inst_freq  Output: instantaneous frequency dφ/dt (length N, optionally NULL)
 */
void  analytic_signal(const double *x, size_t N,
                      fourier_complex_t *z, double *envelope,
                      double *phase, double *inst_freq);

/* ─── L8: Fractional Fourier Transform (FrFT) ──────────────────────────── */

/**
 * @brief  Compute the Fractional Fourier Transform (FrFT) of order α.
 *
 *         The FrFT generalizes the ordinary Fourier transform to a
 *         continuous-order rotation in the time-frequency plane:
 *           α = 1 → ordinary Fourier transform
 *           α = 0 → identity operator
 *           α = 2 → parity operator f(t) → f(-t)
 *           α = 3 → inverse Fourier transform
 *
 *         Transform kernel:
 *           F^α{f}(u) = ∫ K_α(u, t)·f(t) dt
 *
 *         Application domains: chirp detection, time-varying filtering,
 *         quantum mechanics, optics.
 *
 * @param  x      Input signal (length N)
 * @param  N      Signal length (power of 2 recommended)
 * @param  alpha  Fractional order α (0 ≤ α ≤ 4 wraps to mod 4)
 * @param  output Output: FrFT result (complex, length N)
 *
 * Complexity: O(N·log N) via eigendecomposition of DFT matrix
 *
 * Reference: Ozaktas, Zalevsky, Kutay, "The Fractional Fourier Transform"
 *            (2001), Wiley.
 */
void  fractional_fourier_transform(const double *x, size_t N,
                                    double alpha,
                                    fourier_complex_t *output);

#endif /* FOURIER_CONVOLUTION_H */
