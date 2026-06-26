/**
 * @file sysid_frequency.h
 * @brief Frequency-domain system identification methods (L5 Algorithms, L8 Advanced)
 *
 * Knowledge points:
 *   L2 - Bode plot analysis for transfer function identification
 *   L2 - Frequency response from impulse response
 *   L4 - Coherence-based model validation
 *   L5 - Empirical Transfer Function Estimate (ETFE)
 *   L5 - Spectral analysis (Welch's averaged periodogram)
 *   L5 - Frequency-domain LS fitting (Levy, Sanathanan-Koerner iterative)
 *   L5 - Downsampled frequency response for efficiency
 *   L8 - Local Polynomial Method (LPM) for nonparametric FRF estimation
 *   L8 - Time-frequency analysis via Short-Time Fourier Transform
 *
 * References:
 *   Pintelon & Schoukens (2012) "System Identification: A Frequency Domain Approach"
 *   Ljung (1999) Ch. 6 (Nonparametric Time- and Frequency-Domain Methods)
 *   Schoukens, Vandersteen, Barbe, Pintelon (2009) "Nonparametric preprocessing
 *     in system identification: A powerful tool"
 *   Sanathanan & Koerner (1963) "Transfer function synthesis as a ratio of
 *     two complex polynomials"
 *   Welch, P.D. (1967) "The use of Fast Fourier Transform for estimation of
 *     power spectra"
 */

#ifndef SYSID_FREQUENCY_H
#define SYSID_FREQUENCY_H

#include "sysid_types.h"
#include "sysid_models.h"
#include <complex.h>

/*---------------------------------------------------------------------------
 * L2/L4: Coherence estimation
 *
 * Magnitude-squared coherence measures the linear relationship
 * between input u and output y at each frequency:
 *
 *   gamma^2(omega) = |Phi_yu(omega)|^2 / (Phi_uu(omega) * Phi_yy(omega))
 *
 * where Phi_yu(omega) is the cross-spectral density and Phi_uu(omega),
 * Phi_yy(omega) are the auto-spectral densities.
 *
 * Interpretation:
 *   gamma^2 ~ 1.0 : strong linear relationship
 *   gamma^2 ~ 0.0 : no linear relationship (noise or nonlinearity)
 *   gamma^2 < 0.8 : caution - nonlinearities or high noise
 *
 * This is a fundamental L4 validation tool (Bendat & Piersol, 2011).
 *---------------------------------------------------------------------------*/

/**
 * Estimate coherence between input u and output y at specified frequencies.
 *
 * Uses the DFT-based estimator computed from the ratio of cross-spectral
 * magnitude squared to the product of auto-spectra.
 *
 * @param u             Input signal [N]
 * @param y             Output signal [N]
 * @param N             Number of time samples
 * @param omega         Frequency points [Nf] (rad/sample, 0 to pi)
 * @param Nf            Number of frequency points
 * @param coherence     Output: estimated coherence [Nf] (0.0 to 1.0)
 * @param smooth_window Smoothing window length (0 = no smoothing)
 * @return              0 on success, -1 on error
 */
int sysid_coherence(const double *u, const double *y, int N,
                     const double *omega, int Nf,
                     double *coherence, int smooth_window);

/*---------------------------------------------------------------------------
 * L5: Welch's averaged periodogram for Power Spectral Density estimation
 *
 * Divides the signal into overlapping segments, applies a window function
 * (Hanning), computes the periodogram of each segment, and averages them.
 *
 * Variance reduction factor ~ 1 / (number of segments).
 *
 * Reference: Welch (1967) "The use of Fast Fourier Transform for estimation
 * of power spectra", IEEE Trans. Audio Electroacoust.
 *---------------------------------------------------------------------------*/

/**
 * Compute power spectral density (PSD) via Welch's method.
 *
 * @param x          Input signal [N]
 * @param N          Number of samples
 * @param window_len Length of each segment (e.g., N/8 for 8 averages)
 * @param overlap    Number of overlapping points between segments
 * @param freq       Output: frequency bins [Nfft] (Hz if Ts in seconds)
 * @param psd        Output: power spectral density [Nfft]
 * @param Nfft       Number of FFT bins
 * @param Ts         Sampling period [seconds]
 * @return           0 on success, -1 on error
 */
int sysid_welch_psd(const double *x, int N, int window_len, int overlap,
                     double *freq, double *psd, int Nfft, double Ts);

/*---------------------------------------------------------------------------
 * L5/L8: Sanathanan-Koerner (SK) iterative frequency-domain LS
 *
 * The Sanathanan-Koerner (1963) method fits a transfer function
 *   G(s) = N(s) / D(s)
 * to measured frequency response data.
 *
 * At each iteration t:
 *   theta_{t+1} = argmin sum_k |D_t(j*omega_k) G_meas(j*omega_k) -
 *                                N_{t+1}(j*omega_k)|^2 / |D_t(j*omega_k)|^2
 *
 * where D_t is the denominator from the previous iteration.
 * As t -> infinity, converges to true nonlinear least squares.
 *---------------------------------------------------------------------------*/

/**
 * Fit a transfer function to frequency response data using SK iteration.
 *
 * @param fd         Frequency domain data (omega, complex G values)
 * @param num_order  Numerator polynomial order (m)
 * @param den_order  Denominator polynomial order (n, den[0]=1 fixed)
 * @param num        Output: numerator coefficients [num_order+1]
 * @param den        Output: denominator coefficients [den_order+1], den[0]=1
 * @param max_iter   Maximum number of SK iterations
 * @param tol        Convergence tolerance for denominator change
 * @return           Number of iterations used, or -1 on error
 */
int sysid_sk_iteration(const sysid_freq_data_t *fd,
                        int num_order, int den_order,
                        double *num, double *den,
                        int max_iter, double tol);

/*---------------------------------------------------------------------------
 * L8: Local Polynomial Method (LPM) for nonparametric FRF estimation
 *
 * LPM (Schoukens et al., 2009) estimates the frequency response function
 * G(j*omega_k) and its derivatives using a local polynomial approximation
 * over a narrow frequency window, reducing leakage errors.
 *
 * The key idea: for a small frequency band around omega_k, approximate
 *   G(omega_k + nu) ~ G(omega_k) + G'(omega_k)*nu + ...
 * using least squares over 2W+1 neighboring frequencies.
 *---------------------------------------------------------------------------*/

/**
 * Estimate nonparametric FRF via Local Polynomial Method.
 *
 * @param u            Input signal [N]
 * @param y            Output signal [N]
 * @param N            Number of samples
 * @param omega        Target frequencies [Nf] (rad/sample)
 * @param Nf           Number of frequency points
 * @param G            Output: estimated FRF G(j*omega_f) [Nf]
 * @param poly_order   Local polynomial order (typically 2)
 * @param window_half  Half-window width in frequency bins
 * @return             0 on success, -1 on error
 */
int sysid_lpm_frf(const double *u, const double *y, int N,
                   const double *omega, int Nf,
                   double complex *G, int poly_order, int window_half);

/*---------------------------------------------------------------------------
 * L2/L6: Transfer function fitting from Bode plot data
 *
 * Classic "graphical" identification methods that provide good initial
 * estimates for iterative PEM refinement.
 *---------------------------------------------------------------------------*/

/**
 * Fit a first-order transfer function to Bode magnitude data:
 *   G(s) = K / (tau*s + 1)
 *
 * Estimation:
 *   1. K = |G(j*omega_min)| - DC gain from lowest frequency
 *   2. tau = 1 / omega_{-3dB} - time constant from cutoff
 *   3. Fallback: use phase method (-45 deg at omega = 1/tau)
 *
 * @param fd   Frequency-domain data
 * @param K    Output: estimated DC gain
 * @param tau  Output: estimated time constant [seconds]
 * @return     0 on success, -1 on error
 */
int sysid_bode_fit_1st_order(const sysid_freq_data_t *fd,
                              double *K, double *tau);

/**
 * Fit a second-order transfer function to Bode magnitude data:
 *   G(s) = K * omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
 *
 * Estimation:
 *   1. K = |G(j*omega_min)| - DC gain
 *   2. omega_n = frequency at resonance peak
 *   3. zeta from resonance peak: Mr = 1/(2*zeta*sqrt(1-zeta^2))
 *      Fallback: zeta = 1.0 (critically damped)
 *
 * @param fd       Frequency-domain data
 * @param K        Output: DC gain
 * @param omega_n  Output: natural frequency [rad/s]
 * @param zeta     Output: damping ratio (0 < zeta)
 * @return         0 on success, -1 on error
 */
int sysid_bode_fit_2nd_order(const sysid_freq_data_t *fd,
                              double *K, double *omega_n, double *zeta);

/*---------------------------------------------------------------------------
 * L5: Downsampled frequency response for efficient computation
 *
 * When the full sampling rate is much higher than needed for the frequency
 * range of interest, decimation reduces computation.
 *---------------------------------------------------------------------------*/

/**
 * Compute frequency response using decimated (downsampled) data.
 *
 * @param u          Input signal [N]
 * @param y          Output signal [N]
 * @param N          Number of original samples
 * @param freq       Frequency bins [Nf] (Hz, relative to downsampled rate)
 * @param G          Output: complex frequency response [Nf]
 * @param Nf         Number of frequency points
 * @param Ts         Original sampling period [seconds]
 * @param decimation Decimation factor (>= 1)
 * @return           0 on success, -1 on error
 */
int sysid_downsampled_frf(const double *u, const double *y, int N,
                           double *freq, double complex *G, int Nf,
                           double Ts, int decimation);

/*---------------------------------------------------------------------------
 * L8: Short-Time Fourier Transform for time-varying system analysis
 *
 * STFT reveals how the frequency content changes over time, enabling
 * detection of time-varying dynamics (parameter drift, faults).
 *
 * Definition:
 *   STFT{x}(t, f) = sum_k x(k) * w(k - t) * e^{-j*2*pi*f*k}
 *
 * where w(.) is a window function (Hanning).
 *---------------------------------------------------------------------------*/

/**
 * Compute Short-Time Fourier Transform of signal x.
 *
 * @param x          Input signal [N]
 * @param N          Number of samples
 * @param window_len Window length (time resolution)
 * @param overlap    Overlap between windows
 * @param Nfreq      Number of frequency bins per time slice
 * @param tfr        Output: time-frequency representation [n_times x Nfreq]
 * @param Ts         Sampling period [seconds]
 * @return           0 on success, -1 on error
 */
int sysid_stft(const double *x, int N, int window_len, int overlap,
                int Nfreq, double *tfr, double Ts);

/*---------------------------------------------------------------------------
 * L2: Frequency response from impulse response
 *
 * The frequency response G(e^{j*omega}) is the DTFT of g(k):
 *   G(e^{j*omega}) = sum_{k=0}^{inf} g(k) * e^{-j*omega*k}
 *
 * Fundamental relationship (Oppenheim & Schafer, 2010).
 *---------------------------------------------------------------------------*/

/**
 * Compute frequency response from truncated impulse response.
 *
 * @param impulse    Impulse response coefficients [n_impulse]
 * @param n_impulse  Number of impulse response samples
 * @param omega      Frequency points [Nf] (rad/sample)
 * @param Nf         Number of frequency points
 * @param G          Output: complex frequency response [Nf]
 * @return           0 on success, -1 on error
 */
int sysid_frf_from_impulse(const double *impulse, int n_impulse,
                            const double *omega, int Nf,
                            double complex *G);

#endif /* SYSID_FREQUENCY_H */
