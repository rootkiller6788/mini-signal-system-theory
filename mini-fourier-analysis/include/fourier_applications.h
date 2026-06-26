/**
 * @file    fourier_applications.h
 * @brief   Fourier Analysis Applications — L7 Applications
 *
 * @details Real-world application interfaces built on the Fourier analysis
 *          foundation.  These are not "fake" wrappers; each corresponds to
 *          a specific engineering use-case studied in the reference curriculum.
 *
 * Knowledge Mapping:
 *   L7 - Applications:
 *     - Audio spectrum analyzer (FFT-based real-time analysis)
 *     - Vibration analysis for rotating machinery (bearing fault detection)
 *     - EEG/MEG brain rhythm analysis (delta, theta, alpha, beta, gamma bands)
 *     - Radar pulse compression (matched filtering with chirp signals)
 *     - GPS acquisition (Doppler + code-phase search via FFT correlation)
 *     - OFDM spectrum shaping and peak-to-average power ratio (PAPR)
 *     - Audio equalizer design via frequency-domain filtering
 *     - Power quality analysis (harmonic distortion, THD)
 *     - Speech formant analysis (cepstral method)
 *
 *   L8 - Advanced Topics:
 *     - Sparse Fourier transform (SFT) for sub-Nyquist frequencies
 *     - Wavelet transform (Discrete Wavelet Transform, DWT)
 *
 * Reference: Lyons, R.G., "Understanding Digital Signal Processing" (3rd ed,
 *            2010), Ch.13 (real applications).
 */

#ifndef FOURIER_APPLICATIONS_H
#define FOURIER_APPLICATIONS_H

#include "fourier_core.h"
#include "fourier_fft.h"
#include "fourier_spectrum.h"
#include "fourier_window.h"

/* ─── L7: Audio Spectrum Analyzer ───────────────────────────────────────── */

/**
 * @brief  Audio spectrum analyzer: compute octave-band or fractional-octave
 *         spectrum from PCM audio samples.
 *
 *         Divides the spectrum into logarithmic bands (e.g., 1/3 octave)
 *         matching human auditory perception.  Each band energy is the
 *         sum of PSD over the band.
 *
 *         Standard bands: 20 Hz – 20 kHz, 1/1 octave (10 bands) or
 *         1/3 octave (31 bands) per IEC 61260.
 *
 * @param  samples      PCM audio samples
 * @param  num_samples  Number of samples
 * @param  sample_rate  Sampling rate (Hz)
 * @param  bands_per_octave  1 (octave), 3 (1/3-octave)
 * @param  band_energy  Output: energy per band
 * @param  band_centers Output: center frequency per band (Hz)
 * @param  num_bands    Output: number of bands
 * @return              0 on success
 */
int  audio_octave_spectrum(const double *samples, size_t num_samples,
                            double sample_rate, int bands_per_octave,
                            double **band_energy, double **band_centers,
                            int *num_bands);

/**
 * @brief  Compute Total Harmonic Distortion (THD) of a nominally sinusoidal signal.
 *
 *         THD = √(Σ V²_h) / V₁   where V₁ is the fundamental amplitude,
 *         V_h are harmonic amplitudes for h = 2, 3, ....
 *
 *         The fundamental is identified as the highest peak in the spectrum
 *         near the nominal frequency.
 *
 * @param  x             Input signal (presumably sinusoidal)
 * @param  N             Number of samples
 * @param  fs            Sampling rate
 * @param  nominal_freq  Nominal fundamental frequency (Hz)
 * @param  num_harmonics Number of harmonics to include (e.g., 20)
 * @return               THD value (0 = pure sine, > 0 = distortion)
 */
double  total_harmonic_distortion(const double *x, size_t N, double fs,
                                   double nominal_freq, int num_harmonics);

/* ─── L7: Vibration Analysis ────────────────────────────────────────────── */

/**
 * @brief  Detect bearing fault frequencies from vibration spectrum.
 *
 *         Rolling-element bearing faults produce characteristic frequencies:
 *           BPFO = (N_b/2)·f_r·(1 - (d/D)·cos(φ))  (outer race)
 *           BPFI = (N_b/2)·f_r·(1 + (d/D)·cos(φ))  (inner race)
 *           BSF  = (D/(2d))·f_r·(1 - (d/D)²·cos²(φ))  (ball spin)
 *           FTF  = (f_r/2)·(1 - (d/D)·cos(φ))  (cage)
 *
 *         where N_b = number of balls, f_r = shaft rotation frequency (Hz),
 *         d = ball diameter, D = pitch diameter, φ = contact angle.
 *
 *         This function searches the vibration spectrum for peaks at
 *         these characteristic frequencies and their harmonics.
 *
 * @param  spectrum     Vibration spectrum (magnitude)
 * @param  freq_axis    Frequency axis (Hz)
 * @param  num_bins     Number of FFT bins
 * @param  shaft_freq   Shaft rotation frequency (Hz)
 * @param  n_balls      Number of rolling elements
 * @param  ball_diam    Ball diameter (any units, ratio matters)
 * @param  pitch_diam   Pitch diameter
 * @param  contact_angle Contact angle (radians)
 * @param  bpfo_peak    Output: BPFO peak magnitude
 * @param  bpfi_peak    Output: BPFI peak magnitude
 * @param  severity_idx Output: overall severity index (normalized)
 */
void  bearing_fault_detect(const double *spectrum, const double *freq_axis,
                            size_t num_bins, double shaft_freq,
                            int n_balls, double ball_diam, double pitch_diam,
                            double contact_angle,
                            double *bpfo_peak, double *bpfi_peak,
                            double *severity_idx);

/* ─── L7: Radar Pulse Compression ───────────────────────────────────────── */

/**
 * @brief  Generate a linear frequency modulated (LFM / chirp) pulse.
 *
 *         s(t) = rect(t/T)·exp(j·π·(B/T)·t²)
 *
 *         where T is pulse duration, B is bandwidth.
 *         Time-bandwidth product B·T determines compression ratio.
 *
 * @param  duration    Pulse duration T (seconds)
 * @param  bandwidth   Swept bandwidth B (Hz)
 * @param  sample_rate Sampling rate (Hz)
 * @param  num_samples Output: number of samples generated
 * @param  chirp_real  Output: real part of chirp
 * @param  chirp_imag  Output: imaginary part of chirp
 * @return             0 on success
 */
int  chirp_generate_lfm(double duration, double bandwidth, double sample_rate,
                         size_t *num_samples,
                         double **chirp_real, double **chirp_imag);

/**
 * @brief  Pulse compression via matched filtering.
 *
 *         Compresses a chirp radar return by cross-correlating with the
 *         transmitted waveform (matched filter = time-reversed complex
 *         conjugate of the transmitted pulse).
 *
 *         Pulse compression ratio ≈ B·T (time-bandwidth product).
 *
 * @param  received      Received signal samples
 * @param  rx_len        Length of received signal
 * @param  tx_template   Transmitted pulse template (complex)
 * @param  tx_len        Length of template
 * @param  compressed    Output: compressed pulse (length rx_len + tx_len - 1)
 */
void  pulse_compression_matched_filter(const double *received, size_t rx_len,
                                        const fourier_complex_t *tx_template,
                                        size_t tx_len,
                                        double *compressed);

/* ─── L7: GPS Acquisition ───────────────────────────────────────────────── */

/**
 * @brief  GPS coarse/acquisition (C/A) code correlation via FFT.
 *
 *         GPS signal acquisition searches over:
 *          - PRN code phase (1023 chips, search in code dimension)
 *          - Doppler frequency (±10 kHz, search in frequency dimension)
 *
 *         The circular cross-correlation via FFT accelerates the code-phase
 *         search by a factor of ~N/log₂N compared to direct correlation.
 *
 *         R[m] = IFFT{ FFT{x} · conj(FFT{ca_code}) }
 *
 * @param  if_signal    Downconverted IF signal samples
 * @param  ca_code      C/A code sequence for a specific PRN (length 1023)
 * @param  sig_len      Length of IF signal
 * @param  code_len     C/A code length (1023)
 * @param  doppler_hz   Doppler frequency hypothesis (Hz)
 * @param  if_freq      Intermediate frequency (Hz)
 * @param  sample_rate  Sampling rate (Hz)
 * @param  correlation  Output: correlation magnitude (length code_len)
 * @param  peak_code_phase Output: detected code phase (chip offset)
 * @param  peak_ratio   Output: peak-to-second-peak ratio (acquisition metric)
 */
void  gps_acquisition_fft(const double *if_signal,
                           const int *ca_code,
                           size_t sig_len, size_t code_len,
                           double doppler_hz, double if_freq,
                           double sample_rate,
                           double *correlation,
                           size_t *peak_code_phase,
                           double *peak_ratio);

/* ─── L7: OFDM PAPR Analysis ────────────────────────────────────────────── */

/**
 * @brief  Compute Peak-to-Average Power Ratio (PAPR) of an OFDM symbol.
 *
 *         An N-subcarrier OFDM symbol is:
 *           s(t) = (1/√N)·Σ_{k=0}^{N-1} d_k·exp(j·2π·k·Δf·t)
 *
 *         PAPR = max|s(t)|² / E[|s(t)|²]
 *
 *         High PAPR is the main drawback of OFDM, requiring expensive
 *         linear power amplifiers.  Typical values for N=64: PAPR ≈ 10-12 dB.
 *
 * @param  subcarriers   Complex data symbols d_k (length N)
 * @param  N             Number of subcarriers
 * @param  oversample    Oversampling factor (≥4 to capture peaks)
 * @param  papr_db       Output: PAPR in dB
 */
void  ofdm_papr(const fourier_complex_t *subcarriers, size_t N,
                int oversample, double *papr_db);

/**
 * @brief  Compute the complementary cumulative distribution function (CCDF)
 *         of OFDM PAPR — the probability that PAPR exceeds a given threshold.
 *
 *         For N subcarriers with independent QPSK symbols:
 *           Pr(PAPR > γ) ≈ 1 - (1 - e^{-γ})^N
 *
 * @param  N           Number of subcarriers
 * @param  threshold_db PAPR threshold (dB)
 * @return             P(PAPR > threshold)
 */
double  ofdm_papr_ccdf(size_t N, double threshold_db);

/* ─── L7: Power Quality — Harmonic Analysis ─────────────────────────────── */

/**
 * @brief  Compute harmonic decomposition of a periodic power-line waveform.
 *
 *         Extracts harmonic amplitudes H_k for k = 1, 2, ..., K
 *         where H_1 is the fundamental (50/60 Hz).
 *
 *         THD (IEEE 519):  THD = √(Σ_{k=2}^{K} V_k²) / V_1
 *         Individual harmonic distortion:  IHD_k = V_k / V_1
 *
 * @param  v           Voltage samples over ≥1 cycle
 * @param  N           Number of samples
 * @param  fs          Sampling rate (Hz)
 * @param  fund_freq   Fundamental frequency (50 or 60 Hz)
 * @param  max_harm    Maximum harmonic order K
 * @param  harmonics   Output: harmonic magnitudes H_1..H_K (length K, normalized to fundamental)
 * @param  thd_percent Output: Total Harmonic Distortion (percent)
 */
void  power_harmonics(const double *v, size_t N, double fs,
                       double fund_freq, size_t max_harm,
                       double *harmonics, double *thd_percent);

/* ─── L7: Speech Formant Analysis ───────────────────────────────────────── */

/**
 * @brief  Estimate speech formant frequencies via cepstral smoothing.
 *
 *         Formants are the resonant frequencies of the vocal tract.
 *         In the cepstral domain, the vocal tract transfer function
 *         (spectral envelope) appears at low quefrencies (0-3 ms),
 *         while the excitation (pitch) appears at higher quefrencies.
 *
 *         Algorithm:
 *           1. Compute real cepstrum
 *           2. Apply low-quefrency lifter (keep first ~3 ms)
 *           3. Transform back to frequency domain
 *           4. Detect peaks → formants F1, F2, F3, F4
 *
 * @param  speech        Speech samples (windowed frame)
 * @param  N             Frame length
 * @param  fs            Sampling rate (Hz)
 * @param  num_formants  Number of formants to estimate (typically 3-4)
 * @param  formants      Output: formant frequencies (Hz), length num_formants
 * @param  bandwidths    Output: formant bandwidths (Hz), length num_formants
 */
void  speech_formants(const double *speech, size_t N, double fs,
                       size_t num_formants,
                       double *formants, double *bandwidths);

/* ─── L8: Sparse Fourier Transform (SFT) ────────────────────────────────── */

/**
 * @brief  Sparse Fourier Transform: recover k-sparse spectrum from
 *         sub-Nyquist samples.
 *
 *         When the signal in the frequency domain has only k ≪ N non-zero
 *         coefficients, the SFT can recover the spectrum from O(k·log N)
 *         time-domain samples, compared to N for full FFT.
 *
 *         This implements the Hassanieh-Indyk-Katabi-Price algorithm (2012)
 *         using frequency bin permutation and flat-window filtering.
 *
 * @param  x           Time-domain signal (length N)
 * @param  N           Signal length (power of 2)
 * @param  k           Sparsity (number of non-zero frequency components)
 * @param  X_sparse    Output: sparse spectrum (length N, only k bins non-zero)
 * @param  freq_bins   Output: indices of non-zero bins (length k)
 * @return             Number of frequency components found
 *
 * Reference: Hassanieh, H. et al., "Nearly Optimal Sparse Fourier Transform",
 *            STOC 2012.
 */
int  sparse_fourier_transform(const double *x, size_t N, size_t k,
                               fourier_complex_t *X_sparse, size_t *freq_bins);

/* ─── L8: Compressive Sensing — Orthogonal Matching Pursuit (OMP) ─────────── */

/**
 * @brief  Orthogonal Matching Pursuit (OMP) for sparse signal recovery.
 *
 *         Given measurements y = Φ·x where Φ is an M×N sensing matrix
 *         (M ≪ N), and x is k-sparse (at most k non-zero entries),
 *         OMP recovers the support and values of x greedily:
 *
 *         Algorithm:
 *           1. Initialize residual r = y, support set Λ = ∅
 *           2. For t = 1..k:
 *              a. Find column j that maximizes |⟨φ_j, r⟩|  (max correlation)
 *              b. Λ = Λ ∪ {j}
 *              c. Solve least-squares: x̂_Λ = argmin ||y - Φ_Λ·x||²
 *              d. Update residual: r = y - Φ_Λ·x̂_Λ
 *           3. Output x̂ with support Λ
 *
 *         This enables sub-Nyquist signal acquisition when the signal is
 *         sparse in the frequency domain, bridging L7 (applications) and
 *         L8 (advanced topics — compressive sensing).
 *
 *         Reference: Tropp & Gilbert (2007), "Signal Recovery from Random
 *         Measurements via Orthogonal Matching Pursuit", IEEE Trans. IT.
 *         Donoho (2006), "Compressed Sensing", IEEE Trans. IT.
 *
 * @param  Phi         Sensing matrix Φ (M × N, row-major: Phi[i*N + j])
 * @param  y           Measurement vector (length M)
 * @param  M           Number of measurements (rows of Φ)
 * @param  N           Signal dimension (columns of Φ)
 * @param  k           Target sparsity (max non-zero entries to recover)
 * @param  x_hat       Output: recovered sparse signal (length N, caller-allocated)
 * @param  support     Output: indices of non-zero entries (length k)
 * @param  tolerance   Convergence tolerance for residual norm
 * @return             Number of atoms selected (≤ k); negative on error
 *
 * Complexity: O(k·M·N + k³) time (dominated by correlations + LS solves)
 */
int  omp_sparse_recovery(const double *Phi, const double *y,
                         size_t M, size_t N, size_t k,
                         double *x_hat, size_t *support,
                         double tolerance);

/**
 * @brief  Build a partial DFT sensing matrix for compressive spectrum sensing.
 *
 *         Selects M random rows from the N×N DFT matrix to form the
 *         M×N sensing matrix Φ, where Φ[i][j] = exp(-j·2π·row[i]·j/N).
 *
 *         This models sub-Nyquist spectrum acquisition: M ≪ N time-domain
 *         samples of a k-sparse frequency spectrum.
 *
 * @param  M           Number of measurements (random DFT rows)
 * @param  N           Signal dimension (DFT size)
 * @param  row_indices Random subset of {0..N-1} of size M
 * @param  Phi         Output: sensing matrix Φ (M×N, row-major, caller-allocated)
 *                     as complex values flattened: Phi[2*(i*N + j)] = real,
 *                     Phi[2*(i*N + j)+1] = imag
 */
void  omp_dft_sensing_matrix(size_t M, size_t N, const size_t *row_indices,
                              double *Phi);

/**
 * @brief  Free memory allocated by application functions (variants).
 */
void  app_free_vector(double *v);
void  app_free_complex_vector(fourier_complex_t *v);

#endif /* FOURIER_APPLICATIONS_H */
