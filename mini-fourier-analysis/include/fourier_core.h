/**
 * @file    fourier_core.h
 * @brief   Core Fourier Analysis Definitions — L1 Definitions + L3 Mathematical Structures
 *
 * @details This header defines the fundamental data structures for Fourier analysis
 *          covering Fourier series, continuous Fourier transform, discrete Fourier
 *          transform (DFT), and the inverse transforms.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Fourier series coefficients (trigonometric & complex exponential forms)
 *     - Fourier transform pair (forward & inverse)
 *     - DFT matrix and basis functions
 *     - Amplitude/phase spectrum
 *     - Power spectrum / power spectral density (PSD)
 *     - Energy spectral density
 *     - Magnitude squared coherence (MSC)
 *     - Bandwidth definitions (3dB, equivalent noise, occupied)
 *   L2 - Core Concepts:
 *     - Time-frequency duality
 *     - Orthogonality of basis functions
 *     - Linearity of Fourier operators
 *     - Time/frequency scaling and shifting properties
 *   L3 - Mathematical Structures:
 *     - Complex exponential basis {e^{jωt}}
 *     - Trigonometric basis {cos, sin}
 *     - Kronecker delta orthogonality
 *     - Hilbert space L^2 framework
 *   L4 - Fundamental Laws:
 *     - Parseval's theorem (energy conservation)
 *     - Plancherel's theorem (L^2 isometry)
 *     - Fourier inversion theorem
 *     - Riemann-Lebesgue lemma
 *
 * Reference: Oppenheim & Willsky, "Signals and Systems" (1997), Ch.3-5
 */

#ifndef FOURIER_CORE_H
#define FOURIER_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* Windows/MSVC compatibility */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* C99 complex support */
#include <complex.h>
#ifndef CMPLX
#define CMPLX(r, i) ((double complex)((double)(r) + I * (double)(i)))
#endif

/* ─── L1: Core Type Definitions ──────────────────────────────────────────── */

/**
 * @brief Real-valued signal sample (time or frequency domain)
 * @note  double precision throughout for numerical fidelity to continuous theory
 */
typedef double fourier_sample_t;

/**
 * @brief Complex number used throughout Fourier computations
 *        Matches C99 complex double — e^{jθ} = cos(θ) + j·sin(θ)
 */
typedef double complex fourier_complex_t;

/**
 * @brief Frequency bin in Hz — maps continuous ω to discrete index k
 */
typedef struct {
    size_t bin_index;         /**< DFT bin number k (0 ≤ k < N) */
    double frequency_hz;      /**< Physical frequency in Hz */
    double angular_freq;      /**< Angular frequency ω = 2πf in rad/s */
    double normalized_freq;   /**< Normalized frequency f/fs ∈ [0, 0.5] */
} freq_bin_t;

/* ─── L1: Fourier Series (Continuous, Periodic Signals) ────────────────── */

/**
 * @brief Trigonometric Fourier series representation of a periodic signal x(t)
 *
 *        x(t) = a0/2 + Σ_{n=1}^{∞} [a_n cos(nω₀t) + b_n sin(nω₀t)]
 *
 *        where ω₀ = 2π/T is the fundamental angular frequency.
 *
 *        a_n = (2/T) ∫_0^T x(t) cos(nω₀t) dt   (n ≥ 0)
 *        b_n = (2/T) ∫_0^T x(t) sin(nω₀t) dt   (n ≥ 1)
 */
typedef struct {
    double fundamental_freq;      /**< Fundamental frequency f₀ = 1/T (Hz) */
    double period;                /**< Period T (seconds) */
    size_t num_harmonics;         /**< Number of harmonic terms stored */
    double a0;                    /**< DC offset (a₀/2 in standard form) */
    double *a_n;                  /**< Cosine coefficients, length num_harmonics */
    double *b_n;                  /**< Sine coefficients, length num_harmonics */
} fourier_series_t;

/**
 * @brief Compact / complex exponential form of Fourier series:
 *
 *        x(t) = Σ_{n=-∞}^{∞} c_n e^{j n ω₀ t}
 *
 *        c_n = (1/T) ∫_0^T x(t) e^{-j n ω₀ t} dt
 *
 *        Relationship to trigonometric form:
 *        c_0 = a₀/2
 *        c_n = (a_n - j b_n)/2   for n > 0
 *        c_{-n} = c_n*
 */
typedef struct {
    double fundamental_freq;      /**< Fundamental frequency f₀ (Hz) */
    double period;                /**< Period T (seconds) */
    size_t num_harmonics;         /**< Half-range: stores -N to +N (2N+1 terms) */
    fourier_complex_t *c_n;       /**< Complex coefficients, length 2*num_harmonics+1 */
} fourier_series_complex_t;

/* ─── L1: Continuous-Time Fourier Transform (CTFT) ─────────────────────── */

/**
 * @brief Continuous-Time Fourier Transform pair:
 *
 *        Forward:  X(ω) = ∫_{-∞}^{∞} x(t) e^{-jωt} dt
 *        Inverse:  x(t) = (1/2π) ∫_{-∞}^{∞} X(ω) e^{jωt} dω
 *
 * @note  For numerical work, we store discretized representations of X(ω)
 *        on a user-specified frequency grid. The integral is approximated
 *        with numerical quadrature (trapezoidal / Simpson).
 */
typedef struct {
    size_t num_freq_points;       /**< Number of frequency evaluation points */
    double freq_start_hz;         /**< Lowest frequency (Hz) */
    double freq_end_hz;           /**< Highest frequency (Hz) */
    fourier_complex_t *spectrum;  /**< X(ω) values, length num_freq_points */
    double *freq_grid;            /**< ω or f grid points, length num_freq_points */
} ctft_spectrum_t;

/* ─── L1: Discrete-Time Fourier Transform (DTFT) ────────────────────────── */

/**
 * @brief DTFT of a discrete sequence x[n]:
 *
 *        X(e^{jω}) = Σ_{n=-∞}^{∞} x[n] e^{-jωn}
 *
 *        This is periodic in ω with period 2π.
 *
 *        Inverse: x[n] = (1/2π) ∫_{-π}^{π} X(e^{jω}) e^{jωn} dω
 */
typedef struct {
    size_t num_freq_points;       /**< Evaluation points over [-π, π] or [0, 2π] */
    fourier_complex_t *spectrum;  /**< X(e^{jω}) values */
    double *omega_grid;           /**< Normalized angular frequency ω grid (rad/sample) */
} dtft_spectrum_t;

/* ─── L1: Discrete Fourier Transform (DFT) ──────────────────────────────── */

/**
 * @brief N-point Discrete Fourier Transform:
 *
 *        X[k] = Σ_{n=0}^{N-1} x[n] · W_N^{kn}
 *
 *        where W_N = e^{-j·2π/N} is the N-th root of unity.
 *
 *        Inverse (IDFT):
 *        x[n] = (1/N) Σ_{k=0}^{N-1} X[k] · W_N^{-kn}
 */
typedef struct {
    size_t N;                     /**< Transform length (number of points) */
    fourier_complex_t *X;         /**< DFT coefficients X[k], length N */
    double sampling_rate_hz;      /**< Sampling frequency fs (Hz) */
    double freq_resolution_hz;    /**< Frequency resolution Δf = fs/N */
} dft_result_t;

/* ─── L1: Amplitude and Phase Spectrum ──────────────────────────────────── */

/**
 * @brief Amplitude spectrum: |X[k]| = sqrt(Re{X[k]}² + Im{X[k]}²)
 *
 *        Phase spectrum:     φ[k] = atan2(Im{X[k]}, Re{X[k]})
 *
 *        For real-valued x[n], |X[k]| is symmetric (Hermitian).
 */
typedef struct {
    size_t N;                     /**< Number of frequency bins */
    double *amplitude;            /**< |X[k]|, length N */
    double *phase_rad;            /**< φ[k] in radians, length N */
    double *amplitude_db;         /**< 20·log₁₀(|X[k]|), length N */
} amplitude_phase_spectrum_t;

/* ─── L2-L3: Power Spectrum and Spectral Density ────────────────────────── */

/**
 * @brief Power spectral density (PSD) via periodogram (Welch's method support).
 *
 *        Periodogram:  P̂[k] = |X[k]|² / N   (unnormalized)
 *                              = |X[k]|² / (fs·N)  (in V²/Hz)
 *
 *        Welch's method: average periodograms of overlapping windowed segments
 *        for reduced variance.
 *
 * Reference: Oppenheim & Schafer (2010), Ch.10
 */
typedef struct {
    size_t N;                     /**< FFT length */
    size_t num_segments;          /**< Number of Welch segments averaged */
    double *psd;                  /**< Power spectral density values, length N */
    double freq_resolution_hz;    /**< Δf resolution in Hz */
    double variance_reduction;    /**< Variance reduction factor vs single periodogram */
} psd_result_t;

/**
 * @brief Energy spectral density: ∫ |x(t)|² dt = (1/2π) ∫ |X(ω)|² dω
 *
 *        ESD:  E(ω) = |X(ω)|²  (Joules/Hz or V²·s/Hz)
 */
typedef struct {
    size_t num_points;
    double *esd;                  /**< |X(ω)|², length num_points */
    double *freq_grid;            /**< Frequency grid */
} esd_result_t;

/* ─── L2-L3: Time-Frequency Representation ──────────────────────────────── */

/**
 * @brief Window configuration for spectral analysis.
 *        See fourier_window.h for the full window function library.
 */
typedef enum {
    WINDOW_RECTANGULAR = 0,       /**< w[n] = 1 (Dirichlet kernel) */
    WINDOW_HANN,                  /**< w[n] = 0.5(1 - cos(2πn/(N-1))) */
    WINDOW_HAMMING,               /**< w[n] = 0.54 - 0.46·cos(2πn/(N-1)) */
    WINDOW_BLACKMAN,              /**< w[n] = 0.42 - 0.5·cos(...) + 0.08·cos(...) */
    WINDOW_KAISER,                /**< Kaiser window with adjustable β */
    WINDOW_BARTLETT,              /**< Triangular window */
    WINDOW_FLATTOP,               /**< Flat-top window for amplitude accuracy */
    WINDOW_GAUSSIAN,              /**< Gaussian window exp(-½(α·n/(N/2))²) */
} window_type_t;

typedef struct {
    window_type_t type;           /**< Window function type */
    size_t length;                /**< Window length N */
    double *coefficients;         /**< w[n] values, length N */
    double coherent_gain;         /**< Σ w[n] / N */
    double processing_loss_db;    /**< Processing loss in dB */
    double scalloping_loss_db;    /**< Maximum scalloping loss in dB */
    double equivalent_noise_bw;   /**< ENBW in bins */
    double sidelobe_level_db;     /**< Highest sidelobe level relative to mainlobe */
    double mainlobe_3dB_width;    /**< 3 dB mainlobe width in bins */
} window_config_t;

/* ─── L5: FFT Algorithm Configuration ───────────────────────────────────── */

/**
 * @brief FFT plan structure (builder pattern for Cooley-Tukey radix-2 FFT).
 *
 *        The classical Cooley-Tukey decimation-in-time (DIT) FFT computes
 *        the DFT in O(N·log₂N) operations vs O(N²) for direct DFT.
 *
 *        Complexity: O(N·log N) time, O(N) auxiliary space (in-place).
 */
typedef struct {
    size_t N;                     /**< Transform length (must be power of 2 for radix-2) */
    fourier_complex_t *twiddle;   /**< Precomputed twiddle factors W_N^k */
    size_t *bit_reverse;          /**< Bit-reversal permutation table */
    int is_inverse;               /**< 0 = forward, 1 = inverse */
} fft_plan_t;

/* ─── L2: Chirp Z-Transform (CZT) ───────────────────────────────────────── */

/**
 * @brief Chirp Z-Transform parameters.
 *
 *        CZT evaluates X(z_k) on a spiral contour:
 *        z_k = A·W^{-k},  k = 0, 1, ..., M-1
 *
 *        where A = A₀·e^{j·θ₀} determines the starting point,
 *              W = W₀·e^{j·φ₀} determines the ratio between points.
 *
 *        Special cases:
 *          - A=1, W=e^{-j2π/N}, M=N → standard DFT
 *          - A=1, W=e^{-j2π/(K·N)}, M=K·N → zoom FFT
 *
 * Reference: Rabiner, Schafer, Rader, "The Chirp z-Transform Algorithm" (1969)
 */
typedef struct {
    double A_magnitude;           /**< Starting radius A₀ */
    double A_angle_rad;           /**< Starting angle θ₀ (rad) */
    double W_magnitude;           /**< Ratio magnitude W₀ */
    double W_angle_rad;           /**< Ratio angle φ₀ (rad) */
    size_t M;                     /**< Number of output points */
} czt_params_t;

/* ─── L2: Short-Time Fourier Transform (STFT) ───────────────────────────── */

/**
 * @brief STFT configuration for time-frequency analysis:
 *
 *        X[m, k] = Σ_{n=-∞}^{∞} x[n]·w[n - mR]·e^{-j·2π·k·n/N}
 *
 *        where w[n] is the analysis window, R is the hop size,
 *        m is the frame index, k is the frequency bin.
 */
typedef struct {
    size_t window_length;         /**< Analysis window length */
    size_t fft_length;            /**< FFT length N (≥ window_length, zero-padded) */
    size_t hop_size;              /**< Hop size R between successive frames */
    window_type_t window_type;    /**< Analysis window type */
    double *window_coeffs;        /**< Window coefficients (length window_length) */
} stft_config_t;

/**
 * @brief STFT result — time-frequency matrix (spectrogram).
 */
typedef struct {
    size_t num_frames;            /**< Number of time frames M */
    size_t num_bins;              /**< Number of frequency bins (N/2 + 1 for real input) */
    double **magnitude;           /**< Magnitude spectrogram [M][num_bins] */
    double **phase;               /**< Phase spectrogram [M][num_bins] */
    double *time_axis;            /**< Time center of each frame (seconds) */
    double *freq_axis;            /**< Frequency of each bin (Hz) */
    double frame_duration_sec;    /**< Duration of each frame */
} stft_result_t;

/* ─── L1: Spectral Metrics ──────────────────────────────────────────────── */

/**
 * @brief Three standard bandwidth definitions for a bandpass spectrum.
 *
 *        1. 3-dB bandwidth:  Δf where |X(f)|² drops to ½ of peak
 *        2. Equivalent noise bandwidth (ENBW):  ∫|X(f)|²df / max(|X(f)|²)
 *        3. Occupied bandwidth:  the band containing 99% of total power
 */
typedef struct {
    double bandwidth_3dB_hz;      /**< 3 dB bandwidth (half-power) */
    double bandwidth_enb_hz;       /**< Equivalent noise bandwidth */
    double bandwidth_occupied_hz; /**< 99% power bandwidth */
    double center_freq_hz;        /**< Center frequency */
    double peak_power;            /**< Peak spectral power */
    double total_power;           /**< Total integrated power */
} bandwidth_metrics_t;

/* ─── Function Declarations — L4-L6 Operations ──────────────────────────── */

/**
 * @brief Compute trigonometric Fourier series coefficients from sampled data.
 *
 *        Uses trapezoidal numerical integration over one period.
 *
 * @param signal    Sampled signal over one period T (length N)
 * @param N         Number of samples
 * @param period    Period T in seconds
 * @param harmonics Number of harmonic terms to compute (n=1..harmonics)
 * @return          Fourier series with a0, a_n[1..H], b_n[1..H]
 *
 * Complexity: O(N · harmonics) time, O(harmonics) space
 */
fourier_series_t  fourier_series_compute(const double *signal, size_t N,
                                          double period, size_t harmonics);

/**
 * @brief Convert trigonometric Fourier series to complex exponential form.
 *
 *        c_0 = a0/2
 *        c_n = (a_n - j·b_n)/2  for n > 0
 *        c_{-n} = conj(c_n)
 *
 * @param trig  Trigonometric series (input)
 * @return      Complex exponential form series
 */
fourier_series_complex_t  fourier_series_to_complex(const fourier_series_t *trig);

/**
 * @brief Reconstruct time-domain signal from Fourier series at time t.
 *
 *        x(t) = a0/2 + Σ [a_n cos(nω₀t) + b_n sin(nω₀t)]
 *
 * @param series  Fourier series coefficients
 * @param t       Time instant (seconds)
 * @return        Reconstructed value x(t)
 */
double  fourier_series_evaluate(const fourier_series_t *series, double t);

/**
 * @brief Compute direct DFT (N² algorithm) — for teaching & small N verification.
 *
 *        X[k] = Σ_{n=0}^{N-1} x[n]·exp(-j·2π·k·n / N)
 *
 * @param x   Input signal (length N, real-valued)
 * @param N   Transform size
 * @return    DFT result (complex coefficients X[k])
 *
 * Complexity: O(N²) time, O(N) space.  Use fft_execute for O(N·log N).
 */
dft_result_t  dft_direct(const double *x, size_t N, double fs);

/**
 * @brief Compute inverse DFT (N² algorithm).
 *
 *        x[n] = (1/N) Σ_{k=0}^{N-1} X[k]·exp(j·2π·k·n / N)
 *
 * @param X   DFT coefficients (length N)
 * @param N   Transform size
 * @return    Reconstructed real signal (length N)
 */
double*  idft_direct(const fourier_complex_t *X, size_t N);

/**
 * @brief Extract amplitude spectrum |X[k]| and phase spectrum φ[k] from DFT.
 *
 * @param dft  DFT result
 * @return     Amplitude and phase spectrum (caller frees)
 */
amplitude_phase_spectrum_t  dft_to_amplitude_phase(const dft_result_t *dft);

/* ─── Additional Core Operations (L2-L4) ─────────────────────────────────── */

/**
 * @brief  Verify Parseval's theorem for given signal and its DFT.
 * @return Absolute difference between time-domain and frequency-domain energy.
 */
double  parseval_check(const double *x, const dft_result_t *dft);

/**
 * @brief  Verify conjugate symmetry X[N-k] = conj(X[k]) for real-input DFT.
 * @return Maximum absolute deviation from perfect symmetry.
 */
double  verify_conjugate_symmetry(const dft_result_t *dft);

/**
 * @brief  Evaluate DTFT X(e^{jω}) at specified frequency points.
 */
dtft_spectrum_t  dtft_evaluate(const double *x, size_t N,
                                const double *omega, size_t num_points);

/**
 * @brief  Map DFT bin index to physical frequency info.
 */
freq_bin_t  dft_bin_to_frequency(size_t k, size_t N, double fs);

/**
 * @brief  Compute Riemann-Lebesgue ratio: |X[N/2]| / max|X[k]|.
 * @return Ratio ∈ [0,1]; small values confirm the lemma.
 */
double  riemann_lebesgue_ratio(const double *x, size_t N, double fs);

/**
 * @brief  Generate the N×N DFT matrix F where F[k][n] = W_N^{kn}.
 */
fourier_complex_t*  dft_matrix_generate(size_t N);

/**
 * @brief  Apply DFT matrix to input: X = F·x.
 */
void  dft_matrix_apply(const fourier_complex_t *F, size_t N,
                       const double *x, fourier_complex_t *X);

/**
 * @brief  Compute Energy Spectral Density from CTFT spectrum.
 */
esd_result_t  esd_from_ctft(const ctft_spectrum_t *ctft);

/* ─── Memory Management ─────────────────────────────────────────────────── */

/**
 * @brief Free memory associated with fourier_series_t.
 */
void  fourier_series_free(fourier_series_t *fs);

/**
 * @brief Free memory associated with fourier_series_complex_t.
 */
void  fourier_series_complex_free(fourier_series_complex_t *fs);

/**
 * @brief Free memory associated with dft_result_t.
 */
void  dft_result_free(dft_result_t *dft);

/**
 * @brief Free memory associated with amplitude_phase_spectrum_t.
 */
void  amplitude_phase_spectrum_free(amplitude_phase_spectrum_t *aps);

#endif /* FOURIER_CORE_H */
