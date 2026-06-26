/**
 * @file    fourier_core.c
 * @brief   Core Fourier transform implementations — L1 Defs + L4 Theorems
 *
 * @details Implements the fundamental Fourier analysis operations:
 *          trigonometric and complex Fourier series, direct DFT/IDFT,
 *          amplitude/phase spectrum extraction, and discrete-time
 *          Fourier transform evaluation.
 *
 *          All functions use double precision for numerical accuracy
 *          and include boundary condition checks (null pointers, zero
 *          size, division by zero).
 *
 * Knowledge covered:
 *   L1: Fourier series, DFT, amplitude/phase spectrum
 *   L3: Complex exponential basis, orthogonality
 *   L4: Fourier inversion, Parseval checking
 *
 * Reference: Oppenheim & Willsky (1997), Ch.3-5.
 *            Oppenheim & Schafer (2010), Ch.2, 8.
 */

#include "fourier_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ─── Internal helpers ──────────────────────────────────────────────────── */

/** Safe malloc with null check */
static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "fourier_core: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── L1: Fourier Series ────────────────────────────────────────────────── */

/**
 * Compute trigonometric Fourier series coefficients via trapezoidal
 * numerical integration over one period.
 *
 * a0   = (2/T)·∫₀^T x(t) dt       → but we store a0/2 in the struct,
 *        so a0 = (1/T)·∫₀^T x(t) dt
 *        Actually following the standard form:
 *          x(t) = a0/2 + Σ [a_n cos(...) + b_n sin(...)]
 *        Then:
 *          a_n = (2/T) ∫₀^T x(t) cos(n ω₀ t) dt
 *          b_n = (2/T) ∫₀^T x(t) sin(n ω₀ t) dt
 */
fourier_series_t fourier_series_compute(const double *signal, size_t N,
                                         double period, size_t harmonics)
{
    fourier_series_t fs;
    memset(&fs, 0, sizeof(fs));

    if (!signal || N < 2 || period <= 0.0 || harmonics == 0) {
        fs.num_harmonics = 0;
        return fs;
    }

    fs.fundamental_freq = 1.0 / period;
    fs.period = period;
    fs.num_harmonics = harmonics;

    double dt = period / (double)N;
    double omega0 = 2.0 * M_PI * fs.fundamental_freq;

    /* Allocate coefficient arrays */
    fs.a_n = (double*)safe_malloc((harmonics + 1) * sizeof(double));
    fs.b_n = (double*)safe_malloc((harmonics + 1) * sizeof(double));
    /* Index 0 is unused for a_n[0], b_n[0]; a0 is separate */

    /* a0 = (2/T)·∫₀^T x(t) dt  → but standard form uses a0/2 */
    double sum = 0.0;
    for (size_t i = 0; i < N; i++) {
        /* Trapezoidal: first and last half-weighted */
        double w = (i == 0 || i == N - 1) ? 0.5 : 1.0;
        sum += w * signal[i];
    }
    fs.a0 = (2.0 / period) * sum * dt;  /* This is the standard a₀ */

    /* Compute a_n, b_n for n = 1..harmonics */
    for (size_t n = 1; n <= harmonics; n++) {
        double sum_cos = 0.0, sum_sin = 0.0;
        for (size_t i = 0; i < N; i++) {
            double t = (double)i * dt;
            double w = (i == 0 || i == N - 1) ? 0.5 : 1.0;
            double theta = (double)n * omega0 * t;
            sum_cos += w * signal[i] * cos(theta);
            sum_sin += w * signal[i] * sin(theta);
        }
        fs.a_n[n] = (2.0 / period) * sum_cos * dt;
        fs.b_n[n] = (2.0 / period) * sum_sin * dt;
    }

    return fs;
}

/**
 * Convert trigonometric Fourier series to complex exponential form.
 *
 * c₀ = a₀/2
 * c_n = (a_n - j·b_n)/2   for n > 0
 * c_{-n} = c_n* = (a_n + j·b_n)/2
 */
fourier_series_complex_t fourier_series_to_complex(const fourier_series_t *trig)
{
    fourier_series_complex_t cfs;
    memset(&cfs, 0, sizeof(cfs));

    if (!trig || trig->num_harmonics == 0) return cfs;

    cfs.fundamental_freq = trig->fundamental_freq;
    cfs.period = trig->period;
    cfs.num_harmonics = trig->num_harmonics;

    /* c_n for n = -H .. +H: 2*H + 1 elements */
    size_t H = trig->num_harmonics;
    size_t total = 2 * H + 1;
    cfs.c_n = (fourier_complex_t*)safe_malloc(total * sizeof(fourier_complex_t));

    /* c_0 at index H (center) */
    cfs.c_n[H] = trig->a0 / 2.0;

    for (size_t n = 1; n <= H; n++) {
        /* c_n = (a_n - j·b_n)/2 */
        cfs.c_n[H + n] = (trig->a_n[n] - I * trig->b_n[n]) / 2.0;
        /* c_{-n} = conjugate = (a_n + j·b_n)/2 */
        cfs.c_n[H - n] = (trig->a_n[n] + I * trig->b_n[n]) / 2.0;
    }

    return cfs;
}

/**
 * Reconstruct x(t) from Fourier series at a specific time instant.
 *
 * x(t) = a₀/2 + Σ_{n=1}^{H} [a_n·cos(n·ω₀·t) + b_n·sin(n·ω₀·t)]
 */
double fourier_series_evaluate(const fourier_series_t *series, double t)
{
    if (!series || series->num_harmonics == 0) return 0.0;

    double omega0 = 2.0 * M_PI * series->fundamental_freq;
    double result = series->a0 / 2.0;

    for (size_t n = 1; n <= series->num_harmonics; n++) {
        double theta = (double)n * omega0 * t;
        result += series->a_n[n] * cos(theta) + series->b_n[n] * sin(theta);
    }

    return result;
}

/* ─── L1: Direct DFT (O(N²)) ───────────────────────────────────────────── */

/**
 * Direct DFT computation: X[k] = Σ x[n]·exp(-j·2π·k·n / N)
 *
 * This is the pedagogical O(N²) implementation.  For production use,
 * see fourier_fft.c for the O(N·log N) Cooley-Tukey algorithm.
 *
 * Theorem: The DFT matrix W_N^{kn} = exp(-j·2π·k·n/N) is unitary
 * (up to scaling by 1/√N), so the transform is invertible.
 */
dft_result_t dft_direct(const double *x, size_t N, double fs)
{
    dft_result_t result;
    memset(&result, 0, sizeof(result));

    if (!x || N == 0) return result;

    result.N = N;
    result.sampling_rate_hz = fs;
    result.freq_resolution_hz = (N > 0 && fs > 0) ? fs / (double)N : 0.0;
    result.X = (fourier_complex_t*)safe_malloc(N * sizeof(fourier_complex_t));

    for (size_t k = 0; k < N; k++) {
        fourier_complex_t Xk = 0.0;
        double angle_base = -2.0 * M_PI * (double)k / (double)N;
        for (size_t n = 0; n < N; n++) {
            double angle = angle_base * (double)n;
            Xk += x[n] * (cos(angle) + I * sin(angle));
        }
        result.X[k] = Xk;
    }

    return result;
}

/**
 * Direct inverse DFT: x[n] = (1/N)·Σ X[k]·exp(j·2π·k·n / N)
 *
 * For real-valued input, the DFT is conjugate symmetric:
 *   X[N-k] = conj(X[k])  for k = 1..N-1
 * and the IDFT yields a purely real sequence (to machine precision).
 */
double* idft_direct(const fourier_complex_t *X, size_t N)
{
    if (!X || N == 0) return NULL;

    double *x = (double*)safe_malloc(N * sizeof(double));

    for (size_t n = 0; n < N; n++) {
        fourier_complex_t sum = 0.0;
        double angle_base = 2.0 * M_PI * (double)n / (double)N;
        for (size_t k = 0; k < N; k++) {
            double angle = angle_base * (double)k;
            sum += X[k] * (cos(angle) + I * sin(angle));
        }
        /* Take real part; for real-input DFT, imag is ~machine epsilon */
        x[n] = creal(sum) / (double)N;
    }

    return x;
}

/* ─── L1: Amplitude/Phase Spectrum ──────────────────────────────────────── */

/**
 * Extract amplitude spectrum |X[k]| and phase spectrum φ[k].
 *
 * For real-valued x[n], |X[k]| = |X[N-k]| (even symmetry)
 * and φ[k] = -φ[N-k] (odd symmetry).
 */
amplitude_phase_spectrum_t dft_to_amplitude_phase(const dft_result_t *dft)
{
    amplitude_phase_spectrum_t aps;
    memset(&aps, 0, sizeof(aps));

    if (!dft || dft->N == 0) return aps;

    size_t N = dft->N;
    aps.N = N;
    aps.amplitude    = (double*)safe_malloc(N * sizeof(double));
    aps.phase_rad    = (double*)safe_malloc(N * sizeof(double));
    aps.amplitude_db = (double*)safe_malloc(N * sizeof(double));

    for (size_t k = 0; k < N; k++) {
        double re = creal(dft->X[k]);
        double im = cimag(dft->X[k]);
        double mag = sqrt(re * re + im * im);
        aps.amplitude[k] = mag;
        aps.phase_rad[k] = atan2(im, re);
        /* 20·log₁₀(|X|), handle zero with a floor */
        if (mag > 1e-20) {
            aps.amplitude_db[k] = 20.0 * log10(mag);
        } else {
            aps.amplitude_db[k] = -400.0;  /* effectively -∞ */
        }
    }

    return aps;
}

/* ─── L4: Parseval Verification (Energy Conservation) ───────────────────── */

/**
 * Verify Parseval's theorem for a given signal.
 *
 * Parseval:  Σ_{n=0}^{N-1} |x[n]|² = (1/N)·Σ_{k=0}^{N-1} |X[k]|²
 *
 * The difference should be near machine epsilon for exact DFT.
 * Floating-point error accumulates over the summation: for N > 10^4,
 * expect difference ~ N·ε·Σ|x|².
 *
 * This function is deliberately implemented here (not in fourier_convolution.c)
 * to be available for basic DFT verification without convolution dependency.
 */
double parseval_check(const double *x, const dft_result_t *dft)
{
    if (!x || !dft || dft->N == 0) return -1.0;

    size_t N = dft->N;
    double energy_time = 0.0;
    double energy_freq = 0.0;

    for (size_t n = 0; n < N; n++) {
        energy_time += x[n] * x[n];
    }

    for (size_t k = 0; k < N; k++) {
        double re = creal(dft->X[k]);
        double im = cimag(dft->X[k]);
        energy_freq += re * re + im * im;
    }
    energy_freq /= (double)N;

    return fabs(energy_time - energy_freq);
}

/* ─── L2: Frequency Bin Mapping ─────────────────────────────────────────── */

/**
 * Map DFT bin index k to physical frequency.
 *
 * For k = 0..N/2:        f_k = k·fs/N       (positive frequencies)
 * For k = N/2+1..N-1:    f_k = (k-N)·fs/N   (negative frequencies,
 *                           but for real signals, these mirror positive)
 *
 * Convenience function used frequently in spectral display and analysis.
 */
freq_bin_t dft_bin_to_frequency(size_t k, size_t N, double fs)
{
    freq_bin_t bin;
    memset(&bin, 0, sizeof(bin));

    if (N == 0) return bin;

    bin.bin_index = k;
    bin.frequency_hz = (k <= N / 2) ? (double)k * fs / (double)N
                                     : (double)((int64_t)k - (int64_t)N) * fs / (double)N;
    bin.angular_freq = 2.0 * M_PI * bin.frequency_hz;
    bin.normalized_freq = bin.frequency_hz / fs;

    return bin;
}

/* ─── L4: Riemann-Lebesgue Lemma Verification ───────────────────────────── */

/**
 * Verify the Riemann-Lebesgue lemma numerically for a given integrable
 * function sampled at N points.
 *
 * The lemma states:  lim_{|ω|→∞} ∫ f(t)·e^{-jωt} dt = 0
 *
 * We test this by evaluating the DFT at increasing frequencies and
 * checking that |X[k]| → 0 as |k| → N/2.
 *
 * Returns the ratio of |X[N/2]| to max|X[k]| — should be small for
 * continuous functions of bounded variation.
 */
double riemann_lebesgue_ratio(const double *x, size_t N, double fs)
{
    if (!x || N < 4) return -1.0;

    dft_result_t dft = dft_direct(x, N, fs);
    if (!dft.X) return -1.0;

    /* Get amplitude spectrum */
    amplitude_phase_spectrum_t aps = dft_to_amplitude_phase(&dft);

    double max_amp = 0.0;
    for (size_t k = 0; k < N; k++) {
        if (aps.amplitude[k] > max_amp) max_amp = aps.amplitude[k];
    }

    /* Evaluate at highest frequency bin N/2 (Nyquist) */
    double high_amp = (N >= 2) ? aps.amplitude[N / 2] : 0.0;
    double ratio = (max_amp > 1e-20) ? high_amp / max_amp : 0.0;

    amplitude_phase_spectrum_free(&aps);
    dft_result_free(&dft);

    return ratio;
}

/* ─── L2: Discrete-Time Fourier Transform (DTFT) Evaluation ─────────────── */

/**
 * Evaluate DTFT of a finite sequence x[n] at specified frequency points.
 *
 * X(e^{jω}) = Σ_{n=0}^{N-1} x[n]·e^{-jωn}
 *
 * The DTFT is a continuous function of ω, periodic with period 2π.
 * This function samples it on a user-specified grid.
 */
dtft_spectrum_t dtft_evaluate(const double *x, size_t N,
                                const double *omega, size_t num_points)
{
    dtft_spectrum_t spec;
    memset(&spec, 0, sizeof(spec));

    if (!x || N == 0 || !omega || num_points == 0) return spec;

    spec.num_freq_points = num_points;
    spec.spectrum = (fourier_complex_t*)safe_malloc(num_points * sizeof(fourier_complex_t));
    spec.omega_grid = (double*)safe_malloc(num_points * sizeof(double));

    for (size_t m = 0; m < num_points; m++) {
        spec.omega_grid[m] = omega[m];
        fourier_complex_t X = 0.0;
        for (size_t n = 0; n < N; n++) {
            double theta = -omega[m] * (double)n;
            X += x[n] * (cos(theta) + I * sin(theta));
        }
        spec.spectrum[m] = X;
    }

    return spec;
}

/* ─── L3: Energy Spectral Density ───────────────────────────────────────── */

/**
 * Compute Energy Spectral Density from a continuous-time Fourier transform
 * approximation.  ESD(ω) = |X(ω)|².
 *
 * If ctft_spectrum contains the CTFT evaluated on a frequency grid,
 * this gives the ESD at those points.
 */
esd_result_t esd_from_ctft(const ctft_spectrum_t *ctft)
{
    esd_result_t esd;
    memset(&esd, 0, sizeof(esd));

    if (!ctft || ctft->num_freq_points == 0) return esd;

    esd.num_points = ctft->num_freq_points;
    esd.esd = (double*)safe_malloc(ctft->num_freq_points * sizeof(double));
    esd.freq_grid = (double*)safe_malloc(ctft->num_freq_points * sizeof(double));

    for (size_t i = 0; i < ctft->num_freq_points; i++) {
        double re = creal(ctft->spectrum[i]);
        double im = cimag(ctft->spectrum[i]);
        esd.esd[i] = re * re + im * im;
        esd.freq_grid[i] = ctft->freq_grid[i];
    }

    return esd;
}

/* ─── L2: Conjugate Symmetry Check ─────────────────────────────────────── */

/**
 * Verify conjugate symmetry of DFT for real-valued input:
 *   X[N-k] = conj(X[k])  for k = 1..N-1
 *
 * Returns the maximum absolute deviation from perfect symmetry.
 * For real-valued x[n], this should be near machine epsilon * ||x||.
 */
double verify_conjugate_symmetry(const dft_result_t *dft)
{
    if (!dft || dft->N < 2) return -1.0;

    size_t N = dft->N;
    double max_err = 0.0;

    for (size_t k = 1; k < N; k++) {
        size_t k_conj = N - k;
        fourier_complex_t diff = dft->X[k] - conj(dft->X[k_conj]);
        double err = cabs(diff);
        if (err > max_err) max_err = err;
    }

    return max_err;
}

/* ─── L3: DFT Matrix Generation ─────────────────────────────────────────── */

/**
 * Generate the N×N DFT matrix F where F[k][n] = W_N^{kn}.
 *
 * This is the matrix representation of the DFT: X = F·x.
 * Useful for teaching and verification purposes.
 *
 * The DFT matrix is:
 *   - Symmetric: F^T = F
 *   - Unitary up to scaling: F·F* = N·I
 *   - Vandermonde in the N roots of unity
 */
fourier_complex_t* dft_matrix_generate(size_t N)
{
    if (N == 0) return NULL;

    fourier_complex_t *F = (fourier_complex_t*)safe_malloc(N * N * sizeof(fourier_complex_t));

    for (size_t k = 0; k < N; k++) {
        for (size_t n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * (double)(k * n) / (double)N;
            F[k * N + n] = cos(angle) + I * sin(angle);
        }
    }

    return F;
}

/**
 * Apply DFT matrix to input: X = F·x  (matrix-vector multiply).
 * This is equivalent to dft_direct but provides the matrix interpretation.
 */
void dft_matrix_apply(const fourier_complex_t *F, size_t N,
                       const double *x, fourier_complex_t *X)
{
    if (!F || !x || !X || N == 0) return;

    for (size_t k = 0; k < N; k++) {
        fourier_complex_t sum = 0.0;
        for (size_t n = 0; n < N; n++) {
            sum += F[k * N + n] * x[n];
        }
        X[k] = sum;
    }
}

/* ─── Memory Management ─────────────────────────────────────────────────── */

void fourier_series_free(fourier_series_t *fs) {
    if (fs) {
        free(fs->a_n);
        free(fs->b_n);
        fs->a_n = NULL;
        fs->b_n = NULL;
        fs->num_harmonics = 0;
    }
}

void fourier_series_complex_free(fourier_series_complex_t *fs) {
    if (fs) {
        free(fs->c_n);
        fs->c_n = NULL;
        fs->num_harmonics = 0;
    }
}

void dft_result_free(dft_result_t *dft) {
    if (dft) {
        free(dft->X);
        dft->X = NULL;
        dft->N = 0;
    }
}

void amplitude_phase_spectrum_free(amplitude_phase_spectrum_t *aps) {
    if (aps) {
        free(aps->amplitude);
        free(aps->phase_rad);
        free(aps->amplitude_db);
        aps->amplitude = NULL;
        aps->phase_rad = NULL;
        aps->amplitude_db = NULL;
        aps->N = 0;
    }
}

void psd_result_free(psd_result_t *psd) {
    if (psd) {
        free(psd->psd);
        psd->psd = NULL;
        psd->N = 0;
    }
}

void esd_result_free(esd_result_t *esd) {
    if (esd) {
        free(esd->esd);
        free(esd->freq_grid);
        esd->esd = NULL;
        esd->freq_grid = NULL;
        esd->num_points = 0;
    }
}
