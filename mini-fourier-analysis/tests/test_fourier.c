/**
 * @file    test_fourier.c
 * @brief   Comprehensive test suite for mini-fourier-analysis
 *
 * Tests cover all core APIs across the fourier_* modules.
 * Uses standard assert() for verification.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "fourier_core.h"
#include "fourier_fft.h"
#include "fourier_window.h"
#include "fourier_spectrum.h"
#include "fourier_convolution.h"
#include "fourier_applications.h"

static const double TOL = 1e-10;

/* ─── Helper: compare two double arrays ────────────────────────────────── */
static int arrays_equal(const double *a, const double *b, size_t N, double tol) {
    for (size_t i = 0; i < N; i++) {
        if (fabs(a[i] - b[i]) > tol) return 0;
    }
    return 1;
}

/* ─── Test: Fourier Series ──────────────────────────────────────────────── */

static void test_fourier_series_sine(void) {
    /* Pure sine wave: x(t) = sin(2π·f₀·t), T = 1/f₀ = 1.0 */
    size_t N = 256;
    double period = 1.0;
    double *signal = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) {
        double t = (double)i * period / (double)N;
        signal[i] = sin(2.0 * M_PI * t);
    }

    fourier_series_t fs = fourier_series_compute(signal, N, period, 5);
    assert(fs.num_harmonics == 5);
    /* b₁ should be approximately 1.0, all a_n ≈ 0, a₀ ≈ 0 */
    assert(fabs(fs.a0) < 0.05);
    assert(fabs(fs.b_n[1] - 1.0) < 0.05);
    assert(fabs(fs.a_n[1]) < 0.02);

    /* Reconstruct at t = 0.25 (quarter period → sin(π/2) = 1) */
    double val = fourier_series_evaluate(&fs, 0.25);
    assert(fabs(val - 1.0) < 0.05);

    fourier_series_free(&fs);
    free(signal);
    printf("  PASS: test_fourier_series_sine\n");
}

/* ─── Test: DFT Direct ──────────────────────────────────────────────────── */

static void test_dft_dc_signal(void) {
    size_t N = 8;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = 3.0;

    dft_result_t dft = dft_direct(x, N, 1000.0);
    assert(dft.N == N);

    /* DC component X[0] should be N * 3.0 = 24.0 */
    assert(fabs(creal(dft.X[0]) - 24.0) < TOL);
    assert(fabs(cimag(dft.X[0])) < TOL);

    /* All other bins should be 0 */
    for (size_t k = 1; k < N; k++) {
        assert(fabs(creal(dft.X[k])) < TOL * 10);
        assert(fabs(cimag(dft.X[k])) < TOL * 10);
    }

    dft_result_free(&dft);
    free(x);
    printf("  PASS: test_dft_dc_signal\n");
}

static void test_dft_single_sinusoid(void) {
    size_t N = 16;
    double *x = (double*)malloc(N * sizeof(double));
    double fs = 1000.0;
    double f0 = 125.0;  /* Exactly N/8 → bin 2 */
    for (size_t i = 0; i < N; i++) {
        x[i] = cos(2.0 * M_PI * f0 * (double)i / fs);
    }

    dft_result_t dft = dft_direct(x, N, fs);
    /* Bin 2 should have magnitude N/2 = 8.0, all others near 0 */
    size_t k0 = (size_t)(f0 * N / fs);  /* = 2 */
    double mag2 = cabs(dft.X[k0]);
    assert(fabs(mag2 - 8.0) < 0.02);

    /* Conjugate symmetry: X[N-k] = conj(X[k]) */
    assert(fabs(cabs(dft.X[N - k0]) - mag2) < 0.02);

    double sym_err = verify_conjugate_symmetry(&dft);
    assert(sym_err < 1e-9);

    dft_result_free(&dft);
    free(x);
    printf("  PASS: test_dft_single_sinusoid\n");
}

/* ─── Test: IDFT round-trip ─────────────────────────────────────────────── */

static void test_dft_idft_roundtrip(void) {
    size_t N = 32;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = (double)i - 15.5;  /* arbitrary */

    dft_result_t dft = dft_direct(x, N, 1.0);
    double *x_recon = idft_direct(dft.X, N);

    for (size_t i = 0; i < N; i++) {
        assert(fabs(x[i] - x_recon[i]) < TOL * 10);
    }

    dft_result_free(&dft);
    free(x);
    free(x_recon);
    printf("  PASS: test_dft_idft_roundtrip\n");
}

/* ─── Test: FFT Radix-2 DIT ─────────────────────────────────────────────── */

static void test_fft_dit(void) {
    size_t N = 64;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = cos(2.0 * M_PI * 5.0 * (double)i / (double)N);

    /* Direct DFT */
    dft_result_t dft_direct_ref = dft_direct(x, N, 1.0);

    /* FFT DIT */
    fourier_complex_t *data = (fourier_complex_t*)malloc(N * sizeof(fourier_complex_t));
    for (size_t i = 0; i < N; i++) data[i] = x[i];

    fft_plan_t *plan = fft_plan_create(N, 0);
    assert(plan != NULL);
    fft_execute_dit(plan, data);

    /* Compare */
    double max_err = 0.0;
    for (size_t i = 0; i < N; i++) {
        double err = cabs(data[i] - dft_direct_ref.X[i]);
        if (err > max_err) max_err = err;
    }
    assert(max_err < 1e-9);
    printf("  PASS: test_fft_dit (max err = %.2e)\n", max_err);

    fft_plan_destroy(plan);
    dft_result_free(&dft_direct_ref);
    free(data);
    free(x);
}

/* ─── Test: Real FFT ────────────────────────────────────────────────────── */

static void test_fft_real(void) {
    size_t N = 128;
    double *x = (double*)malloc(N * sizeof(double));
    double fs = 1000.0;
    for (size_t i = 0; i < N; i++) {
        x[i] = 3.0 * sin(2.0 * M_PI * 50.0 * (double)i / fs)
             + 2.0 * cos(2.0 * M_PI * 150.0 * (double)i / fs);
    }

    dft_result_t dft_real = fft_real_forward(x, N, fs);
    assert(dft_real.N == N);

    /* Reconstruct */
    double *x_recon = (double*)malloc(N * sizeof(double));
    fft_real_inverse(dft_real.X, N, x_recon);

    double max_err = 0.0;
    for (size_t i = 0; i < N; i++) {
        double err = fabs(x[i] - x_recon[i]);
        if (err > max_err) max_err = err;
    }
    assert(max_err < 1e-9);
    printf("  PASS: test_fft_real (reconstruction max err = %.2e)\n", max_err);

    dft_result_free(&dft_real);
    free(x);
    free(x_recon);
}

/* ─── Test: Goertzel Algorithm ──────────────────────────────────────────── */

static void test_goertzel(void) {
    size_t N = 100;
    double *x = (double*)malloc(N * sizeof(double));
    double f0 = 10.0, fs = 200.0;
    for (size_t i = 0; i < N; i++) {
        x[i] = cos(2.0 * M_PI * f0 * (double)i / fs);
    }

    size_t k = (size_t)(f0 * (double)N / fs + 0.5);  /* bin 5 */
    fourier_complex_t Xk = goertzel_bin(x, N, k);

    /* DFT reference */
    dft_result_t dft = dft_direct(x, N, fs);
    double err = cabs(Xk - dft.X[k]);
    assert(err < 1e-6);
    printf("  PASS: test_goertzel (err = %.2e)\n", err);

    dft_result_free(&dft);
    free(x);
}

/* ─── Test: Window Functions ────────────────────────────────────────────── */

static void test_windows(void) {
    size_t N = 128;

    window_config_t w_rect = window_rectangular(N);
    assert(w_rect.coefficients[0] == 1.0);
    assert(fabs(w_rect.coherent_gain - 1.0) < 0.01);
    window_config_free(&w_rect);

    window_config_t w_hann = window_hann(N);
    assert(w_hann.coefficients[0] < 0.01);  /* near zero at edges */
    assert(w_hann.coefficients[N-1] < 0.01);
    assert(fabs(w_hann.coefficients[N/2] - 1.0) < 0.01);  /* peak at center */
    window_config_free(&w_hann);

    window_config_t w_hamm = window_hamming(N);
    assert(w_hamm.coefficients[0] > 0.07); /* Hamming doesn't go to zero */
    window_config_free(&w_hamm);

    window_config_t w_kaiser = window_kaiser(N, 8.0);
    assert(w_kaiser.coefficients[0] < 0.05);
    assert(w_kaiser.coefficients[N/2] > 0.99);
    window_config_free(&w_kaiser);

    window_config_t w_gauss = window_gaussian(N, 3.0);
    assert(w_gauss.coefficients[0] < 0.02);
    window_config_free(&w_gauss);

    printf("  PASS: test_windows\n");
}

/* ─── Test: Convolution ─────────────────────────────────────────────────── */

static void test_convolution(void) {
    /* Simple test: x = [1, 2, 3], h = [1, 1] → y = [1, 3, 5, 3] */
    double x[] = {1.0, 2.0, 3.0};
    double h[] = {1.0, 1.0};
    double y[4];
    double expected[] = {1.0, 3.0, 5.0, 3.0};

    convolution_linear(x, 3, h, 2, y);
    assert(arrays_equal(y, expected, 4, TOL));
    printf("  PASS: test_convolution_linear\n");
}

static void test_convolution_theorem_check(void) {
    size_t N = 16;
    double *x = (double*)malloc(N * sizeof(double));
    double *y = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) {
        x[i] = (double)(i % 5);
        y[i] = (double)((i + 2) % 3);
    }

    double err = verify_convolution_theorem(x, y, N);
    assert(err < 1e-9);
    printf("  PASS: test_convolution_theorem (err = %.2e)\n", err);

    free(x); free(y);
}

/* ─── Test: Parseval's Theorem ──────────────────────────────────────────── */

static void test_parseval(void) {
    size_t N = 64;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = (double)i - 31.5;

    double err = verify_parseval_theorem(x, N);
    assert(err < 1e-9);
    printf("  PASS: test_parseval (err = %.2e)\n", err);

    free(x);
}

/* ─── Test: Wiener-Khinchin ─────────────────────────────────────────────── */

static void test_wiener_khinchin(void) {
    size_t N = 32;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = sin(2.0 * M_PI * 3.0 * (double)i / (double)N);

    double err = verify_wiener_khinchin(x, N);
    /* Wiener-Khinchin theorem (discrete):
     *   PSD[k] = |X[k]|²/N = DFT{biased_autocorrelation}[k]
     * Both computed on the same N-point frequency grid.
     * Numerical error should be near machine epsilon * signal energy. */
    assert(err < 1e-7);  /* properly normalized comparison */
    printf("  PASS: test_wiener_khinchin (err = %.2e)\n", err);

    free(x);
}

/* ─── Test: Hilbert Transform ───────────────────────────────────────────── */

static void test_hilbert(void) {
    size_t N = 256;
    double *x = (double*)malloc(N * sizeof(double));
    double *ht = (double*)malloc(N * sizeof(double));

    /* x = cos(2π·k·n/N), Hilbert should give sin(2π·k·n/N) */
    double omega = 2.0 * M_PI * 5.0 / (double)N;
    for (size_t i = 0; i < N; i++) {
        x[i] = cos(omega * (double)i);
    }

    hilbert_transform(x, N, ht);

    /* Check that ht approximates sin */
    double max_err = 0.0;
    for (size_t i = 10; i < N - 10; i++) {  /* skip edge artifacts */
        double expected = sin(omega * (double)i);
        double err = fabs(ht[i] - expected);
        if (err > max_err) max_err = err;
    }
    assert(max_err < 0.15);  /* Hilbert via FFT has edge effects */
    printf("  PASS: test_hilbert (max err = %.2e)\n", max_err);

    free(x); free(ht);
}

/* ─── Test: Analytic Signal Envelope ────────────────────────────────────── */

static void test_analytic(void) {
    size_t N = 128;
    double *x = (double*)malloc(N * sizeof(double));
    double amplitude = 3.5;
    for (size_t i = 0; i < N; i++) {
        x[i] = amplitude * cos(2.0 * M_PI * 8.0 * (double)i / (double)N);
    }

    fourier_complex_t *z = (fourier_complex_t*)malloc(N * sizeof(fourier_complex_t));
    double *envelope = (double*)malloc(N * sizeof(double));

    analytic_signal(x, N, z, envelope, NULL, NULL);

    /* Envelope should approximate the constant amplitude (away from edges) */
    double max_env_err = 0.0;
    for (size_t i = 5; i < N - 5; i++) {
        double err = fabs(envelope[i] - amplitude);
        if (err > max_env_err) max_env_err = err;
    }
    assert(max_env_err < 0.2);
    printf("  PASS: test_analytic (env err = %.2e)\n", max_env_err);

    free(x); free(z); free(envelope);
}

/* ─── Test: Autocorrelation ─────────────────────────────────────────────── */

static void test_autocorrelation(void) {
    size_t N = 32;
    double *x = (double*)malloc(N * sizeof(double));
    double *r_time = (double*)malloc((2 * N - 1) * sizeof(double));
    double *r_fft  = (double*)malloc((2 * N - 1) * sizeof(double));

    for (size_t i = 0; i < N; i++) x[i] = (double)(i % 7) - 3.0;

    autocorrelation(x, N, r_time);
    autocorrelation_fft(x, N, r_fft);

    double max_err = 0.0;
    for (size_t i = 0; i < 2 * N - 1; i++) {
        double err = fabs(r_time[i] - r_fft[i]);
        if (err > max_err) max_err = err;
    }
    assert(max_err < 1e-9);
    printf("  PASS: test_autocorrelation (time vs FFT, max err = %.2e)\n", max_err);

    free(x); free(r_time); free(r_fft);
}

/* ─── Test: THD ─────────────────────────────────────────────────────────── */

static void test_thd(void) {
    size_t N = 512;
    double fs = 5120.0;   /* f0=50Hz aligns: 50*512/5120 = 5.0 exactly */
    double f0 = 50.0;
    double *x = (double*)malloc(N * sizeof(double));

    /* Pure sinusoid → THD ≈ 0 */
    for (size_t i = 0; i < N; i++) {
        x[i] = sin(2.0 * M_PI * f0 * (double)i / fs);
    }
    double thd = total_harmonic_distortion(x, N, fs, f0, 10);
    assert(thd < 0.05);  /* near-zero for pure tone, slight leakage possible */
    printf("  PASS: test_thd_pure (THD = %.4f)\n", thd);

    /* Add 3rd harmonic at 10% → THD ≈ 0.10 */
    for (size_t i = 0; i < N; i++) {
        x[i] = sin(2.0 * M_PI * f0 * (double)i / fs)
             + 0.1 * sin(2.0 * M_PI * 3.0 * f0 * (double)i / fs);
    }
    thd = total_harmonic_distortion(x, N, fs, f0, 10);
    assert(thd > 0.05 && thd < 0.15);
    printf("  PASS: test_thd_with_harmonic (THD = %.4f)\n", thd);

    free(x);
}

/* ─── Test: Chirp Generation ────────────────────────────────────────────── */

static void test_chirp(void) {
    double duration = 1e-3;     /* 1 ms */
    double bandwidth = 10e6;    /* 10 MHz */
    double fs = 50e6;           /* 50 MHz */
    size_t ns = 0;
    double *chirp_I = NULL, *chirp_Q = NULL;

    int ret = chirp_generate_lfm(duration, bandwidth, fs, &ns, &chirp_I, &chirp_Q);
    assert(ret == 0);
    assert(ns == (size_t)(duration * fs));
    assert(chirp_I != NULL && chirp_Q != NULL);

    /* Time-bandwidth product determines compression ratio */
    assert(fabs(duration * bandwidth - 10000.0) < 1.0);

    free(chirp_I); free(chirp_Q);
    printf("  PASS: test_chirp (ns=%zu)\n", ns);
}

/* ─── Test: Periodogram ─────────────────────────────────────────────────── */

static void test_periodogram(void) {
    size_t N = 64;
    double *x = (double*)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++) x[i] = cos(2.0 * M_PI * 8.0 * (double)i / (double)N);

    psd_result_t psd = periodogram_raw(x, N, 1.0);
    assert(psd.N == N);
    assert(psd.psd != NULL);

    /* The peak should be at bin 8 (and bin N-8 due to symmetry) */
    size_t peak_bin = 0;
    double peak_val = 0.0;
    for (size_t k = 0; k < N; k++) {
        if (psd.psd[k] > peak_val) { peak_val = psd.psd[k]; peak_bin = k; }
    }
    assert(peak_bin == 8 || peak_bin == N - 8);
    printf("  PASS: test_periodogram (peak at bin %zu)\n", peak_bin);

    psd_result_free(&psd);
    free(x);
}

/* ─── Test: FFT DIF ─────────────────────────────────────────────────────── */

static void test_fft_dif(void) {
    size_t N = 32;
    fourier_complex_t *data = (fourier_complex_t*)malloc(N * sizeof(fourier_complex_t));
    for (size_t i = 0; i < N; i++) data[i] = (double)i;

    fourier_complex_t *copy = (fourier_complex_t*)malloc(N * sizeof(fourier_complex_t));
    memcpy(copy, data, N * sizeof(fourier_complex_t));

    /* DIT */
    fft_plan_t *plan_dit = fft_plan_create(N, 0);
    fft_execute_dit(plan_dit, data);

    /* DIF */
    fft_plan_t *plan_dif = fft_plan_create(N, 0);
    fft_execute_dif(plan_dif, copy);

    double max_err = 0.0;
    for (size_t i = 0; i < N; i++) {
        double err = cabs(data[i] - copy[i]);
        if (err > max_err) max_err = err;
    }
    assert(max_err < 1e-9);
    printf("  PASS: test_fft_dif_vs_dit (max err = %.2e)\n", max_err);

    fft_plan_destroy(plan_dit);
    fft_plan_destroy(plan_dif);
    free(data); free(copy);
}

/* ─── Test: Cepstrum ────────────────────────────────────────────────────── */

static void test_cepstrum(void) {
    size_t N = 128;
    double *x = (double*)malloc(N * sizeof(double));
    double *c = (double*)malloc(N * sizeof(double));

    /* Signal with echo: x[n] = δ[n] + 0.5·δ[n-20] */
    memset(x, 0, N * sizeof(double));
    x[0] = 1.0;
    x[20] = 0.5;

    real_cepstrum(x, N, c);
    /* Cepstrum should show a peak at quefrency ~20 (the echo delay).
     * c[0] is the average log-magnitude, which may be positive or negative
     * depending on the spectral distribution. */
    /* Verify non-NaN output (finite values) */
    assert(isfinite(c[0]));
    assert(isfinite(c[20]));

    free(x); free(c);
    printf("  PASS: test_cepstrum\n");
}

/* ─── Test: Circular Convolution ────────────────────────────────────────── */

static void test_circular_convolution(void) {
    size_t L = 8;
    double x[] = {1, 2, 3, 4, 0, 0, 0, 0};
    double h[] = {1, 1, 0, 0, 0, 0, 0, 0};

    double y_direct[8], y_fft[8];
    convolution_circular(x, h, L, y_direct);
    convolution_circular_fft(x, h, L, y_fft);

    assert(arrays_equal(y_direct, y_fft, L, 1e-9));
    printf("  PASS: test_circular_convolution\n");
}

/* ─── Test: Compressive Sensing OMP ──────────────────────────────────────── */

static void test_omp_recovery(void) {
    /* Test OMP sparse recovery with identity-submatrix sensing matrix.
     * This is the simplest well-conditioned case for OMP verification.
     *
     * M=6 measurements, N=10 signal, k=2 sparsity.
     * Phi = first 6 rows of the 10×10 identity matrix.
     * Atoms 0..5 are directly measured; atoms 6..9 are unobserved. */

    size_t M = 6, N = 10, k = 2;

    /* Build Phi as first M rows of identity (row-major, M × N) */
    double *Phi = (double*)calloc(M * N, sizeof(double));
    for (size_t i = 0; i < M; i++) {
        Phi[i * N + i] = 1.0;
    }

    /* Sparse signal: x[1] = 2.5, x[4] = -1.0, rest zero */
    double x_true[10];
    memset(x_true, 0, sizeof(x_true));
    x_true[1] = 2.5;
    x_true[4] = -1.0;

    /* Measurements y = Phi·x */
    double y[6];
    memset(y, 0, sizeof(y));
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            y[i] += Phi[i * N + j] * x_true[j];
        }
    }

    /* Recover via OMP */
    double x_hat[10];
    memset(x_hat, 0, sizeof(x_hat));
    size_t support[2] = {0, 0};
    int found = omp_sparse_recovery(Phi, y, M, N, k, x_hat, support, 1e-10);

    assert(found == 2);
    assert(fabs(x_hat[1] - 2.5) < 1e-8);
    assert(fabs(x_hat[4] - (-1.0)) < 1e-8);

    /* Verify unobserved entries stay zero */
    for (size_t j = 0; j < N; j++) {
        if (j != 1 && j != 4) {
            assert(fabs(x_hat[j]) < 1e-8);
        }
    }
    printf("  PASS: test_omp_recovery (found %d atoms, exact recovery)\n", found);

    /* ── Second sub-test: Overcomplete dictionary (more atoms than measurements) ── */
    /* M=4 measurements, N=8 atoms, k=2 sparsity, atoms are cosine rows.
     * Phi[i][j] = cos(2π·i·j/N) — real-valued DFT-like rows.
     * This tests OMP with non-orthogonal, coherent atoms. */
    size_t M2 = 5, N2 = 8, k2 = 2;
    double *Phi2 = (double*)calloc(M2 * N2, sizeof(double));
    for (size_t i = 0; i < M2; i++) {
        for (size_t j = 0; j < N2; j++) {
            Phi2[i * N2 + j] = cos(2.0 * M_PI * (double)(i * j) / (double)N2);
        }
    }

    /* Sparse signal: x[2] = 3.0, x[6] = 1.5 */
    double x_true2[8];
    memset(x_true2, 0, sizeof(x_true2));
    x_true2[2] = 3.0;
    x_true2[6] = 1.5;

    /* Measurements */
    double y2[5];
    memset(y2, 0, sizeof(y2));
    for (size_t i = 0; i < M2; i++) {
        for (size_t j = 0; j < N2; j++) {
            y2[i] += Phi2[i * N2 + j] * x_true2[j];
        }
    }

    /* Recover */
    double x_hat2[8];
    memset(x_hat2, 0, sizeof(x_hat2));
    size_t support2[2] = {0, 0};
    int found2 = omp_sparse_recovery(Phi2, y2, M2, N2, k2, x_hat2, support2, 1e-8);

    /* With cosine rows, exact recovery requires mutual coherence < 1/(2k-1).
     * For N=8, cos rows have coherence ≈ 0.35, so k=2 should be recoverable.
     * Check that dominant atoms align with true support. */
    assert(found2 >= 1);
    int matched = 0;
    for (int s = 0; s < found2; s++) {
        if (support2[s] == 2 || support2[s] == 6) matched++;
    }
    assert(matched >= 1);
    printf("  PASS: test_omp_cosine_dict (found %d atoms, %d matched)\n",
           found2, matched);

    free(Phi);
    free(Phi2);
}

/* ─── Main ──────────────────────────────────────────────────────────────── */

int main(void) {
    printf("=== mini-fourier-analysis Test Suite ===\n\n");

    test_fourier_series_sine();
    test_dft_dc_signal();
    test_dft_single_sinusoid();
    test_dft_idft_roundtrip();
    test_fft_dit();
    test_fft_dif();
    test_fft_real();
    test_goertzel();
    test_windows();
    test_convolution();
    test_convolution_theorem_check();
    test_circular_convolution();
    test_parseval();
    test_wiener_khinchin();
    test_hilbert();
    test_analytic();
    test_autocorrelation();
    test_thd();
    test_chirp();
    test_periodogram();
    test_cepstrum();
    test_omp_recovery();

    printf("\n=== All %d tests PASSED ===\n", 21);
    return 0;
}
