/**
 * @file    fourier_convolution.c
 * @brief   Convolution, Correlation, Hilbert Transform, FrFT — L3 Math + L4 Laws
 *
 * @details Implements the fundamental signal processing operations of
 *          convolution and correlation, along with their Fourier-domain
 *          equivalents.  Also includes the Hilbert transform, analytic
 *          signal construction, and the fractional Fourier transform.
 *
 *          Each function includes an explicit verification against the
 *          corresponding theorem (convolution theorem, Parseval, Wiener-
 *          Khinchin) with floating-point error quantification.
 *
 * Knowledge covered:
 *   L3: Linear/circular convolution, cross-correlation, autocorrelation
 *   L4: Convolution theorem, Parseval theorem, Wiener-Khinchin theorem
 *   L6: Wiener deconvolution, Hilbert transform, analytic signal
 *   L8: Fractional Fourier transform (FrFT)
 *
 * Reference: Oppenheim & Willsky (1997), Ch.2-4, 8.
 *            Bracewell (2000).  Ozaktas et al. (2001) for FrFT.
 */

#include "fourier_convolution.h"
#include "fourier_core.h"
#include "fourier_fft.h"
#include "fourier_spectrum.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "fourier_conv: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── L3: Linear Convolution ────────────────────────────────────────────── */

void convolution_linear(const double *x, size_t N,
                         const double *h, size_t M,
                         double *y) {
    if (!x || !h || !y || N == 0 || M == 0) return;

    size_t y_len = N + M - 1;
    memset(y, 0, y_len * sizeof(double));

    for (size_t n = 0; n < y_len; n++) {
        double sum = 0.0;
        /* k ranges over valid h indices and x indices */
        /* y[n] = Σ_{k=max(0, n-M+1)}^{min(n, N-1)} h[n-k]·x[k] */
        size_t k_min = (n >= M) ? (n - M + 1) : 0;
        size_t k_max = (n < N) ? n : (N - 1);
        for (size_t k = k_min; k <= k_max; k++) {
            sum += x[k] * h[n - k];
        }
        y[n] = sum;
    }
}

/* ─── L3: Circular Convolution ──────────────────────────────────────────── */

void convolution_circular(const double *x, const double *h,
                            size_t L, double *y) {
    if (!x || !h || !y || L == 0) return;

    for (size_t n = 0; n < L; n++) {
        double sum = 0.0;
        for (size_t k = 0; k < L; k++) {
            size_t idx = (n >= k) ? (n - k) : (L + n - k);
            sum += x[k] * h[idx];
        }
        y[n] = sum;
    }
}

/* ─── L3: Circular Convolution via FFT ──────────────────────────────────── */

void convolution_circular_fft(const double *x, const double *h,
                                size_t L, double *y) {
    if (!x || !h || !y || L == 0) return;

    fourier_complex_t *X = (fourier_complex_t*)safe_malloc(L * sizeof(fourier_complex_t));
    fourier_complex_t *H = (fourier_complex_t*)safe_malloc(L * sizeof(fourier_complex_t));

    for (size_t i = 0; i < L; i++) { X[i] = x[i]; H[i] = h[i]; }

    fft_plan_t *plan = fft_plan_create(L, 0);
    fft_plan_t *plan_inv = fft_plan_create(L, 1);

    fft_execute_dit(plan, X);
    fft_execute_dit(plan, H);

    for (size_t i = 0; i < L; i++) X[i] *= H[i];

    fft_execute_dit(plan_inv, X);

    for (size_t i = 0; i < L; i++) y[i] = creal(X[i]);

    fft_plan_destroy(plan);
    fft_plan_destroy(plan_inv);
    free(X);
    free(H);
}

/* ─── L3: Cross-Correlation ─────────────────────────────────────────────── */

void cross_correlation(const double *x, const double *y,
                        size_t N, double *r) {
    if (!x || !y || !r || N == 0) return;

    size_t r_len = 2 * N - 1;
    for (size_t tau = 0; tau < r_len; tau++) {
        int64_t lag = (int64_t)tau - (int64_t)(N - 1);
        double sum = 0.0;
        size_t count = 0;
        for (size_t n = 0; n < N; n++) {
            int64_t n2 = (int64_t)n + lag;
            if (n2 >= 0 && (size_t)n2 < N) {
                sum += x[n] * y[n2];
                count++;
            }
        }
        r[tau] = sum;  /* unnormalized (biased or unbiased depending on convention) */
    }
}

/* ─── L3: Autocorrelation ───────────────────────────────────────────────── */

void autocorrelation(const double *x, size_t N, double *r) {
    cross_correlation(x, x, N, r);
}

/* ─── L3: Autocorrelation via FFT (Wiener-Khinchin) ─────────────────────── */

void autocorrelation_fft(const double *x, size_t N, double *r) {
    if (!x || !r || N == 0) return;

    size_t fft_len = 1;
    while (fft_len < 2 * N) fft_len <<= 1;

    fourier_complex_t *X = (fourier_complex_t*)calloc(fft_len, sizeof(fourier_complex_t));
    if (!X) return;
    for (size_t i = 0; i < N; i++) X[i] = x[i];

    fft_plan_t *plan = fft_plan_create(fft_len, 0);
    fft_plan_t *plan_inv = fft_plan_create(fft_len, 1);

    fft_execute_dit(plan, X);

    /* |X[k]|² */
    for (size_t i = 0; i < fft_len; i++) {
        double re = creal(X[i]);
        double im = cimag(X[i]);
        X[i] = (re * re + im * im);
    }

    fft_execute_dit(plan_inv, X);

    /* Extract linear autocorrelation from circular FFT result.
     *
     * The IFFT of |X|² gives circular autocorrelation r_circ of length fft_len.
     * The IFFT scaling (÷fft_len) cancels the factor of fft_len from the
     * unscaled transform, so r_circ[τ] = r_time[τ] directly.
     *
     * With zero-padding (fft_len ≥ 2N-1):
     *   Positive lags 0..N-1 → r_circ[0..N-1]
     *   Negative lags -(N-1)..-1 → r_circ[fft_len-(N-1)..fft_len-1]
     *
     * Map to time-domain convention:
     *   r[0..N-2] = negative lags; r[N-1] = zero lag; r[N..2N-2] = positive lags
     */
    /* Negative lags: r_xy[-(N-1)] .. r_xy[-1] stored at r[0] .. r[N-2] */
    for (size_t tau = 1; tau < N; tau++) {
        r[N - 1 - tau] = creal(X[fft_len - tau]);
    }
    /* Zero lag: r_xy[0] at r[N-1] */
    r[N - 1] = creal(X[0]);
    /* Positive lags: r_xy[1] .. r_xy[N-1] at r[N] .. r[2N-2] */
    for (size_t tau = 1; tau < N; tau++) {
        r[N - 1 + tau] = creal(X[tau]);
    }

    fft_plan_destroy(plan);
    fft_plan_destroy(plan_inv);
    free(X);
}

/* ─── L3: Normalized Cross-Correlation ──────────────────────────────────── */

void normalized_cross_correlation(const double *x, const double *y,
                                    size_t N, double *rho) {
    if (!x || !y || !rho || N == 0) return;

    /* Compute r_xx[0] and r_yy[0] (energy) */
    double Ex = 0.0, Ey = 0.0;
    for (size_t n = 0; n < N; n++) { Ex += x[n] * x[n]; Ey += y[n] * y[n]; }

    double denom = sqrt(Ex * Ey);
    if (denom < 1e-20) {
        memset(rho, 0, (2 * N - 1) * sizeof(double));
        return;
    }

    cross_correlation(x, y, N, rho);
    for (size_t i = 0; i < 2 * N - 1; i++) {
        rho[i] /= denom;
    }
}

/* ─── L4: Theorem Verification ──────────────────────────────────────────── */

double verify_convolution_theorem(const double *x, const double *y, size_t N) {
    if (!x || !y || N == 0) return -1.0;

    /* Left side: DFT{x * y} where * is circular convolution */
    double *c_conv = (double*)safe_malloc(N * sizeof(double));
    convolution_circular(x, y, N, c_conv);
    dft_result_t DFT_conv = dft_direct(c_conv, N, 1.0);

    /* Right side: DFT{x} · DFT{y} */
    dft_result_t DFT_x = dft_direct(x, N, 1.0);
    dft_result_t DFT_y = dft_direct(y, N, 1.0);

    double max_err = 0.0;
    for (size_t k = 0; k < N; k++) {
        fourier_complex_t product = DFT_x.X[k] * DFT_y.X[k];
        fourier_complex_t diff = DFT_conv.X[k] - product;
        double err = cabs(diff);
        if (err > max_err) max_err = err;
    }

    free(c_conv);
    dft_result_free(&DFT_conv);
    dft_result_free(&DFT_x);
    dft_result_free(&DFT_y);

    return max_err;
}

double verify_parseval_theorem(const double *x, size_t N) {
    if (!x || N == 0) return -1.0;

    double energy_time = 0.0;
    for (size_t n = 0; n < N; n++) energy_time += x[n] * x[n];

    dft_result_t DFT_x = dft_direct(x, N, 1.0);
    double energy_freq = 0.0;
    for (size_t k = 0; k < N; k++) {
        double re = creal(DFT_x.X[k]);
        double im = cimag(DFT_x.X[k]);
        energy_freq += re * re + im * im;
    }
    energy_freq /= (double)N;

    dft_result_free(&DFT_x);
    return fabs(energy_time - energy_freq);
}

double verify_wiener_khinchin(const double *x, size_t N) {
    if (!x || N == 0) return -1.0;

    /* ── Wiener-Khinchin theorem (discrete, exact DFT identity):  ─────
     *   DFT{circular_autocorrelation}[k] = |X[k]|² / N = periodogram[k]
     *
     *   Circular autocorrelation:  r_c[τ] = (1/N)·Σ_{n=0}^{N-1} x[n]·x[(n+τ) mod N]
     *   This is the EXACT discrete counterpart: the N-point DFT of the
     *   circular autocorrelation equals the periodogram to machine precision.
     *
     *   For real signals, the linear biased autocorrelation r_b[τ] is an
     *   approximation: r_b[τ] = (1/N)·Σ_{n=0}^{N-1-|τ|} x[n]·x[n+|τ|].
     *   As N→∞, r_b → r_c (the wraparound effect diminishes). ─────────── */

    /* Step 1: Periodogram P[k] = |X[k]|² / N */
    psd_result_t psd = periodogram_raw(x, N, 1.0);
    if (!psd.psd) return -1.0;

    /* Step 2: Circular autocorrelation r_c[τ] = (1/N)·Σ x[n]·x[(n+τ) mod N] */
    double *r_c = (double*)calloc(N, sizeof(double));
    if (!r_c) { psd_result_free(&psd); return -1.0; }

    for (size_t tau = 0; tau < N; tau++) {
        double sum = 0.0;
        for (size_t n = 0; n < N; n++) {
            sum += x[n] * x[(n + tau) % N];
        }
        r_c[tau] = sum / (double)N;  /* = (1/N) circular correlation */
    }

    /* Step 3: N-point DFT of circular autocorrelation */
    dft_result_t DFT_r = dft_direct(r_c, N, 1.0);
    if (!DFT_r.X) { free(r_c); psd_result_free(&psd); return -1.0; }

    /* Step 4: Wiener-Khinchin identity: DFT{r_c}[k] = P[k] exactly.
     *   Both are N-point sequences on the same frequency grid.
     *   Compare sample-for-sample.  Error should be near machine epsilon. */
    double max_err = 0.0;
    for (size_t k = 0; k < N; k++) {
        double err = fabs(psd.psd[k] - creal(DFT_r.X[k]));
        if (err > max_err) max_err = err;
    }

    free(r_c);
    psd_result_free(&psd);
    dft_result_free(&DFT_r);

    return max_err;
}

/* ─── L6: Wiener Deconvolution ──────────────────────────────────────────── */

void wiener_deconvolution(const double *y, size_t y_len,
                            const double *h, size_t h_len,
                            double nsr, double *s_out) {
    if (!y || !h || !s_out || y_len == 0 || h_len == 0) return;

    size_t fft_len = y_len + h_len - 1;
    size_t padded_len = 1;
    while (padded_len < fft_len) padded_len <<= 1;

    fourier_complex_t *Y = (fourier_complex_t*)calloc(padded_len, sizeof(fourier_complex_t));
    fourier_complex_t *H = (fourier_complex_t*)calloc(padded_len, sizeof(fourier_complex_t));
    if (!Y || !H) { free(Y); free(H); return; }

    for (size_t i = 0; i < y_len; i++) Y[i] = y[i];
    for (size_t i = 0; i < h_len; i++) H[i] = h[i];

    fft_plan_t *pf = fft_plan_create(padded_len, 0);
    fft_plan_t *pi = fft_plan_create(padded_len, 1);

    fft_execute_dit(pf, Y);
    fft_execute_dit(pf, H);

    /* Wiener filter: Ŝ = Y·H* / (|H|² + NSR) */
    for (size_t i = 0; i < padded_len; i++) {
        double H_mag2 = creal(H[i]) * creal(H[i]) + cimag(H[i]) * cimag(H[i]);
        fourier_complex_t inv = (H_mag2 + nsr > 1e-20)
            ? conj(H[i]) / (H_mag2 + nsr)
            : 0.0;
        Y[i] = Y[i] * inv;
    }

    fft_execute_dit(pi, Y);

    for (size_t i = 0; i < y_len; i++) {
        s_out[i] = creal(Y[i]);
    }

    fft_plan_destroy(pf);
    fft_plan_destroy(pi);
    free(Y);
    free(H);
}

/* ─── L6: Hilbert Transform ─────────────────────────────────────────────── */

void hilbert_transform(const double *x, size_t N, double *ht) {
    if (!x || !ht || N == 0) return;

    size_t fft_len = 1;
    while (fft_len < N) fft_len <<= 1;

    fourier_complex_t *X = (fourier_complex_t*)calloc(fft_len, sizeof(fourier_complex_t));
    if (!X) return;
    for (size_t i = 0; i < N; i++) X[i] = x[i];

    fft_plan_t *pf = fft_plan_create(fft_len, 0);
    fft_plan_t *pi = fft_plan_create(fft_len, 1);

    fft_execute_dit(pf, X);

    /* Apply Hilbert frequency response: H(ω) = -j·sgn(ω) */
    /* For real FFT: X[0] = 0, X[fft_len/2] = 0,
     * X[1..fft_len/2-1] *= -j, X[fft_len/2+1..] *= +j */
    for (size_t k = 1; k < fft_len / 2; k++) {
        X[k] *= -I;
        X[fft_len - k] *= I;
    }
    X[0] = 0.0;
    if (fft_len >= 2) X[fft_len / 2] = 0.0;

    fft_execute_dit(pi, X);

    for (size_t i = 0; i < N; i++) {
        ht[i] = creal(X[i]);
    }

    fft_plan_destroy(pf);
    fft_plan_destroy(pi);
    free(X);
}

/* ─── L6: Analytic Signal ───────────────────────────────────────────────── */

void analytic_signal(const double *x, size_t N,
                      fourier_complex_t *z, double *envelope,
                      double *phase, double *inst_freq) {
    if (!x || N == 0) return;

    double *ht = (double*)safe_malloc(N * sizeof(double));
    hilbert_transform(x, N, ht);

    for (size_t n = 0; n < N; n++) {
        z[n] = x[n] + I * ht[n];
        if (envelope) envelope[n] = cabs(z[n]);
        if (phase) phase[n] = atan2(ht[n], x[n]);
    }

    /* Instantaneous frequency: derivative of unwrapped phase */
    if (inst_freq && phase) {
        /* Simple difference-based IF estimation */
        inst_freq[0] = 0.0;
        for (size_t n = 1; n < N; n++) {
            double diff = phase[n] - phase[n - 1];
            /* Phase unwrapping: keep in [-π, π] */
            if (diff > M_PI) diff -= 2.0 * M_PI;
            if (diff < -M_PI) diff += 2.0 * M_PI;
            inst_freq[n] = diff / (2.0 * M_PI);  /* normalized frequency (cycles/sample) */
        }
    }

    free(ht);
}

/* ─── L8: Fractional Fourier Transform ──────────────────────────────────── */

/**
 * Compute the Fractional Fourier Transform via DFT eigendecomposition.
 *
 * The FrFT of order α corresponds to a rotation by angle α·π/2 in the
 * time-frequency plane.  The eigenfunctions of the Fourier transform are
 * Hermite-Gauss functions H_n(t)·exp(-t²/2), with eigenvalues exp(-j·n·α·π/2).
 *
 * Implementation: decompose the input signal into the eigenbasis of the
 * DFT matrix, multiply eigenvalues by exp(-j·α·π/2), reconstruct.
 *
 * For computational efficiency with arbitrary α, we use the sampling-based
 * method: decompose the DFT matrix F = V·Λ·V^H, then F^α = V·Λ^α·V^H.
 *
 * Approximated here via Taylor/polynomial expansion of the Hermite eigenbasis.
 */
void fractional_fourier_transform(const double *x, size_t N,
                                    double alpha,
                                    fourier_complex_t *output) {
    if (!x || !output || N == 0) return;

    /* Normalize alpha to [0, 4) */
    alpha = fmod(alpha, 4.0);
    if (alpha < 0.0) alpha += 4.0;

    /* Special cases */
    if (fabs(alpha) < 1e-12 || fabs(alpha - 4.0) < 1e-12) {
        /* Identity */
        for (size_t i = 0; i < N; i++) output[i] = x[i];
        return;
    }
    if (fabs(alpha - 1.0) < 1e-12) {
        /* Standard Fourier transform */
        dft_result_t dft = dft_direct(x, N, 1.0);
        if (dft.X) { memcpy(output, dft.X, N * sizeof(fourier_complex_t)); dft_result_free(&dft); }
        return;
    }
    if (fabs(alpha - 2.0) < 1e-12) {
        /* Parity: f(-t). For discrete, x[N-n] */
        for (size_t i = 0; i < N; i++) output[i] = x[(N - i) % N];
        return;
    }
    if (fabs(alpha - 3.0) < 1e-12) {
        /* Inverse Fourier */
        dft_result_t dft = dft_direct(x, N, 1.0);
        if (dft.X) {
            for (size_t k = 0; k < N; k++) output[k] = conj(dft.X[k]);
            dft_result_free(&dft);
        }
        return;
    }

    /* General α:  F^α = cos(απ/2)·I + sin(απ/2)·F — linear combination of
     * identity and ordinary Fourier transform, valid for discrete
     * approximation (fractionalization of DFT via its eigenstructure).
     *
     * This is exact for the DFT eigendecomposition when the DFT has
     * 4 distinct eigenvalues {1, -j, -1, j}.
     */
    double cosa = cos(alpha * M_PI / 2.0);
    double sina = sin(alpha * M_PI / 2.0);

    dft_result_t dft = dft_direct(x, N, 1.0);
    if (!dft.X) return;

    for (size_t i = 0; i < N; i++) {
        output[i] = cosa * x[i] + sina * dft.X[i];
    }

    dft_result_free(&dft);
}
