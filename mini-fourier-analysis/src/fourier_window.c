/**
 * @file    fourier_window.c
 * @brief   Window Functions for Spectral Analysis — L2 + L5
 *
 * @details Implements all major window functions used in DFT-based
 *          spectral analysis.  Each window balances the trade-off
 *          between mainlobe width (frequency resolution) and sidelobe
 *          level (spectral leakage suppression).
 *
 *          Reference: Harris, F.J., "On the Use of Windows for Harmonic
 *          Analysis with the Discrete Fourier Transform", Proc. IEEE 66(1),
 *          1978, pp.51-83.
 *
 * Knowledge covered:
 *   L5: Hann, Hamming, Blackman, Blackman-Harris, Nuttall, flat-top,
 *       Kaiser, Bartlett, Gaussian, Dolph-Chebyshev, Tukey, Lanczos
 *   L2: Coherent gain, ENBW, processing loss, scalloping loss
 *   L4: Dirichlet kernel, uncertainty principle
 */

#include "fourier_window.h"
#include "fourier_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ─── Internal helpers ──────────────────────────────────────────────────── */

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "fourier_window: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/** Modified Bessel function of the first kind I₀(x) — power series */
static double bessel_I0(double x) {
    double sum = 1.0;
    double term = 1.0;
    double x2 = x * x / 4.0;
    for (int k = 1; k < 50; k++) {
        term *= x2 / ((double)k * (double)k);
        sum += term;
        if (term < 1e-15 * sum) break;
    }
    return sum;
}

/** sinc function: sin(π·x)/(π·x), with limit 1 at x=0 */
static double sinc_normalized(double x) {
    if (fabs(x) < 1e-12) return 1.0;
    double arg = M_PI * x;
    return sin(arg) / arg;
}

/** Initialize window_config_t common fields */
static void window_init(window_config_t *wc, window_type_t type, size_t N) {
    memset(wc, 0, sizeof(*wc));
    wc->type = type;
    wc->length = N;
    wc->coefficients = (double*)safe_malloc(N * sizeof(double));
}

/* ─── L5: Window Generation ─────────────────────────────────────────────── */

window_config_t window_rectangular(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_RECTANGULAR, N);
    for (size_t i = 0; i < N; i++) wc.coefficients[i] = 1.0;
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_hann(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_HANN, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        wc.coefficients[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)(N - 1)));
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_hamming(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_HAMMING, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        wc.coefficients[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)i / (double)(N - 1));
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_blackman(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_BLACKMAN, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        double a = 2.0 * M_PI * (double)i / (double)(N - 1);
        wc.coefficients[i] = 0.42 - 0.5 * cos(a) + 0.08 * cos(2.0 * a);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_blackman_harris4(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_BLACKMAN, N);  /* reuse type enum; note: new type could be added */

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        double a = 2.0 * M_PI * (double)i / (double)N;
        wc.coefficients[i] = 0.35875
                           - 0.48829 * cos(a)
                           + 0.14128 * cos(2.0 * a)
                           - 0.01168 * cos(3.0 * a);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_nuttall4(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_BLACKMAN, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        double a = 2.0 * M_PI * (double)i / (double)N;
        wc.coefficients[i] = 0.3635819
                           - 0.4891775 * cos(a)
                           + 0.1365995 * cos(2.0 * a)
                           - 0.0106411 * cos(3.0 * a);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_flattop(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_FLATTOP, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        double a = 2.0 * M_PI * (double)i / (double)(N - 1);
        wc.coefficients[i] = 0.21557895
                           - 0.41663158 * cos(a)
                           + 0.277263158 * cos(2.0 * a)
                           - 0.083578947 * cos(3.0 * a)
                           + 0.006947368 * cos(4.0 * a);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_kaiser(size_t N, double beta) {
    window_config_t wc;
    window_init(&wc, WINDOW_KAISER, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    double denom = bessel_I0(beta);
    if (denom < 1e-20) denom = 1.0;  /* safety */

    for (size_t i = 0; i < N; i++) {
        double x = (2.0 * (double)i / (double)(N - 1)) - 1.0;
        double arg = beta * sqrt(1.0 - x * x);
        wc.coefficients[i] = bessel_I0(arg) / denom;
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_bartlett(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_BARTLETT, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    double half = (double)(N - 1) / 2.0;
    for (size_t i = 0; i < N; i++) {
        if ((double)i <= half) {
            wc.coefficients[i] = (2.0 * (double)i) / (double)(N - 1);
        } else {
            wc.coefficients[i] = 2.0 - (2.0 * (double)i) / (double)(N - 1);
        }
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_gaussian(size_t N, double alpha) {
    window_config_t wc;
    window_init(&wc, WINDOW_GAUSSIAN, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    double center = (double)(N - 1) / 2.0;
    double sigma = (double)(N - 1) / (2.0 * alpha);
    for (size_t i = 0; i < N; i++) {
        double x = ((double)i - center) / sigma;
        wc.coefficients[i] = exp(-0.5 * x * x);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_dolph_chebyshev(size_t N, double sidelobe_db) {
    window_config_t wc;
    window_init(&wc, WINDOW_RECTANGULAR, N);  /* no enum for DC yet */

    if (N < 2 || sidelobe_db <= 0.0) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    /* Dolph-Chebyshev window via frequency-domain sampling & IDFT.
     * The window is designed such that its DTFT has equiripple sidelobes
     * at the specified level.
     *
     * Steps:
     * 1. Define M = N (odd-length equivalent)
     * 2. Compute x₀ = cosh(acosh(10^{sidelobe_db/20})/(M-1))
     * 3. W(ω_k) = (-1)^k · T_{M-1}(x₀·cos(π·k/(M-1))) / 10^{sidelobe_db/20}
     *    for k = 0..M-1
     * 4. IDFT to get time-domain window
     */

    int M = (int)N;
    double atten_linear = pow(10.0, sidelobe_db / 20.0);

    /* x₀ = cosh(acosh(R)/(M-1)) where R = 10^{sidelobe/20} */
    double R = atten_linear;
    double acosh_R = log(R + sqrt(R * R - 1.0));  /* acosh(x) = ln(x + √(x²-1)) */
    double x0 = cosh(acosh_R / (double)(M - 1));

    /* Compute frequency response samples */
    fourier_complex_t *W_freq = (fourier_complex_t*)safe_malloc(N * sizeof(fourier_complex_t));
    for (int k = 0; k < M; k++) {
        double omega = M_PI * (double)k / (double)(M - 1);
        double xn = x0 * cos(omega);
        /* T_{M-1}(xn) clipped to avoid blow-up */
        double T_val;
        if (fabs(xn) <= 1.0) {
            T_val = cos((double)(M - 1) * acos(xn));
        } else {
            /* cosh formula for |x| > 1 */
            double acosh_val = log(fabs(xn) + sqrt(xn * xn - 1.0));
            T_val = cosh((double)(M - 1) * acosh_val);
            if (xn < 0 && ((M - 1) % 2 == 1)) T_val = -T_val;
        }
        /* W[k] = (-1)^k · T_val / R */
        double sign = (k % 2 == 0) ? 1.0 : -1.0;
        W_freq[k] = sign * T_val / R;
    }

    /* IDFT of W_freq (size N) to get time-domain window.
     * Use N-point DFT (treating samples 0..N-1 as spectrum).
     * But the DC window is derived from DFT of frequency sampling.
     * For an N-point symmetric window from Chebyshev design:
     * w[n] = (1/N)·Σ_{k=0}^{N-1} W[k]·exp(j·2π·k·(n - (N-1)/2)/N)
     *
     * Actually compute via DFT of W_freq treated as spectrum. */
    fourier_complex_t *w_complex = (fourier_complex_t*)safe_malloc(N * sizeof(fourier_complex_t));
    /* X[k] = W_freq[k]; compute x[n] via IDFT */
    for (size_t n = 0; n < N; n++) {
        fourier_complex_t sum = 0.0;
        for (size_t k = 0; k < N; k++) {
            double angle = 2.0 * M_PI * (double)k * (double)n / (double)N;
            sum += W_freq[k] * (cos(angle) + I * sin(angle));
        }
        w_complex[n] = sum / (double)N;
    }

    /* Extract real part and circular-shift to center */
    /* The resulting window should be symmetric */
    for (size_t i = 0; i < N; i++) {
        size_t idx = (i + N / 2) % N;
        wc.coefficients[i] = creal(w_complex[idx]);
    }

    /* Normalize to peak = 1 */
    double wmax = 0.0;
    for (size_t i = 0; i < N; i++) {
        if (wc.coefficients[i] > wmax) wmax = wc.coefficients[i];
    }
    if (wmax > 1e-20) {
        for (size_t i = 0; i < N; i++) wc.coefficients[i] /= wmax;
    }

    free(W_freq);
    free(w_complex);
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_tukey(size_t N, double alpha) {
    window_config_t wc;
    window_init(&wc, WINDOW_HANN, N);

    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    if (N < 2 || alpha <= 0.0) {
        for (size_t i = 0; i < N; i++) wc.coefficients[i] = 1.0;
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    double taper_len = alpha * (double)(N - 1) / 2.0;
    for (size_t i = 0; i < N; i++) {
        if ((double)i < taper_len) {
            wc.coefficients[i] = 0.5 * (1.0 + cos(M_PI * ((double)i / taper_len - 1.0)));
        } else if ((double)i >= (double)(N - 1) - taper_len) {
            double idx = (double)i - ((double)(N - 1) - taper_len);
            wc.coefficients[i] = 0.5 * (1.0 + cos(M_PI * (idx / taper_len)));
        } else {
            wc.coefficients[i] = 1.0;
        }
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

window_config_t window_lanczos(size_t N) {
    window_config_t wc;
    window_init(&wc, WINDOW_RECTANGULAR, N);

    if (N < 2) {
        if (N == 1) wc.coefficients[0] = 1.0;
        window_compute_metrics(wc.coefficients, N, &wc);
        return wc;
    }

    for (size_t i = 0; i < N; i++) {
        double x = 2.0 * (double)i / (double)(N - 1) - 1.0;
        wc.coefficients[i] = sinc_normalized(x);
    }
    window_compute_metrics(wc.coefficients, N, &wc);
    return wc;
}

/* ─── L2: Window Metrics Computation ────────────────────────────────────── */

void window_compute_metrics(double *coeffs, size_t N, window_config_t *config) {
    if (!coeffs || N == 0 || !config) return;

    /* Coherent gain: Σ w[n] / N */
    double sum_w = 0.0, sum_w2 = 0.0;
    for (size_t i = 0; i < N; i++) {
        sum_w += coeffs[i];
        sum_w2 += coeffs[i] * coeffs[i];
    }
    config->coherent_gain = sum_w / (double)N;

    /* Equivalent noise bandwidth: N·Σ w²[n] / (Σ w[n])² */
    if (sum_w != 0.0) {
        config->equivalent_noise_bw = (double)N * sum_w2 / (sum_w * sum_w);
    } else {
        config->equivalent_noise_bw = INFINITY;
    }

    /* Processing loss (due to coherent gain + ENBW):
     * PL = 10·log₁₀(ENBW) - 20·log₁₀(coherent_gain) */
    if (config->coherent_gain > 1e-20 && isfinite(config->equivalent_noise_bw)) {
        config->processing_loss_db = 10.0 * log10(config->equivalent_noise_bw)
                                     - 20.0 * log10(config->coherent_gain);
    }

    /* Scalloping loss: worst-case amplitude error when a sinusoid falls
     * exactly between two DFT bins (half-bin offset).
     * Approximate via the window's DTFT at ω = π/N. */
    fourier_complex_t W_half = 0.0;
    double omega_half = M_PI / (double)N;
    for (size_t n = 0; n < N; n++) {
        double angle = -omega_half * (double)n;
        W_half += coeffs[n] * (cos(angle) + I * sin(angle));
    }
    double W_peak = sum_w;  /* DC gain */
    if (W_peak > 1e-20) {
        config->scalloping_loss_db = 20.0 * log10(W_peak / cabs(W_half));
    }

    /* Sidelobe level and mainlobe width: approximate via oversampled
     * frequency response (zero-padded FFT equivalent).
     * We sample at 16× oversampling to find sidelobe peak. */
    size_t oversample = 16;
    size_t M = N * oversample;
    double *W_amp = (double*)malloc(M * sizeof(double));
    if (W_amp) {
        /* Evaluate DTFT at M points over [0, π] */
        double max_mainlobe = 0.0;
        int in_mainlobe = 1;
        double peak_after_main = 0.0;
        double bw_low = -1.0, bw_high = -1.0;

        for (size_t k = 0; k <= M / 2; k++) {
            double omega = 2.0 * M_PI * (double)k / (double)M;
            fourier_complex_t W_val = 0.0;
            for (size_t n = 0; n < N; n++) {
                double angle = -omega * (double)n;
                W_val += coeffs[n] * (cos(angle) + I * sin(angle));
            }
            W_amp[k] = cabs(W_val);

            /* Find mainlobe maximum */
            if (W_amp[k] > max_mainlobe) max_mainlobe = W_amp[k];

            /* Simple mainlobe tracking: scan until first valley */
            if (in_mainlobe && k > 0 && W_amp[k] < W_amp[k - 1] && W_amp[k] < 0.1 * max_mainlobe) {
                in_mainlobe = 0;
            }
            if (!in_mainlobe && W_amp[k] > peak_after_main) {
                peak_after_main = W_amp[k];
            }

            /* 3 dB width */
            if (max_mainlobe > 1e-20) {
                if (bw_low < 0 && W_amp[k] >= max_mainlobe / sqrt(2.0)) bw_low = omega;
                if (bw_low >= 0 && bw_high < 0 && W_amp[k] < max_mainlobe / sqrt(2.0)) {
                    bw_high = omega;
                }
            }
        }

        if (max_mainlobe > 1e-20) {
            /* Sidelobe level relative to mainlobe */
            if (peak_after_main > 1e-30) {
                config->sidelobe_level_db = 20.0 * log10(peak_after_main / max_mainlobe);
            } else {
                config->sidelobe_level_db = -400.0;
            }

            /* Mainlobe 3 dB width in bins: bw_normalized × N / (2π) */
            if (bw_low >= 0 && bw_high > bw_low) {
                double bw_rad = 2.0 * (bw_high - bw_low);
                config->mainlobe_3dB_width = bw_rad * (double)N / (2.0 * M_PI);
            } else {
                config->mainlobe_3dB_width = config->equivalent_noise_bw;  /* fallback */
            }
        }

        free(W_amp);
    }
}

/* ─── L2: Window Application ────────────────────────────────────────────── */

void window_apply(const double *x, const double *window, size_t N, double *y) {
    if (!x || !window || !y || N == 0) return;
    for (size_t i = 0; i < N; i++) {
        y[i] = x[i] * window[i];
    }
}

/* ─── L2: Window Frequency Response ─────────────────────────────────────── */

void window_frequency_response(const double *window, size_t N,
                                 size_t freq_points,
                                 fourier_complex_t *W_out, double *omega_grid) {
    if (!window || N == 0 || !W_out || !omega_grid || freq_points == 0) return;

    for (size_t m = 0; m < freq_points; m++) {
        double omega = M_PI * (double)m / (double)(freq_points - 1);
        omega_grid[m] = omega;
        fourier_complex_t W = 0.0;
        for (size_t n = 0; n < N; n++) {
            double angle = -omega * (double)n;
            W += window[n] * (cos(angle) + I * sin(angle));
        }
        W_out[m] = W;
    }
}

/* ─── L5: Kaiser Beta Estimation ────────────────────────────────────────── */

double window_kaiser_beta_from_attenuation(double attenuation_db) {
    if (attenuation_db > 50.0) {
        return 0.1102 * (attenuation_db - 8.7);
    } else if (attenuation_db >= 21.0) {
        return 0.5842 * pow(attenuation_db - 21.0, 0.4)
               + 0.07886 * (attenuation_db - 21.0);
    } else {
        return 0.0;
    }
}

/* ─── L5: Window Length Recommendation ──────────────────────────────────── */

size_t window_recommend_length(window_type_t window_type,
                                 double mainlobe_width_hz, double fs) {
    if (fs <= 0.0 || mainlobe_width_hz <= 0.0) return 256;

    /* Mainlobe 3dB width in bins for each window type (approximate) */
    double K;
    switch (window_type) {
        case WINDOW_RECTANGULAR: K = 0.89; break;
        case WINDOW_HANN:        K = 1.44; break;
        case WINDOW_HAMMING:     K = 1.30; break;
        case WINDOW_BLACKMAN:    K = 1.68; break;
        case WINDOW_BARTLETT:    K = 1.28; break;
        case WINDOW_FLATTOP:     K = 3.72; break;
        default:                 K = 1.50; break;
    }

    double N_double = K * fs / mainlobe_width_hz;
    size_t N = (size_t)ceil(N_double);
    if (N < 4) N = 4;

    /* Round up to even (better for most applications) */
    if (N % 2 != 0) N++;
    return N;
}

/* ─── Memory Management ─────────────────────────────────────────────────── */

void window_config_free(window_config_t *wc) {
    if (wc) {
        free(wc->coefficients);
        wc->coefficients = NULL;
        wc->length = 0;
    }
}
