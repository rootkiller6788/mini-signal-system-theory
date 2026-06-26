/**
 * @file convolution.c
 * @brief L2/L5 Convolution Algorithm Implementations
 *
 * Implements all convolution variants: CT numerical, DT direct, circular,
 * overlap-add, overlap-save, deconvolution, and correlation.
 *
 * Knowledge Implementation:
 *   L2-I01: CT convolution (trapezoidal numerical integration)
 *   L2-I02: CT convolution fast (subsampled IR)
 *   L2-I03: DT linear convolution (direct O(N^2))
 *   L2-I04: DT linear convolution (FFT-based O(N log N))
 *   L2-I05: DT circular convolution (direct)
 *   L2-I06: DT circular convolution (FFT-based via DFT multiplication)
 *   L2-I07: Overlap-add block convolution
 *   L2-I08: Overlap-save block convolution
 *   L2-I09: Deconvolution via frequency-domain division
 *   L2-I10: Cross-correlation
 *   L2-I11: Auto-correlation
 *   L2-I12: Normalized cross-correlation (Pearson-like)
 *   L2-I13: Step response from impulse response (cumulative sum)
 *   L2-I14: Impulse response from step response (difference)
 *   L2-I15: Convolution commutativity verification
 *   L2-I16: Zero-padding utility
 *   L2-I17: DFT computation (simple O(N^2) for internal use)
 *   L2-I18: IDFT computation (simple O(N^2) for internal use)
 *
 * References:
 *   Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd Ed., Ch. 8
 *   Proakis & Manolakis, Digital Signal Processing, 4th Ed., Ch. 5
 */

#include "convolution.h"
#include "system_defs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 *  INTERNAL: Simple DFT / IDFT (O(N^2)) for convolution
 *  NOTE: Production code would use FFT. This DFT is for completeness
 *  of the convolution algorithms.
 * ======================================================================== */

static int dft_direct(const double *x, int N, double complex *X)
{
    if (!x || !X || N <= 0) return -1;

    for (int k = 0; k < N; k++) {
        X[k] = 0.0 + 0.0 * I;
        for (int n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * (double)k * (double)n / (double)N;
            X[k] += x[n] * (cos(angle) + sin(angle) * I);
        }
    }
    return 0;
}

static int idft_direct(const double complex *X, int N, double *x)
{
    if (!X || !x || N <= 0) return -1;

    for (int n = 0; n < N; n++) {
        double complex sum = 0.0 + 0.0 * I;
        for (int k = 0; k < N; k++) {
            double angle = 2.0 * M_PI * (double)k * (double)n / (double)N;
            sum += X[k] * (cos(angle) + sin(angle) * I);
        }
        x[n] = creal(sum) / (double)N;
    }
    return 0;
}

static int next_power_of_two(int n)
{
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

/* ========================================================================
 *  L2-I01: CT Convolution — Direct Trapezoidal Integration
 *
 *  Theorem (Convolution integral):
 *    y(t) = x(t) * h(t) = integral_{-inf}^{inf} x(tau) h(t-tau) dtau
 *
 *  For causal h(t) with x(t)=0 for t<0:
 *    y(t) = integral_{0}^{t} x(tau) h(t-tau) dtau
 *
 *  Complexity: O(N_x * N_h), where N_x = len(x), N_h = len(h)
 * ======================================================================== */

int ct_convolution_direct(const ct_signal_t *x, const ct_signal_t *h,
                          ct_signal_t *result)
{
    if (!x || !h || !result) return -1;
    if (!x->samples || !h->samples) return -1;
    if (x->sample_rate != h->sample_rate) return -1;

    double dt = 1.0 / x->sample_rate;
    size_t Nx = x->num_samples;
    size_t Nh = h->num_samples;
    size_t Ny = Nx + Nh - 1;

    if (Ny > result->num_samples) return -1;

    double duration = (double)(Ny - 1) * dt;
    result->t_start = x->t_start + h->t_start;
    result->t_end = result->t_start + duration;

    /* For each output time t_y[k], compute convolution integral */
    for (size_t k = 0; k < Ny; k++) {
        double t = result->t_start + (double)k * dt;
        double sum = 0.0;

        /* integration over tau where both x(tau) and h(t-tau) are non-zero */
        for (size_t i = 0; i < Nx; i++) {
            double tau = x->t_start + (double)i * dt;
            double t_minus_tau = t - tau;

            /* Find h(t-tau) by interpolation */
            double h_idx_f = (t_minus_tau - h->t_start) / dt;
            if (h_idx_f < 0.0 || h_idx_f >= (double)(Nh - 1)) continue;

            size_t h_idx = (size_t)h_idx_f;
            double frac = h_idx_f - (double)h_idx;
            double h_val;
            if (h_idx + 1 < Nh) {
                h_val = h->samples[h_idx] * (1.0 - frac) +
                        h->samples[h_idx + 1] * frac;
            } else {
                h_val = h->samples[Nh - 1];
            }

            sum += x->samples[i] * h_val * dt;
        }
        result->samples[k] = sum;
    }
    return 0;
}

/* ========================================================================
 *  L2-I02: CT Convolution — Fast (Subsampled IR)
 *
 *  When h(t) varies slowly relative to x(t), we can decimate h(t)
 *  to reduce computation. Uses piecewise-constant approximation.
 * ======================================================================== */

int ct_convolution_fast(const ct_signal_t *x, const ct_signal_t *h,
                        ct_signal_t *result, int decimation_factor)
{
    if (!x || !h || !result || decimation_factor < 1) return -1;

    /* Build decimated h */
    size_t Nh_dec = (h->num_samples + decimation_factor - 1) / decimation_factor;
    ct_signal_t h_dec = ct_signal_alloc(Nh_dec, h->t_start, h->t_end,
                                         h->sample_rate / decimation_factor);
    if (!h_dec.samples) return -1;

    for (size_t i = 0; i < Nh_dec; i++) {
        size_t src_idx = i * decimation_factor;
        if (src_idx < h->num_samples) {
            h_dec.samples[i] = h->samples[src_idx];
        }
    }

    int ret = ct_convolution_direct(x, &h_dec, result);
    ct_signal_free(&h_dec);
    return ret;
}

/* ========================================================================
 *  L2-I03: DT Linear Convolution — Direct O(N^2)
 *
 *  Definition: y[n] = sum_{k=0}^{N-1} x[k] * h[n-k]
 *  Output length: Nx + Nh - 1
 * ======================================================================== */

int dt_convolution_direct(const dt_signal_t *x, const dt_signal_t *h,
                          dt_signal_t *y)
{
    if (!x || !h || !y) return -1;
    if (!x->samples || !h->samples) return -1;

    size_t Nx = x->num_samples;
    size_t Nh = h->num_samples;
    size_t Ny = Nx + Nh - 1;

    if (y->num_samples < Ny) return -1;

    /* Zero output */
    memset(y->samples, 0, y->num_samples * sizeof(double));

    /* Convolution sum */
    for (size_t n = 0; n < Ny; n++) {
        double sum = 0.0;
        size_t k_start = (n >= Nh) ? (n - Nh + 1) : 0;
        size_t k_end = (n < Nx) ? n : (Nx - 1);

        for (size_t k = k_start; k <= k_end; k++) {
            sum += x->samples[k] * h->samples[n - k];
        }
        y->samples[n] = sum;
    }
    return 0;
}

/* ========================================================================
 *  L2-I04: DT Linear Convolution — FFT-based O(N log N)
 *
 *  Uses DFT multiplication: y = IDFT{ DFT{x_zp} * DFT{h_zp} }
 *  where x_zp and h_zp are zero-padded to length N >= Nx+Nh-1.
 *
 *  Theorem (Convolution Theorem):
 *    x[n] * h[n] <-> X(e^{jw}) * H(e^{jw}) in frequency domain
 * ======================================================================== */

int dt_convolution_fft(const dt_signal_t *x, const dt_signal_t *h,
                       dt_signal_t *y)
{
    if (!x || !h || !y) return -1;
    if (!x->samples || !h->samples) return -1;

    size_t Nx = x->num_samples;
    size_t Nh = h->num_samples;
    size_t Nfft = (size_t)next_power_of_two((int)(Nx + Nh - 1));

    if (Nfft > SIG_MAX_SAMPLES) return -1;

    /* Zero-pad inputs */
    double *x_pad = (double *)calloc(Nfft, sizeof(double));
    double *h_pad = (double *)calloc(Nfft, sizeof(double));
    double complex *X = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *H = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *Y = (double complex *)calloc(Nfft, sizeof(double complex));

    if (!x_pad || !h_pad || !X || !H || !Y) {
        free(x_pad); free(h_pad); free(X); free(H); free(Y);
        return -1;
    }

    memcpy(x_pad, x->samples, Nx * sizeof(double));
    memcpy(h_pad, h->samples, Nh * sizeof(double));

    /* DFT */
    dft_direct(x_pad, (int)Nfft, X);
    dft_direct(h_pad, (int)Nfft, H);

    /* Multiply in frequency domain */
    for (size_t k = 0; k < Nfft; k++) {
        Y[k] = X[k] * H[k];
    }

    /* IDFT */
    double *y_full = (double *)calloc(Nfft, sizeof(double));
    if (!y_full) {
        free(x_pad); free(h_pad); free(X); free(H); free(Y);
        return -1;
    }

    idft_direct(Y, (int)Nfft, y_full);

    /* Copy valid portion */
    size_t Ny = Nx + Nh - 1;
    if (y->num_samples < Ny) Ny = y->num_samples;
    memcpy(y->samples, y_full, Ny * sizeof(double));

    free(x_pad); free(h_pad); free(X); free(H); free(Y); free(y_full);
    return 0;
}

/* ========================================================================
 *  L2-I05: DT Circular Convolution — Direct O(N^2)
 *
 *  y[n] = sum_{k=0}^{N-1} x[k] * h[(n-k) mod N]
 *  Both x and h must have the same length N.
 * ======================================================================== */

int dt_circular_convolution_direct(const dt_signal_t *x, const dt_signal_t *h,
                                    dt_signal_t *y)
{
    if (!x || !h || !y) return -1;
    if (x->num_samples != h->num_samples) return -1;

    size_t N = x->num_samples;
    if (y->num_samples < N) return -1;

    for (size_t n = 0; n < N; n++) {
        double sum = 0.0;
        for (size_t k = 0; k < N; k++) {
            size_t idx = (n >= k) ? (n - k) : (n - k + N);
            sum += x->samples[k] * h->samples[idx];
        }
        y->samples[n] = sum;
    }
    return 0;
}

/* ========================================================================
 *  L2-I06: DT Circular Convolution — FFT-based
 * ======================================================================== */

int dt_circular_convolution_fft(const dt_signal_t *x, const dt_signal_t *h,
                                 dt_signal_t *y)
{
    if (!x || !h || !y) return -1;
    if (x->num_samples != h->num_samples) return -1;

    int N = (int)x->num_samples;
    if (y->num_samples < (size_t)N) return -1;

    double complex *X = (double complex *)calloc((size_t)N, sizeof(double complex));
    double complex *H = (double complex *)calloc((size_t)N, sizeof(double complex));
    double complex *Y = (double complex *)calloc((size_t)N, sizeof(double complex));

    if (!X || !H || !Y) {
        free(X); free(H); free(Y);
        return -1;
    }

    dft_direct(x->samples, N, X);
    dft_direct(h->samples, N, H);

    for (int k = 0; k < N; k++) {
        Y[k] = X[k] * H[k];
    }

    idft_direct(Y, N, y->samples);

    free(X); free(H); free(Y);
    return 0;
}

/* ========================================================================
 *  L2-I07: Overlap-Add Block Convolution
 *
 *  Algorithm:
 *    1. Split x[n] into blocks of length L
 *    2. Zero-pad each block to N = L + M - 1 (where M = len(h))
 *    3. FFT-convolve each block with h[n]
 *    4. Add overlapping regions
 *
 *  Complexity: O(Nx * log M) when L >> M
 * ======================================================================== */

int dt_overlap_add(const dt_signal_t *x, const dt_signal_t *h,
                   int block_size, dt_signal_t *y)
{
    if (!x || !h || !y || block_size < 1) return -1;

    size_t Nx = x->num_samples;
    size_t M = h->num_samples;
    size_t Nfft = (size_t)next_power_of_two(block_size + (int)M - 1);

    if (Nfft > SIG_MAX_SAMPLES) return -1;

    size_t Ny = Nx + M - 1;
    if (y->num_samples < Ny) return -1;
    memset(y->samples, 0, y->num_samples * sizeof(double));

    double *block = (double *)calloc(Nfft, sizeof(double));
    double *y_block = (double *)calloc(Nfft, sizeof(double));
    double complex *Xf = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *Hf = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *Yf = (double complex *)calloc(Nfft, sizeof(double complex));

    if (!block || !y_block || !Xf || !Hf || !Yf) {
        free(block); free(y_block); free(Xf); free(Hf); free(Yf);
        return -1;
    }

    /* Pre-compute DFT of h (zero-padded to Nfft) */
    double *h_pad = (double *)calloc(Nfft, sizeof(double));
    if (!h_pad) {
        free(block); free(y_block); free(Xf); free(Hf); free(Yf);
        return -1;
    }
    memcpy(h_pad, h->samples, M * sizeof(double));
    dft_direct(h_pad, (int)Nfft, Hf);
    free(h_pad);

    /* Process blocks */
    for (size_t offset = 0; offset < Nx; offset += (size_t)block_size) {
        /* Extract block */
        memset(block, 0, Nfft * sizeof(double));
        size_t copy_len = (offset + (size_t)block_size <= Nx)
                          ? (size_t)block_size : (Nx - offset);
        memcpy(block, x->samples + offset, copy_len * sizeof(double));

        /* DFT -> multiply -> IDFT */
        dft_direct(block, (int)Nfft, Xf);
        for (size_t k = 0; k < Nfft; k++) {
            Yf[k] = Xf[k] * Hf[k];
        }
        idft_direct(Yf, (int)Nfft, y_block);

        /* Add to output (overlap-add) */
        for (size_t i = 0; i < Nfft && (offset + i) < Ny; i++) {
            y->samples[offset + i] += y_block[i];
        }
    }

    free(block); free(y_block); free(Xf); free(Hf); free(Yf);
    return 0;
}

/* ========================================================================
 *  L2-I08: Overlap-Save Block Convolution
 *
 *  Algorithm:
 *    1. Use blocks of length N (FFT size)
 *    2. Each block overlaps previous by M-1 samples
 *    3. After IDFT, discard first M-1 samples (circular wrap-around)
 *    4. Concatenate valid portions
 * ======================================================================== */

int dt_overlap_save(const dt_signal_t *x, const dt_signal_t *h,
                    int block_size, dt_signal_t *y)
{
    if (!x || !h || !y || block_size < 1) return -1;

    size_t Nx = x->num_samples;
    size_t M = h->num_samples;
    size_t N = (size_t)next_power_of_two(block_size);

    if (N < M || N > SIG_MAX_SAMPLES) return -1;

    size_t valid_per_block = N - M + 1;
    size_t Ny = Nx + M - 1;
    if (y->num_samples < Ny) return -1;
    memset(y->samples, 0, y->num_samples * sizeof(double));

    double *block_buf = (double *)calloc(N, sizeof(double));
    double *y_full = (double *)calloc(N, sizeof(double));
    double complex *Xf = (double complex *)calloc(N, sizeof(double complex));
    double complex *Hf = (double complex *)calloc(N, sizeof(double complex));
    double complex *Yf = (double complex *)calloc(N, sizeof(double complex));

    if (!block_buf || !y_full || !Xf || !Hf || !Yf) {
        free(block_buf); free(y_full); free(Xf); free(Hf); free(Yf);
        return -1;
    }

    /* Pre-compute Hf */
    double *h_pad = (double *)calloc(N, sizeof(double));
    if (!h_pad) {
        free(block_buf); free(y_full); free(Xf); free(Hf); free(Yf);
        return -1;
    }
    memcpy(h_pad, h->samples, M * sizeof(double));
    dft_direct(h_pad, (int)N, Hf);
    free(h_pad);

    size_t y_idx = 0;
    size_t x_idx = 0;
    int first_block = 1;

    while (x_idx < Nx + M - 1) {
        memset(block_buf, 0, N * sizeof(double));

        if (first_block) {
            /* First block: prepend M-1 zeros, then N-M+1 samples of x */
            size_t take = (valid_per_block < Nx) ? valid_per_block : Nx;
            memcpy(block_buf + M - 1, x->samples, take * sizeof(double));
            x_idx = take;
            first_block = 0;
        } else {
            /* Subsequent blocks: last M-1 from previous + new samples */
            size_t take = (x_idx + valid_per_block <= Nx)
                          ? valid_per_block : (Nx - x_idx);
            /* Shift: copy last M-1 of previous block (handled by circular buffer logic) */
            /* For simplicity, just take next chunk starting from x_idx */
            if (x_idx >= M - 1) {
                memcpy(block_buf, x->samples + x_idx - (M - 1),
                       (take + M - 1) * sizeof(double));
            } else {
                memcpy(block_buf + M - 1 - x_idx, x->samples,
                       (take + x_idx) * sizeof(double));
            }
            x_idx += take;
        }

        dft_direct(block_buf, (int)N, Xf);
        for (size_t k = 0; k < N; k++) {
            Yf[k] = Xf[k] * Hf[k];
        }
        idft_direct(Yf, (int)N, y_full);

        /* Save valid portion (discard first M-1 samples) */
        for (size_t i = M - 1; i < N && y_idx < Ny; i++) {
            y->samples[y_idx++] = y_full[i];
        }
    }

    free(block_buf); free(y_full); free(Xf); free(Hf); free(Yf);
    return 0;
}

/* ========================================================================
 *  L2-I09: Deconvolution — Frequency-Domain Division
 *
 *  Given y = x * h and h, estimate x via X(f) = Y(f) / H(f).
 *  Regularization epsilon prevents division by zero or noise amplification.
 *
 *  Wiener deconvolution would include SNR-based regularization.
 * ======================================================================== */

int dt_deconvolution_freq(const dt_signal_t *y, const dt_signal_t *h,
                          double epsilon, dt_signal_t *x_est)
{
    if (!y || !h || !x_est || epsilon < 0.0) return -1;

    size_t Ny = y->num_samples;
    size_t Nh = h->num_samples;

    /* Estimated x length */
    size_t Nx_est = (Ny >= Nh) ? (Ny - Nh + 1) : Ny;
    if (Nx_est == 0) Nx_est = 1;

    size_t Nfft = (size_t)next_power_of_two((int)Ny);
    if (Nfft > SIG_MAX_SAMPLES) return -1;

    double complex *Yf = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *Hf = (double complex *)calloc(Nfft, sizeof(double complex));
    double complex *Xf = (double complex *)calloc(Nfft, sizeof(double complex));

    if (!Yf || !Hf || !Xf) {
        free(Yf); free(Hf); free(Xf);
        return -1;
    }

    double *y_pad = (double *)calloc(Nfft, sizeof(double));
    double *h_pad = (double *)calloc(Nfft, sizeof(double));
    if (!y_pad || !h_pad) {
        free(Yf); free(Hf); free(Xf); free(y_pad); free(h_pad);
        return -1;
    }

    memcpy(y_pad, y->samples, Ny * sizeof(double));
    memcpy(h_pad, h->samples, Nh * sizeof(double));

    dft_direct(y_pad, (int)Nfft, Yf);
    dft_direct(h_pad, (int)Nfft, Hf);

    /* Division with regularization */
    for (size_t k = 0; k < Nfft; k++) {
        double mag_sq = creal(Hf[k]) * creal(Hf[k]) +
                        cimag(Hf[k]) * cimag(Hf[k]);
        if (mag_sq < epsilon * epsilon) {
            Xf[k] = 0.0 + 0.0 * I;  /* Suppress near-zero frequencies */
        } else {
            Xf[k] = Yf[k] / Hf[k];
        }
    }

    double *x_full = (double *)calloc(Nfft, sizeof(double));
    if (!x_full) {
        free(Yf); free(Hf); free(Xf); free(y_pad); free(h_pad);
        return -1;
    }

    idft_direct(Xf, (int)Nfft, x_full);

    if (x_est->num_samples < Nx_est) {
        free(Yf); free(Hf); free(Xf); free(y_pad); free(h_pad); free(x_full);
        return -1;
    }
    memcpy(x_est->samples, x_full, Nx_est * sizeof(double));

    free(Yf); free(Hf); free(Xf); free(y_pad); free(h_pad); free(x_full);
    return 0;
}

/* ========================================================================
 *  L2-I10: Cross-Correlation R_{xy}[m] = sum x[n] * y[n+m]
 *
 *  Related to convolution: R_{xy}[m] = x[-m] * y[m]
 *  Output length: Nx + Ny - 1, centered around lag = 0
 * ======================================================================== */

int dt_cross_correlation(const dt_signal_t *x, const dt_signal_t *y,
                         dt_signal_t *r_xy)
{
    if (!x || !y || !r_xy) return -1;
    if (!x->samples || !y->samples) return -1;

    size_t Nx = x->num_samples;
    size_t Ny = y->num_samples;
    size_t Nr = Nx + Ny - 1;

    if (r_xy->num_samples < Nr) return -1;
    memset(r_xy->samples, 0, r_xy->num_samples * sizeof(double));

    /* lag m ranges from -(Nx-1) to +(Ny-1) */
    /* Store as r_xy[k] where k = m + (Nx-1), so k in [0, Nx+Ny-2] */
    for (size_t k = 0; k < Nr; k++) {
        int m = (int)k - (int)(Nx - 1);
        double sum = 0.0;

        for (size_t n = 0; n < Nx; n++) {
            int idx = (int)n + m;
            if (idx >= 0 && idx < (int)Ny) {
                sum += x->samples[n] * y->samples[idx];
            }
        }
        r_xy->samples[k] = sum;
    }
    return 0;
}

/* ========================================================================
 *  L2-I11: Auto-Correlation R_{xx}[m] = sum x[n] * x[n+m]
 *  Symmetric: R_{xx}[-m] = R_{xx}[m]
 * ======================================================================== */

int dt_auto_correlation(const dt_signal_t *x, dt_signal_t *r_xx)
{
    return dt_cross_correlation(x, x, r_xx);
}

/* ========================================================================
 *  L2-I12: Normalized Cross-Correlation (Correlation Coefficient)
 *
 *  rho_{xy}[m] = R_{xy}[m] / sqrt(R_{xx}[0] * R_{yy}[0])
 *  Bounded: |rho| <= 1 by Cauchy-Schwarz inequality
 * ======================================================================== */

int dt_cross_correlation_normalized(const dt_signal_t *x, const dt_signal_t *y,
                                     dt_signal_t *rho_xy)
{
    if (!x || !y || !rho_xy) return -1;

    /* Compute Rxy */
    int ret = dt_cross_correlation(x, y, rho_xy);
    if (ret != 0) return ret;

    /* Compute Rxx[0] and Ryy[0] (energies) */
    double Ex = dt_signal_energy(x);
    double Ey = dt_signal_energy(y);

    if (Ex < 1e-30 || Ey < 1e-30) {
        /* Zero signal -> all zeros */
        memset(rho_xy->samples, 0, rho_xy->num_samples * sizeof(double));
        return 0;
    }

    double norm_factor = sqrt(Ex * Ey);
    for (size_t i = 0; i < rho_xy->num_samples; i++) {
        rho_xy->samples[i] /= norm_factor;
    }
    return 0;
}

/* ========================================================================
 *  L2-I13: Step Response from Impulse Response
 *
 *  CT: s(t) = integral_0^t h(tau) dtau
 *  DT: s[n] = sum_{k=0}^{n} h[k]  (cumulative sum)
 * ======================================================================== */

int dt_step_from_impulse(const dt_signal_t *h, dt_signal_t *s)
{
    if (!h || !s) return -1;
    if (s->num_samples < h->num_samples) return -1;

    double accum = 0.0;
    for (size_t n = 0; n < h->num_samples; n++) {
        accum += h->samples[n];
        s->samples[n] = accum;
    }
    return 0;
}

/* ========================================================================
 *  L2-I14: Impulse Response from Step Response
 *
 *  CT: h(t) = d[s(t)]/dt
 *  DT: h[0] = s[0]; h[n] = s[n] - s[n-1] for n >= 1
 * ======================================================================== */

int dt_impulse_from_step(const dt_signal_t *s, dt_signal_t *h)
{
    if (!s || !h) return -1;
    if (h->num_samples < s->num_samples) return -1;
    if (s->num_samples == 0) return -1;

    h->samples[0] = s->samples[0];
    for (size_t n = 1; n < s->num_samples; n++) {
        h->samples[n] = s->samples[n] - s->samples[n - 1];
    }
    return 0;
}

/* ========================================================================
 *  L2-I15: Convolution Commutativity Check
 *
 *  Property: x * h = h * x
 *  Returns max absolute difference between the two convolution orders.
 * ======================================================================== */

double dt_convolution_commutativity_check(const dt_signal_t *x,
                                           const dt_signal_t *h)
{
    if (!x || !h) return -1.0;

    size_t N = x->num_samples + h->num_samples - 1;
    dt_signal_t y1 = dt_signal_alloc(N);
    dt_signal_t y2 = dt_signal_alloc(N);

    if (!y1.samples || !y2.samples) {
        dt_signal_free(&y1);
        dt_signal_free(&y2);
        return -1.0;
    }

    dt_convolution_direct(x, h, &y1);
    dt_convolution_direct(h, x, &y2);

    double max_diff = dt_signal_max_abs_diff(&y1, &y2);

    dt_signal_free(&y1);
    dt_signal_free(&y2);
    return max_diff;
}

/* ========================================================================
 *  L2-I16: Zero-Padding Utility
 * ======================================================================== */

int dt_signal_zero_pad(const dt_signal_t *src, size_t new_len,
                       dt_signal_t *dst)
{
    if (!src || !dst || new_len < src->num_samples) return -1;
    if (dst->num_samples < new_len) return -1;

    memcpy(dst->samples, src->samples, src->num_samples * sizeof(double));
    memset(dst->samples + src->num_samples, 0,
           (new_len - src->num_samples) * sizeof(double));
    return 0;
}
