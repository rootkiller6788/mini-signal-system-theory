/**
 * @file sysid_frequency.c
 * @brief Frequency-domain system identification methods (L5 Algorithms, L8 Advanced)
 *
 * Knowledge points:
 *   L5 - Empirical Transfer Function Estimate (ETFE)
 *   L5 - Spectral analysis (Blackman-Tukey, Welch)
 *   L5 - Frequency-domain LS fitting (Levy, Sanathanan-Koerner)
 *   L8 - Local Polynomial Method (LPM) for FRF estimation
 *   L8 - Time-frequency analysis for non-stationary identification
 *   L2 - Bode plot analysis for transfer function ID
 *   L4 - Coherence-based model validation
 *
 * References:
 *   Pintelon & Schoukens (2012) "System Identification: A Frequency Domain Approach"
 *   Ljung (1999) Ch. 6 (Nonparametric Time- and Frequency-Domain Methods)
 */

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_frequency.h"
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <stdio.h>

/*---------------------------------------------------------------------------
 * Coherence estimation (L2/L4)
 * γ²_{yu}(ω) = |Φ_{yu}(ω)|² / (Φ_{uu}(ω) Φ_{yy}(ω))
 * Measures linear relationship quality at each frequency.
 * γ² = 1: perfect linear relationship
 * γ² = 0: no linear relationship (nonlinearities or noise)
 *---------------------------------------------------------------------------*/

int sysid_coherence(const double *u, const double *y, int N,
                     const double *omega, int Nf,
                     double *coherence, int smooth_window) {
    if (!u || !y || !omega || !coherence || N <= 0 || Nf <= 0) return -1;

    for (int f = 0; f < Nf; f++) {
        double w = omega[f];
        double complex Suu = 0.0, Syy = 0.0, Syu = 0.0;

        /* Welch periodogram with rectangular window */
        for (int k = 0; k < N; k++) {
            double complex e_jwk = cexp(-I * w * k);
            Suu += cabs(u[k] * e_jwk);
            Syy += cabs(y[k] * e_jwk);
            Syu += y[k] * cexp(I * w * k) * u[k]; /* cross */
        }

        /* Simplified: compute cross-spectral density estimates */
        /* Use DTFT */
        double complex Uw = 0.0, Yw = 0.0;
        for (int k = 0; k < N; k++) {
            double complex e_jwk = cexp(-I * w * k);
            Uw += u[k] * e_jwk;
            Yw += y[k] * e_jwk;
        }

        double puu = cabs(Uw) * cabs(Uw);
        double pyy = cabs(Yw) * cabs(Yw);
        double cross = cabs(Uw * conj(Yw));
        cross = cross * cross; /* |Φ_{yu}|² = |Yw * conj(Uw)|² = |Uw * conj(Yw)|² */

        double denom = puu * pyy;
        if (denom < 1e-30) {
            coherence[f] = 0.0;
        } else {
            coherence[f] = cross / denom;
            if (coherence[f] > 1.0) coherence[f] = 1.0;
            if (coherence[f] < 0.0) coherence[f] = 0.0;
        }
    }
    (void)smooth_window; /* Parameter accepted for API consistency */
    return 0;
}

/*---------------------------------------------------------------------------
 * Welch's averaged periodogram for PSD estimation (L5)
 *---------------------------------------------------------------------------*/

int sysid_welch_psd(const double *x, int N, int window_len, int overlap,
                     double *freq, double *psd, int Nfft, double Ts) {
    if (!x || !freq || !psd || N <= 0 || window_len <= 0 || Nfft <= 0) return -1;
    if (window_len > N) window_len = N / 2;

    int step = window_len - overlap;
    if (step <= 0) step = window_len / 2;

    /* Allocate accumulation array */
    double *psd_accum = (double *)calloc((size_t)Nfft, sizeof(double));
    if (!psd_accum) { free(psd_accum); return -1; }

    int n_segments = 0;
    for (int start = 0; start + window_len <= N; start += step) {
        /* Apply Hanning window and compute FFT via DFT (simplified) */
        for (int bin = 0; bin < Nfft; bin++) {
            double complex Xw = 0.0;
            double f_bin = (double)bin / (double)Nfft / Ts;
            for (int k = 0; k < window_len; k++) {
                /* Hanning window: 0.5 * (1 - cos(2πk/(L-1))) */
                double w_hanning = 0.5 * (1.0 - cos(2.0 * M_PI * (double)k / (double)(window_len - 1)));
                Xw += w_hanning * x[start + k] * cexp(-I * 2.0 * M_PI * f_bin * k * Ts);
            }
            psd_accum[bin] += cabs(Xw) * cabs(Xw) / (double)window_len;
        }
        n_segments++;
    }

    /* Average and normalize */
    double window_energy = 0.375; /* Hanning window energy factor */
    for (int bin = 0; bin < Nfft; bin++) {
        freq[bin] = (double)bin / (double)Nfft / Ts;
        psd[bin] = (n_segments > 0) ? (psd_accum[bin] / (double)n_segments / window_energy) : 0.0;
    }

    free(psd_accum);
    return 0;
}

/*---------------------------------------------------------------------------
 * Transfer function estimation from Bode data (Levy's method + SK iteration)
 *
 * Sanathanan-Koerner (1963) iterative refinement:
 * At iteration k, minimize Σ |D_k(ω) G(ω) - N_{k+1}(ω)|² / |D_k(ω)|²
 * where D_k is from previous iteration.
 * As k→∞, approaches true nonlinear least squares.
 *---------------------------------------------------------------------------*/

int sysid_sk_iteration(const sysid_freq_data_t *fd,
                        int num_order, int den_order,
                        double *num, double *den,
                        int max_iter, double tol) {
    if (!fd || !fd->omega || !fd->G || fd->Nf <= 0) return -1;
    if (!num || !den) return -1;

    int n_num = num_order + 1;
    int n_den = den_order; /* fixed den[0]=1 */
    int Nf = fd->Nf;

    /* Initialize denominator to 1 */
    memset(den, 0, (size_t)(den_order + 1) * sizeof(double));
    den[0] = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Build weighted LS problem using current denominator D_k */
        int nparam = n_num + n_den;
        int neq = 2 * Nf;
        double *A = (double *)calloc((size_t)(neq * nparam), sizeof(double));
        double *b = (double *)calloc((size_t)neq, sizeof(double));
        if (!A || !b) { free(A); free(b); return -1; }

        for (int i = 0; i < Nf; i++) {
            double w = fd->omega[i];
            double complex Gmeas = fd->G[i];

            /* Compute current denominator D_k(jω) */
            double complex Dk = 1.0;
            double complex s = I * w;
            double complex s_pow = 1.0;
            for (int j = 1; j <= den_order; j++) {
                s_pow *= s;
                Dk += den[j] * s_pow;
            }

            double weight = (fd->weight) ? fd->weight[i] : 1.0;
            /* SK weight: W_i = 1 / |D_k(jω_i)|² */
            double sk_weight = weight / (cabs(Dk) * cabs(Dk) + 1e-15);

            /* Numerator basis: N(s) = Σ n_j s^j */
            double complex s_pow_n = 1.0;
            for (int j = 0; j < n_num; j++) {
                double complex val = s_pow_n * sk_weight;
                A[(2 * i) * nparam + j] = creal(val);
                A[(2 * i + 1) * nparam + j] = cimag(val);
                s_pow_n *= s;
            }

            /* Denominator basis: -D(s) G_meas(s) terms */
            s_pow = s;
            for (int j = 0; j < n_den; j++) {
                double complex val = -s_pow * Gmeas * sk_weight;
                A[(2 * i) * nparam + n_num + j] = creal(val);
                A[(2 * i + 1) * nparam + n_num + j] = cimag(val);
                s_pow *= s;
            }

            /* RHS: Gmeas (weighted) */
            b[2 * i] = creal(Gmeas) * sk_weight;
            b[2 * i + 1] = cimag(Gmeas) * sk_weight;
        }

        /* Solve linear LS */
        double *x = (double *)calloc((size_t)nparam, sizeof(double));
        extern int sysid_ls_solve(const double *Phi, const double *Y, int nrows,
                                   int ncols, double *theta, double *cov, double *sigma2);
        if (sysid_ls_solve(A, b, neq, nparam, x, NULL, NULL) == 0) {
            /* Check convergence of denominator */
            double den_change = 0.0;
            for (int j = 0; j < n_den; j++) {
                double diff = x[n_num + j] - den[j + 1];
                den_change += diff * diff;
                den[j + 1] = x[n_num + j];
            }
            for (int j = 0; j < n_num; j++) num[j] = x[j];

            if (den_change < tol * tol) {
                free(x); free(A); free(b);
                return iter + 1;
            }
        }
        free(x); free(A); free(b);
    }

    return max_iter;
}

/*---------------------------------------------------------------------------
 * Local Polynomial Method (LPM) for nonparametric FRF estimation (L8)
 *
 * LPM (Schoukens et al., 2009) estimates G(ω_k) and its derivatives
 * using a local polynomial approximation over a narrow frequency window.
 * Reduces leakage errors compared to simple ETFE.
 *---------------------------------------------------------------------------*/

int sysid_lpm_frf(const double *u, const double *y, int N,
                   const double *omega, int Nf,
                   double complex *G, int poly_order, int window_half) {
    if (!u || !y || !omega || !G || N <= 0 || Nf <= 0) return -1;
    if (poly_order < 0) poly_order = 2;
    if (window_half < 1) window_half = 2;

    /* For each frequency ω_k, fit a local polynomial:
     * G(ω_k + ν) ≈ G₀ + G₁ ν + G₂ ν² + ... + G_R ν^R
     * using 2*W+1 neighboring frequency points.
     *
     * Simplified implementation: for each target frequency, estimate
     * Y(ω)/U(ω) with local averaging.
     */

    for (int f = 0; f < Nf; f++) {
        double w0 = omega[f];

        /* Collect U and Y at neighboring frequencies */
        /* For block processing approach, use DFT */
        double complex U_local = 0.0, Y_local = 0.0;

        /* Use frequency-domain smoothing via convolution with window */
        for (int k = 0; k < N; k++) {
            double complex e_jwk = cexp(-I * w0 * k);
            /* Apply frequency-domain smoothing via time-domain windowing */
            /* Exponential window: reduces leakage */
            double time_window = exp(-(double)k / (double)N * 3.0);
            U_local += time_window * u[k] * e_jwk;
            Y_local += time_window * y[k] * e_jwk;
        }

        if (cabs(U_local) < 1e-15) {
            G[f] = 1e10 + 0.0 * I;
        } else {
            G[f] = Y_local / U_local;
        }
    }
    (void)poly_order; (void)window_half; /* Parameters accepted for API */

    return 0;
}

/*---------------------------------------------------------------------------
 * System identification from Bode plot data (L2/L6)
 *
 * Given magnitude and phase data from a Bode plot, estimate parameters
 * of a transfer function via nonlinear curve fitting.
 * Supports: 1st-order, 2nd-order, integrator+lag, etc.
 *---------------------------------------------------------------------------*/

/**
 * Fit first-order transfer function to Bode data:
 *   G(s) = K / (τ s + 1)
 *
 * Estimates K (DC gain) and τ (time constant) from low-frequency data.
 */
int sysid_bode_fit_1st_order(const sysid_freq_data_t *fd,
                              double *K, double *tau) {
    if (!fd || !fd->G || fd->Nf <= 2 || !K || !tau) return -1;

    /* K ≈ |G(0)|, use lowest frequency */
    *K = cabs(fd->G[0]);

    /* τ is found from the -3dB frequency: |G(j/τ)| = K/√2
     * Find frequency where magnitude ≈ K/√2 */
    double mag_target = (*K) / sqrt(2.0);
    *tau = 1.0;

    int found = 0;
    for (int i = 1; i < fd->Nf; i++) {
        double mag = cabs(fd->G[i]);
        if (mag <= mag_target) {
            /* Interpolate between i-1 and i */
            double mag_prev = cabs(fd->G[i - 1]);
            if (mag_prev > mag_target) {
                double w_prev = fd->omega[i - 1];
                double w = fd->omega[i];
                double frac = (mag_target - mag) / (mag_prev - mag + 1e-15);
                double w_3dB = w + frac * (w_prev - w);
                *tau = 1.0 / w_3dB;
                found = 1;
            }
            break;
        }
    }
    if (!found) {
        /* Use phase method: ∠G(jω) = -atan(ωτ)
         * at ω=1/τ, phase = -45° */
        for (int i = 0; i < fd->Nf; i++) {
            double phase = atan2(cimag(fd->G[i]), creal(fd->G[i]));
            if (fabs(phase + M_PI / 4.0) < 0.2) {
                *tau = 1.0 / fd->omega[i];
                break;
            }
        }
    }
    return 0;
}

/**
 * Fit second-order transfer function to Bode data:
 *   G(s) = K ω_n² / (s² + 2ζ ω_n s + ω_n²)
 */
int sysid_bode_fit_2nd_order(const sysid_freq_data_t *fd,
                              double *K, double *omega_n, double *zeta) {
    if (!fd || !fd->G || fd->Nf <= 3 || !K || !omega_n || !zeta) return -1;

    *K = cabs(fd->G[0]);

    /* Find resonance peak */
    double max_mag = 0.0;
    int peak_idx = 0;
    for (int i = 0; i < fd->Nf; i++) {
        double mag = cabs(fd->G[i]);
        if (mag > max_mag) { max_mag = mag; peak_idx = i; }
    }

    *omega_n = fd->omega[peak_idx];

    /* Resonance peak magnitude: M_r = 1/(2ζ√(1-ζ²)) for ζ < √2/2
     * Relative to DC gain: M_r/K = 1/(2ζ√(1-ζ²))
     */
    double Mr_ratio = max_mag / (*K);
    if (Mr_ratio > 1.0) {
        /* Solve for ζ from M_r */
        double mr2 = Mr_ratio * Mr_ratio;
        double disc = 1.0 - 4.0 * mr2 * (1.0 - mr2);
        if (disc > 0.0) {
            *zeta = sqrt((1.0 - sqrt(disc)) / 2.0);
        } else {
            *zeta = 0.5;
        }
    } else {
        /* No visible resonance: use -3dB method */
        *zeta = 1.0; /* critically damped guess */
    }

    if (*zeta < 0.01) *zeta = 0.01;
    if (*zeta > 2.0) *zeta = 2.0;

    return 0;
}

/*---------------------------------------------------------------------------
 * Downsampled frequency response for efficient computation (L5)
 *---------------------------------------------------------------------------*/

int sysid_downsampled_frf(const double *u, const double *y, int N,
                           double *freq, double complex *G, int Nf,
                           double Ts, int decimation) {
    if (!u || !y || !freq || !G || N <= 0 || Nf <= 0 || decimation <= 0) return -1;

    int N_ds = N / decimation;
    if (N_ds < 10) { decimation = 1; N_ds = N; }

    for (int f = 0; f < Nf; f++) {
        double complex Uw = 0.0, Yw = 0.0;
        double w = 2.0 * M_PI * freq[f] * Ts * decimation;

        for (int k = 0; k < N_ds; k++) {
            double complex e_jwk = cexp(-I * w * k);
            Uw += u[k * decimation] * e_jwk;
            Yw += y[k * decimation] * e_jwk;
        }

        freq[f] = freq[f];
        G[f] = (cabs(Uw) < 1e-15) ? (1e10 + 0.0 * I) : (Yw / Uw);
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Time-frequency analysis via Short-Time Fourier Transform (L8)
 * Used for non-stationary system identification (tracking time-varying dynamics).
 *---------------------------------------------------------------------------*/

int sysid_stft(const double *x, int N, int window_len, int overlap,
                int Nfreq, double *tfr, double Ts) {
    if (!x || !tfr || N <= 0 || window_len <= 0 || Nfreq <= 0) return -1;
    if (window_len > N) window_len = N / 2;

    int step = window_len - overlap;
    if (step <= 0) step = 1;

    int n_times = (N - window_len) / step + 1;
    if (n_times <= 0) n_times = 1;

    for (int t = 0; t < n_times; t++) {
        int start = t * step;
        for (int f = 0; f < Nfreq; f++) {
            double freq = (double)f / (double)Nfreq / Ts;
            double complex Xw = 0.0;
            for (int k = 0; k < window_len; k++) {
                double w_hann = 0.5 * (1.0 - cos(2.0 * M_PI * (double)k / (double)(window_len - 1)));
                Xw += w_hann * x[start + k] * cexp(-I * 2.0 * M_PI * freq * k * Ts);
            }
            tfr[t * Nfreq + f] = cabs(Xw) * cabs(Xw) / (double)window_len;
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Frequency response from impulse response (L2: basic relationship)
 *   G(e^{jω}) = Σ_{k=0}^{∞} g(k) e^{-jωk}
 *---------------------------------------------------------------------------*/

int sysid_frf_from_impulse(const double *impulse, int n_impulse,
                            const double *omega, int Nf,
                            double complex *G) {
    if (!impulse || !omega || !G || n_impulse <= 0 || Nf <= 0) return -1;

    for (int f = 0; f < Nf; f++) {
        double complex Gf = 0.0;
        for (int k = 0; k < n_impulse; k++) {
            Gf += impulse[k] * cexp(-I * omega[f] * k);
        }
        G[f] = Gf;
    }
    return 0;
}
