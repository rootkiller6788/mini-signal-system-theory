/**
 * @file sysid_signals.c
 * @brief Excitation signal generation for system identification (L2, L4, L6)
 */

#include "sysid_signals.h"
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>

/* Simple LCG random number generator: state = (a * state + c) mod m */
static unsigned int lcg_rand(unsigned int *state) {
    *state = (1103515245U * (*state) + 12345U) & 0x7fffffffU;
    return *state;
}

/* Uniform random in [0, 1) */
static double lcg_uniform(unsigned int *state) {
    return (double)lcg_rand(state) / (double)0x80000000U;
}

/*---------------------------------------------------------------------------
 * LFSR / PRBS (L6)
 *---------------------------------------------------------------------------*/

uint32_t sysid_lfsr_step(uint32_t state, const int *taps, int n_taps) {
    uint32_t feedback = 0;
    for (int i = 0; i < n_taps; i++) {
        feedback ^= (state >> (taps[i] - 1)) & 1U;
    }
    state = (state << 1) | feedback;
    /* Mask to keep within bit width (use highest tap to determine mask) */
    /* Find highest tap for mask */
    int max_tap = 0;
    for (int i = 0; i < n_taps; i++) {
        if (taps[i] > max_tap) max_tap = taps[i];
    }
    uint32_t mask = (1U << max_tap) - 1;
    if (max_tap >= 32) mask = 0xFFFFFFFFU;
    state &= mask;
    if (state == 0) state = 1; /* avoid lock-up */
    return state;
}

int sysid_prbs_init(sysid_prbs_t *prbs, int reg_length, const int *taps,
                     int n_taps, int amplitude) {
    if (!prbs || reg_length <= 0 || reg_length > 31 || !taps || n_taps <= 0) return -1;
    prbs->taps = (int *)malloc((size_t)n_taps * sizeof(int));
    if (!prbs->taps) return -1;
    memcpy(prbs->taps, taps, (size_t)n_taps * sizeof(int));
    prbs->n_taps = n_taps;
    prbs->reg_length = reg_length;
    prbs->state = 1; /* initial state, non-zero */
    prbs->amplitude = amplitude;
    return 0;
}

void sysid_prbs_free(sysid_prbs_t *prbs) {
    if (!prbs) return;
    free(prbs->taps);
    prbs->taps = NULL;
}

int sysid_prbs_next(sysid_prbs_t *prbs) {
    if (!prbs) return 0;
    int bit = (int)(prbs->state & 1U);
    prbs->state = sysid_lfsr_step(prbs->state, prbs->taps, prbs->n_taps);
    return (bit) ? prbs->amplitude : -prbs->amplitude;
}

int sysid_prbs_generate(sysid_prbs_t *prbs, double *signal, int N, int T_clk) {
    if (!prbs || !signal || N <= 0 || T_clk <= 0) return -1;
    for (int k = 0; k < N; k++) {
        if (k % T_clk == 0) {
            /* Advance register only at clock rate */
            int val = sysid_prbs_next(prbs);
            /* Fill T_clk samples with same value */
            int end = k + T_clk;
            if (end > N) end = N;
            for (int j = k; j < end; j++) {
                signal[j] = (double)val;
            }
        }
        /* else: already filled by the loop above */
    }
    /* Re-fill to ensure all set */
    int cur_val = 0;
    for (int k = 0; k < N; k++) {
        if (k % T_clk == 0) {
            cur_val = (int)(prbs->state & 1U);
            prbs->state = sysid_lfsr_step(prbs->state, prbs->taps, prbs->n_taps);
        }
        signal[k] = (double)(cur_val ? prbs->amplitude : -prbs->amplitude);
    }
    return 0;
}

int sysid_prbs_generate_std(int reg_length, double *signal, int N,
                             int T_clk, double amplitude) {
    if (reg_length <= 0 || reg_length > 31 || !signal || N <= 0) return -1;

    /* Standard maximal-length tap configurations */
    static const int taps_table[32][4] = {
        {0}, {0}, {2,1}, {3,2}, {4,3}, {5,3}, {6,5}, {7,6},
        {8,6,5,4}, {9,5}, {10,7}, {11,9}, {12,11,8,6},
        {13,12,10,9}, {14,13,11,9}, {15,14}, {16,14,13,11},
        {17,14}, {18,11}, {19,18,17,14}, {20,17},
        {21,19}, {22,21}, {23,18}, {24,23,21,20},
        {25,22}, {26,25,24,20}, {27,26,25,22},
        {28,25}, {29,27}, {30,29,26,24}, {31,28}
    };

    int n_taps = 0;
    for (int i = 0; i < 4 && taps_table[reg_length][i] != 0; i++) n_taps++;
    if (n_taps == 0) return -1;

    sysid_prbs_t prbs;
    if (sysid_prbs_init(&prbs, reg_length, taps_table[reg_length], n_taps, (int)amplitude) != 0)
        return -1;

    /* Set initial state to non-zero */
    prbs.state = 0x5A5AU & ((1U << reg_length) - 1);
    if (prbs.state == 0) prbs.state = 1;

    for (int k = 0; k < N; k++) {
        if (k % T_clk == 0) {
            int bit = (int)(prbs.state & 1U);
            prbs.state = sysid_lfsr_step(prbs.state, prbs.taps, prbs.n_taps);
            int val = bit ? (int)amplitude : -(int)amplitude;
            int end = k + T_clk;
            if (end > N) end = N;
            for (int j = k; j < end; j++) signal[j] = (double)val;
        }
    }

    sysid_prbs_free(&prbs);
    return 0;
}

/*---------------------------------------------------------------------------
 * RBS (L6)
 *---------------------------------------------------------------------------*/

int sysid_rbs_generate(double *signal, int N, double amplitude, unsigned int seed) {
    if (!signal || N <= 0) return -1;
    unsigned int state = seed;
    for (int i = 0; i < N; i++) {
        double u = lcg_uniform(&state);
        signal[i] = (u < 0.5) ? amplitude : -amplitude;
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Multisine (L6)
 *---------------------------------------------------------------------------*/

int sysid_multisine_alloc(sysid_multisine_t *ms, int Nf, double Ts) {
    if (!ms || Nf <= 0) return -1;
    memset(ms, 0, sizeof(sysid_multisine_t));
    ms->frequencies = (double *)calloc((size_t)Nf, sizeof(double));
    ms->amplitudes   = (double *)calloc((size_t)Nf, sizeof(double));
    ms->phases      = (double *)calloc((size_t)Nf, sizeof(double));
    if (!ms->frequencies || !ms->amplitudes || !ms->phases) {
        free(ms->frequencies); free(ms->amplitudes); free(ms->phases);
        return -1;
    }
    ms->Nf = Nf;
    ms->Ts = Ts;
    ms->offset = 0.0;
    return 0;
}

void sysid_multisine_free(sysid_multisine_t *ms) {
    if (!ms) return;
    free(ms->frequencies); free(ms->amplitudes); free(ms->phases);
    ms->frequencies = NULL; ms->amplitudes = NULL; ms->phases = NULL;
    ms->Nf = 0;
}

void sysid_multisine_set_flat_amp(sysid_multisine_t *ms, double amplitude) {
    if (!ms || !ms->amplitudes) return;
    for (int i = 0; i < ms->Nf; i++) {
        ms->amplitudes[i] = amplitude / sqrt((double)ms->Nf);
    }
}

void sysid_multisine_set_schroeder_phase(sysid_multisine_t *ms, unsigned int seed) {
    if (!ms || !ms->phases || ms->Nf <= 0) return;
    unsigned int state = seed;
    /* Schroeder's rule: φ_k = φ_0 - 2π Σ_{j=0}^{k-1} (k-j) P_j
     * where P_j = A_j² / Σ A_i² are normalized powers.
     */
    double total_power = 0.0;
    for (int i = 0; i < ms->Nf; i++) {
        total_power += ms->amplitudes[i] * ms->amplitudes[i];
    }
    if (total_power < 1e-15) total_power = 1.0;

    /* φ_0 random */
    double phase = lcg_uniform(&state) * 2.0 * M_PI;
    ms->phases[0] = phase;

    double cum_power = ms->amplitudes[0] * ms->amplitudes[0] / total_power;
    for (int k = 1; k < ms->Nf; k++) {
        double power_k = ms->amplitudes[k] * ms->amplitudes[k] / total_power;
        /* Simplified: φ_k = -π k(k-1)/Nf (Newman phase, good for flat spectra) */
        /* Schroeder: φ_k = φ_{k-1} - 2π * (k-1) * P_{k-1} */
        phase -= 2.0 * M_PI * (double)(k) * cum_power / (double)ms->Nf;
        ms->phases[k] = phase;
        cum_power += power_k;
    }
}

int sysid_multisine_generate(const sysid_multisine_t *ms, double *signal, int N) {
    if (!ms || !signal || N <= 0) return -1;
    for (int k = 0; k < N; k++) {
        double val = ms->offset;
        double t = k * ms->Ts;
        for (int f = 0; f < ms->Nf; f++) {
            val += ms->amplitudes[f] * sin(2.0 * M_PI * ms->frequencies[f] * t + ms->phases[f]);
        }
        signal[k] = val;
    }
    return 0;
}

double sysid_crest_factor(const double *signal, int N) {
    if (!signal || N <= 0) return 0.0;
    double max_abs = 0.0, rms = 0.0;
    for (int i = 0; i < N; i++) {
        double abs_val = fabs(signal[i]);
        if (abs_val > max_abs) max_abs = abs_val;
        rms += signal[i] * signal[i];
    }
    rms = sqrt(rms / (double)N);
    if (rms < 1e-15) return 1e10;
    return max_abs / rms;
}

int sysid_multisine_logspace(int Nf, double f_min, double f_max, double Ts,
                              double amplitude, double *signal, int N) {
    if (Nf <= 0 || f_min <= 0 || f_max <= f_min || !signal || N <= 0) return -1;

    sysid_multisine_t ms;
    if (sysid_multisine_alloc(&ms, Nf, Ts) != 0) return -1;

    /* Log-spaced frequencies */
    double log_fmin = log(f_min);
    double log_fmax = log(f_max);
    for (int i = 0; i < Nf; i++) {
        ms.frequencies[i] = exp(log_fmin + (log_fmax - log_fmin) * (double)i / (double)(Nf - 1));
    }

    sysid_multisine_set_flat_amp(&ms, amplitude);
    sysid_multisine_set_schroeder_phase(&ms, 12345U);
    int ret = sysid_multisine_generate(&ms, signal, N);
    sysid_multisine_free(&ms);
    return ret;
}

/*---------------------------------------------------------------------------
 * Chirp signals (L6)
 *---------------------------------------------------------------------------*/

int sysid_chirp_linear(double *signal, int N, double f0, double f1,
                        double Ts, double amplitude) {
    if (!signal || N <= 0) return -1;
    for (int k = 0; k < N; k++) {
        double t = k * Ts;
        double T_total = (N - 1) * Ts;
        /* Instantaneous frequency: f(t) = f0 + (f1 - f0) * t / T
         * Phase: φ(t) = 2π ∫ f(τ) dτ = 2π [f0 t + (f1-f0) t²/(2T)] */
        double phase = 2.0 * M_PI * (f0 * t + (f1 - f0) * t * t / (2.0 * T_total));
        signal[k] = amplitude * sin(phase);
    }
    return 0;
}

int sysid_chirp_log(double *signal, int N, double f0, double f1,
                     double Ts, double amplitude) {
    if (!signal || N <= 0 || f0 <= 0 || f1 <= f0) return -1;
    double T_total = (N - 1) * Ts;
    double k_rate = (log(f1) - log(f0)) / T_total;

    for (int k = 0; k < N; k++) {
        double t = k * Ts;
        /* f(t) = f0 * exp(k_rate * t)
         * φ(t) = 2π f0 (exp(k_rate t) - 1) / k_rate */
        if (fabs(k_rate) < 1e-15) {
            signal[k] = amplitude * sin(2.0 * M_PI * f0 * t);
        } else {
            double phase = 2.0 * M_PI * f0 * (exp(k_rate * t) - 1.0) / k_rate;
            signal[k] = amplitude * sin(phase);
        }
    }
    return 0;
}

int sysid_swept_sine(double *signal, int N, const double *inst_freq,
                      double Ts, double amplitude) {
    if (!signal || !inst_freq || N <= 0) return -1;
    double phase = 0.0;
    for (int k = 0; k < N; k++) {
        signal[k] = amplitude * sin(phase);
        phase += 2.0 * M_PI * inst_freq[k] * Ts;
    }
    /* Wrap phase */
    return 0;
}

/*---------------------------------------------------------------------------
 * Step, impulse, noise (L6)
 *---------------------------------------------------------------------------*/

int sysid_step_signal(double *signal, int N, int step_sample, double amplitude) {
    if (!signal || N <= 0) return -1;
    for (int i = 0; i < N; i++) {
        signal[i] = (i >= step_sample) ? amplitude : 0.0;
    }
    return 0;
}

int sysid_impulse_signal(double *signal, int N, double amplitude) {
    if (!signal || N <= 0) return -1;
    signal[0] = amplitude;
    for (int i = 1; i < N; i++) signal[i] = 0.0;
    return 0;
}

int sysid_gaussian_noise(double *signal, int N, double sigma, unsigned int seed) {
    if (!signal || N <= 0) return -1;
    unsigned int state = seed;
    /* Box-Muller transform */
    for (int i = 0; i < N; i += 2) {
        double u1 = lcg_uniform(&state);
        double u2 = lcg_uniform(&state);
        /* Avoid log(0) */
        if (u1 < 1e-10) u1 = 1e-10;
        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;
        signal[i] = sigma * r * cos(theta);
        if (i + 1 < N) signal[i + 1] = sigma * r * sin(theta);
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Persistence of excitation analysis (L4)
 *---------------------------------------------------------------------------*/

int sysid_test_pe_order(const double *u, int N, int order, double *R, double tol) {
    if (!u || N <= order || order <= 0 || !R) return -1;

    /* Compute covariance matrix R(τ) for τ = 0..order-1 */
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            double cov = 0.0;
            int start = (i > j) ? i : j;
            int cnt = 0;
            for (int k = start; k < N; k++) {
                cov += u[k - i] * u[k - j];
                cnt++;
            }
            R[i * order + j] = (cnt > 0) ? cov / (double)cnt : 0.0;
        }
    }

    /* Check positive definiteness via smallest eigenvalue estimate */
    /* Use Gershgorin circle theorem: all eigenvalues ≥ mini(R_ii - Σ_{j≠i}|R_ij|) */
    double min_gershgorin = 1e30;
    for (int i = 0; i < order; i++) {
        double radius = 0.0;
        for (int j = 0; j < order; j++) {
            if (j != i) radius += fabs(R[i * order + j]);
        }
        double center = R[i * order + i];
        double lower_bound = center - radius;
        if (lower_bound < min_gershgorin) min_gershgorin = lower_bound;
    }

    return (min_gershgorin > tol) ? 1 : 0;
}

double sysid_pe_quality(const double *u, int N, int order) {
    if (!u || N <= order || order <= 0) return 0.0;
    double *R = (double *)calloc((size_t)(order * order), sizeof(double));
    if (!R) return 0.0;

    sysid_test_pe_order(u, N, order, R, 0.0);

    /* Compute condition number estimate: max/min singular value via Gershgorin */
    double min_bound = 1e30, max_bound = 0.0;
    for (int i = 0; i < order; i++) {
        double radius = 0.0;
        for (int j = 0; j < order; j++) {
            if (j != i) radius += fabs(R[i * order + j]);
        }
        double center = R[i * order + i];
        double low = center - radius;
        double high = center + radius;
        if (low < min_bound) min_bound = low;
        if (high > max_bound) max_bound = high;
    }

    free(R);
    if (min_bound <= 0.0 || max_bound <= 0.0) return 0.0;
    return min_bound / max_bound; /* 1/cond → quality measure, 1=best */
}

int sysid_input_psd(const double *u, int N, double *freq, double *psd,
                     int Nfft, double Ts) {
    if (!u || N <= 0 || !freq || !psd || Nfft <= 0) return -1;
    /* Periodogram PSD estimate */
    for (int bin = 0; bin < Nfft; bin++) {
        double complex Xw = 0.0;
        double f_bin = (double)bin / (double)Nfft / Ts;
        freq[bin] = f_bin;
        for (int k = 0; k < N && k < Nfft; k++) {
            Xw += u[k] * cexp(-I * 2.0 * M_PI * f_bin * k * Ts);
        }
        psd[bin] = cabs(Xw) * cabs(Xw) / (double)N;
    }
    return 0;
}
