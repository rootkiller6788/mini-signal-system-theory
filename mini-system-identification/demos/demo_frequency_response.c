/**
 * @file demo_frequency_response.c
 * @brief Frequency-domain system identification demo (L6/L7)
 *
 * Demonstrates ETF estimation, coherence analysis from noisy I/O data,
 * and parametric ARX fitting for a 2nd-order system.
 */

#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_signals.h"
#include "sysid_validation.h"
#include "sysid_frequency.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Generate noisy data from a 2nd-order low-pass filter (discretized via ZOH) */
static void generate_sys_data(double *u, double *y, int N, double Ts,
                               double wn, double zeta, double noise_std) {
    double a_exp = exp(-zeta * wn * Ts);
    double wd = wn * sqrt(1.0 - zeta * zeta);
    double cos_wd = cos(wd * Ts), sin_wd = sin(wd * Ts);
    double a1 = -2.0 * a_exp * cos_wd;
    double a2 = a_exp * a_exp;
    double b1 = 1.0 - a_exp * (cos_wd + (zeta * wn / wd) * sin_wd);
    double b2 = a_exp * a_exp + a_exp * ((zeta * wn / wd) * sin_wd - cos_wd);

    /* PRBS excitation via a simple LFSR */
    unsigned int lfsr = 0xACE1u;
    for (int k = 0; k < N; k++) {
        unsigned int bit = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5)) & 1u;
        lfsr = (lfsr >> 1) | (bit << 15);
        u[k] = bit ? 1.0 : -1.0;
    }

    /* Simulate output with noise */
    y[0] = b1 * u[0] + noise_std * ((double)rand() / RAND_MAX - 0.5);
    if (N > 1)
        y[1] = -a1 * y[0] + b1 * u[1] + b2 * u[0]
               + noise_std * ((double)rand() / RAND_MAX - 0.5);
    for (int k = 2; k < N; k++) {
        y[k] = -a1 * y[k-1] - a2 * y[k-2]
               + b1 * u[k-1] + b2 * u[k-2]
               + noise_std * ((double)rand() / RAND_MAX - 0.5);
    }
}

int main(void) {
    printf("==============================================================\n");
    printf("  Frequency-Domain System Identification Demo\n");
    printf("  ETF + Coherence + Parametric ARX Fit\n");
    printf("==============================================================\n\n");

    const double Ts   = 0.001;
    const int    N    = 4096;
    const double wn   = 100.0;   /* rad/s */
    const double zeta = 0.1;
    const double noise_std = 0.05;

    printf("Test system: 2nd-order LPF, wn=%.1f rad/s, zeta=%.2f\n", wn, zeta);
    printf("Sampling: Ts=%.4f s, N=%d (%.1f s of data)\n\n", Ts, N, N*Ts);

    double *u = (double *)malloc((size_t)N * sizeof(double));
    double *y = (double *)malloc((size_t)N * sizeof(double));
    if (!u || !y) { printf("Memory error\n"); return 1; }

    generate_sys_data(u, y, N, Ts, wn, zeta, noise_std);

    /* --- Coherence spectrum --- */
    int Nf = 30;
    double *freq      = (double *)malloc((size_t)Nf * sizeof(double));
    double *coherence = (double *)malloc((size_t)Nf * sizeof(double));
    double complex *Getfe = (double complex *)malloc((size_t)Nf * sizeof(double complex));
    double complex *Gtrue = (double complex *)malloc((size_t)Nf * sizeof(double complex));

    for (int i = 0; i < Nf; i++)
        freq[i] = (i + 1) * (0.5 / (N * Ts));

    printf("--- Coherence Spectrum (smooth_window=16) ---\n");
    sysid_coherence(u, y, N, freq, Nf, coherence, 16);
    double avg_coh = 0.0;
    for (int i = 0; i < Nf; i++) avg_coh += coherence[i];
    avg_coh /= Nf;
    printf("  Average coherence: %.4f\n", avg_coh);
    printf("  Interpretation: gamma^2 near 1.0 = high SNR frequency\n");
    printf("                  gamma^2 near 0.0 = noise-dominated\n\n");

    /* --- ETFE frequency response --- */
    printf("--- ETF Frequency Response (first 15 points) ---\n");
    sysid_etfe(u, y, N, freq, Nf, Getfe);
    for (int i = 0; i < Nf; i++) {
        double w = 2.0 * M_PI * freq[i];
        double complex s = w * I;
        Gtrue[i] = wn * wn / (s*s + 2.0*zeta*wn*s + wn*wn);
    }
    printf("  %10s  %10s  %10s  %10s\n", "Freq[Hz]", "|G_etfe|dB", "|G_true|dB", "Phase_est");
    int show = (Nf < 15) ? Nf : 15;
    for (int i = 0; i < show; i++) {
        double db_est  = 20.0 * log10(cabs(Getfe[i]) + 1e-30);
        double db_true = 20.0 * log10(cabs(Gtrue[i]) + 1e-30);
        double ph_deg  = 180.0 / M_PI * carg(Getfe[i]);
        printf("  %10.3f  %10.2f  %10.2f  %10.2f\n", freq[i], db_est, db_true, ph_deg);
    }

    /* --- Parametric ARX fit --- */
    printf("\n--- Parametric ARX Model Fit (na=2, nb=2, nk=1) ---\n");
    sysid_arx_t model;
    sysid_arx_alloc(&model, 2, 2, 1, Ts);

    /* Build regression and solve LS */
    int nrows;
    double *Phi = (double *)malloc((size_t)(N * 4) * sizeof(double));
    double *Yvec = (double *)malloc((size_t)N * sizeof(double));
    sysid_build_regression_arx(u, y, N, 2, 2, 1, Phi, Yvec, &nrows);

    sysid_param_vector_t pv = {0};
    sysid_param_alloc(&pv, 4);
    double V_loss;
    sysid_ls_solve(Phi, Yvec, nrows, 4, pv.theta, NULL, &V_loss);
    sysid_arx_set_params(&model, &pv);

    /* Extract discrete poles */
    double complex poles[2];
    if (sysid_arx_poles(&model, poles, 2) == 0) {
        printf("  Discrete poles:\n");
        for (int i = 0; i < 2; i++) {
            printf("    z%d = %.4f %+.4f j  |z|=%.4f", i+1,
                   creal(poles[i]), cimag(poles[i]), cabs(poles[i]));
            if (cabs(poles[i]) < 1.0) printf(" (stable)");
            else printf(" (unstable!)");
            printf("\n");
        }
        printf("  Continuous equivalents (s = ln(z)/Ts):\n");
        for (int i = 0; i < 2; i++) {
            double complex sp = clog(poles[i]) / Ts;
            double w_est  = cabs(sp);
            double z_est  = -creal(sp) / (w_est + 1e-30);
            printf("    s%d = %.2f %+.2f j  -> wn=%.1f, zeta=%.3f\n",
                   i+1, creal(sp), cimag(sp), w_est, z_est);
        }
    }

    /* Compute NRMSE fit */
    double *ypred = (double *)malloc((size_t)N * sizeof(double));
    sysid_arx_predict(&model, u, y, ypred, N);
    double y_mean = 0.0, y_var = 0.0, mse = 0.0;
    for (int k = 0; k < N; k++) y_mean += y[k];
    y_mean /= N;
    for (int k = 0; k < N; k++) { double dy = y[k] - y_mean; y_var += dy * dy; }
    for (int k = 0; k < N; k++) { double e = y[k] - ypred[k]; mse += e * e; }
    mse /= N;
    double nrmse = (y_var > 1e-30) ? (1.0 - sqrt(mse / y_var)) * 100.0 : 0.0;
    printf("  NRMSE Fit: %.2f %%\n", nrmse);

    /* --- Cross-validation (5-fold) --- */
    printf("\n--- 5-Fold Cross-Validation ---\n");
    double cv_avg_fit, cv_std_fit;
    if (sysid_cross_validate_arx(u, y, N, 2, 2, 1, 5, &cv_avg_fit, &cv_std_fit) == 0) {
        printf("  Cross-validated NRMSE: %.2f +/- %.2f %%\n", cv_avg_fit, cv_std_fit);
    }

    /* Cleanup */
    sysid_arx_free(&model);
    sysid_param_free(&pv);
    free(u); free(y); free(freq); free(coherence);
    free(Getfe); free(Gtrue); free(ypred); free(Phi); free(Yvec);

    printf("\n===== Demo Complete =====\n");
    return 0;
}
