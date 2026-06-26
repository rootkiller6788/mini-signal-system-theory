/**
 * @file    example_noise_cancel.c
 * @brief   Example: Adaptive Noise Cancellation
 *
 * Demonstrates adaptive noise cancellation: a clean sine wave is
 * corrupted by additive noise. A reference noise signal (correlated
 * with the corrupting noise) is used to cancel the interference.
 *
 * This implements the classic Widrow noise-cancelling architecture
 * (Widrow et al., Proc. IEEE, 1975).
 *
 * Knowledge: L6 (Noise Cancellation), L7 (ECG/Biomedical Application)
 * Build: gcc -I../include -o example_noise_cancel example_noise_cancel.c
 *        ../src/adaptive_filter_core.c ../src/adaptive_lms.c
 *        ../src/adaptive_applications.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adaptive_filter.h"

int main(void) {
    size_t n_samples = 4000;
    size_t order = 32;
    double mu = 0.01;

    /* Generate clean signal: 5 Hz sine wave */
    double *clean = (double *)malloc(n_samples * sizeof(double));
    double *noise = (double *)malloc(n_samples * sizeof(double));
    double *noise_ref = (double *)malloc(n_samples * sizeof(double));
    double *corrupted = (double *)malloc(n_samples * sizeof(double));
    double *cleaned = (double *)malloc(n_samples * sizeof(double));

    size_t n;
    double fs = 1000.0;  /* 1 kHz sampling */

    printf("Adaptive Noise Cancellation Demo\n");
    printf("Clean signal: 5 Hz sine wave\n");
    printf("Noise: filtered white noise (correlated reference available)\n");
    printf("Filter order: %zu, mu=%.3f\n\n", order, mu);

    /* Generate signals */
    for (n = 0; n < n_samples; n++) {
        double t = (double)n / fs;
        clean[n] = sin(2.0 * M_PI * 5.0 * t);
        noise[n] = 0.5 * ((double)rand() / RAND_MAX - 0.5);
        /* Reference: filtered version of noise (simulates sensor pickup) */
        noise_ref[n] = 0.8 * noise[n] +
                       0.2 * (n > 0 ? noise[n-1] : 0.0);
        corrupted[n] = clean[n] + noise_ref[n];
    }

    /* Apply noise cancellation */
    af_app_noise_cancel_lms(corrupted, noise, n_samples, order,
                             mu, cleaned);

    /* Compute SNR improvement */
    double power_clean = 0.0, power_residual = 0.0;
    size_t start = n_samples / 2; /* Steady-state region */
    for (n = start; n < n_samples; n++) {
        power_clean += clean[n] * clean[n];
        double residual = cleaned[n] - clean[n];
        power_residual += residual * residual;
    }
    size_t ss_len = n_samples - start;
    power_clean /= (double)ss_len;
    power_residual /= (double)ss_len;

    double snr_improve = 10.0 * log10(power_clean / (power_residual + 1e-15));
    printf("SNR improvement: %.2f dB\n", snr_improve);
    printf("(Higher is better - noise has been cancelled)\n");

    free(clean); free(noise); free(noise_ref); free(corrupted); free(cleaned);
    return 0;
}
