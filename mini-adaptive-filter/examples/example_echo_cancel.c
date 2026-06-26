/**
 * @file    example_echo_cancel.c
 * @brief   Example: Acoustic Echo Cancellation using NLMS
 *
 * Simulates an acoustic echo path (room impulse response) and uses
 * NLMS to cancel the echo from the microphone signal.
 *
 * The echo path is modeled as an exponentially decaying impulse
 * response, simulating a typical room with reverberation.
 *
 * Knowledge: L6 (Echo Cancellation), L7 (Hands-free telephony)
 * Build: gcc -I../include -o example_echo_cancel example_echo_cancel.c
 *        ../src/adaptive_filter_core.c ../src/adaptive_lms.c
 *        ../src/adaptive_applications.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adaptive_filter.h"

int main(void) {
    size_t n_samples = 8000;
    size_t echo_order = 128;  /* Echo path length */
    double mu = 0.5;           /* NLMS step-size (0 < mu < 2) */
    double fs = 8000.0;        /* 8 kHz (telephony) */

    /* Generate echo path: exponentially decaying */
    double *echo_path = (double *)malloc(echo_order * sizeof(double));
    size_t i;
    for (i = 0; i < echo_order; i++) {
        echo_path[i] = 0.8 * exp(-(double)i / 20.0) *
                       ((i % 10 == 0) ? 1.0 : 0.0); /* Sparse reflections */
    }

    /* Generate far-end speech (simulated as filtered noise) */
    double *far_end = (double *)malloc(n_samples * sizeof(double));
    double *echo = (double *)malloc(n_samples * sizeof(double));
    double *near_end = (double *)malloc(n_samples * sizeof(double));
    double *mic_signal = (double *)malloc(n_samples * sizeof(double));
    double *residual = (double *)malloc(n_samples * sizeof(double));

    size_t n;
    double far_state[128] = {0.0};
    srand(12345);

    printf("Acoustic Echo Cancellation using NLMS\n");
    printf("Echo path: %zu taps (sparse, exponential decay)\n", echo_order);
    printf("NLMS step-size: %.2f\n\n", mu);

    for (n = 0; n < n_samples; n++) {
        /* Far-end signal: bandlimited noise (simulating speech) */
        far_end[n] = 0.7 * ((double)rand() / RAND_MAX - 0.5);
        if (n > 0) far_end[n] = 0.9 * far_end[n] + 0.1 * far_end[n-1];

        /* Generate echo through echo path */
        for (i = echo_order - 1; i > 0; i--) far_state[i] = far_state[i-1];
        far_state[0] = far_end[n];
        echo[n] = 0.0;
        for (i = 0; i < echo_order; i++) echo[n] += echo_path[i] * far_state[i];

        /* Near-end speaker: silence (single-talk scenario) */
        near_end[n] = 0.0;

        /* Microphone picks up echo + near-end */
        mic_signal[n] = echo[n] + near_end[n];
    }

    /* Apply echo cancellation */
    af_app_echo_cancel_nlms(far_end, mic_signal, n_samples,
                             echo_order, mu, residual);

    /* Compute echo return loss enhancement (ERLE) */
    double power_echo = 0.0, power_residual = 0.0;
    size_t start = n_samples / 2; /* Steady-state */
    for (n = start; n < n_samples; n++) {
        power_echo += echo[n] * echo[n];
        power_residual += residual[n] * residual[n];
    }
    size_t ss_len = n_samples - start;
    power_echo /= (double)ss_len;
    power_residual /= (double)ss_len;

    double erle = 10.0 * log10(power_echo / (power_residual + 1e-15));
    printf("ERLE (Echo Return Loss Enhancement): %.2f dB\n", erle);
    printf("(Typical commercial systems achieve 20-30 dB ERLE)\n");

    /* Report first 5 echo path estimate vs true */
    printf("\nEcho path estimates (first 5 taps out of %zu):\n", echo_order);

    free(echo_path); free(far_end); free(echo);
    free(near_end); free(mic_signal); free(residual);
    return 0;
}
