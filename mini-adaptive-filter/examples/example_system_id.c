/**
 * @file    example_system_id.c
 * @brief   Example: System Identification using LMS
 *
 * Demonstrates how to use the adaptive filter library to identify
 * an unknown linear system. The unknown plant is a 5-tap FIR filter,
 * and LMS is used to estimate its impulse response from input-output
 * observations.
 *
 * Knowledge: L6 (System Identification)
 * Build: gcc -I../include -o example_system_id example_system_id.c
 *        ../src/adaptive_filter_core.c ../src/adaptive_lms.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adaptive_filter.h"

int main(void) {
    /* Unknown system: h = [0.5, -0.3, 0.8, -0.1, 0.2] */
    size_t order = 5;
    double h_unknown[5] = {0.5, -0.3, 0.8, -0.1, 0.2};

    /* Training parameters */
    size_t n_samples = 5000;
    double mu = 0.005;
    double *x = (double *)malloc(n_samples * sizeof(double));
    double *d = (double *)malloc(n_samples * sizeof(double));

    /* Generate white noise excitation and system response */
    double state[5] = {0.0};
    size_t n, i;
    srand(42);

    printf("System Identification using LMS\n");
    printf("Unknown system: h = [%.1f, %.1f, %.1f, %.1f, %.1f]\n",
           h_unknown[0], h_unknown[1], h_unknown[2], h_unknown[3], h_unknown[4]);
    printf("Training with %zu samples, mu=%.4f\n\n", n_samples, mu);

    for (n = 0; n < n_samples; n++) {
        x[n] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;

        /* Shift state and compute output */
        for (i = order - 1; i > 0; i--) state[i] = state[i-1];
        state[0] = x[n];
        d[n] = 0.0;
        for (i = 0; i < order; i++) d[n] += h_unknown[i] * state[i];
    }

    /* Run system identification */
    double w_est[5];
    double final_mse;
    af_app_system_identify_lms(x, d, n_samples, order, mu, w_est, &final_mse);

    /* Report results */
    printf("Estimated system: w = [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
           w_est[0], w_est[1], w_est[2], w_est[3], w_est[4]);
    printf("\nEstimation errors:\n");
    for (i = 0; i < order; i++) {
        printf("  w[%zu] - h[%zu] = %+.6f\n", i, i, w_est[i] - h_unknown[i]);
    }
    printf("Final MSE: %.6f\n", final_mse);

    free(x); free(d);
    return 0;
}
