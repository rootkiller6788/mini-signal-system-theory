/**
 * @file example_mass_spring.c
 * @brief L6: Mass-Spring-Damper System Identification
 *
 * Identifies a second-order mechanical system from PRBS excitation data.
 *
 * Physics: m ẍ + c ẋ + k x = F(t)
 * Transfer function: G(s) = 1 / (m s² + c s + k)
 * Natural frequency: ω_n = √(k/m), damping ratio: ζ = c/(2√(mk))
 *
 * Discrete: ARX(2,2) can exactly represent a 2nd-order system.
 *
 * This example:
 *   1. Generates mass-spring-damper response to PRBS force input
 *   2. Identifies ARX model
 *   3. Recovers modal parameters (ω_n, ζ) from identified discrete poles
 *   4. Validates with frequency response comparison
 *
 * Reference: Ewins (2000) "Modal Testing: Theory, Practice and Application"
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_validation.h"
#include "sysid_signals.h"

int main(void) {
    printf("==========================================\n");
    printf("  Mass-Spring-Damper System ID Example\n");
    printf("==========================================\n\n");

    /* Physical parameters */
    double m = 1.0;       /* kg */
    double c = 0.4;       /* N·s/m (damping) */
    double k = 100.0;     /* N/m (stiffness) */
    double omega_n = sqrt(k / m);  /* 10 rad/s */
    double zeta = c / (2.0 * sqrt(m * k)); /* 0.02 (lightly damped) */
    double f_n = omega_n / (2.0 * M_PI);  /* ~1.59 Hz */

    double Ts = 0.005;    /* 5 ms sampling (200 Hz, well above Nyquist) */
    int N = 1000;

    printf("Physical parameters:\n");
    printf("  m=%.2f kg, c=%.2f Ns/m, k=%.1f N/m\n", m, c, k);
    printf("  ω_n=%.2f rad/s, f_n=%.2f Hz, ζ=%.4f\n\n", omega_n, f_n, zeta);
    printf("Sampling: Ts=%.3f s, N=%d (%.1f s of data)\n\n", Ts, N, N*Ts);

    /* Generate PRBS excitation */
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));
    sysid_prbs_generate_std(10, u, N, 3, 10.0); /* ±10N PRBS */

    /* Simulate continuous system via discrete state-space */
    /* Discretize: [x; v]_{k+1} = Φ [x; v]_k + Γ F_k
     * A_cont = [[0, 1], [-k/m, -c/m]]
     */
    double A_c[4] = {0.0, 1.0, -k/m, -c/m};
    double B_c[2] = {0.0, 1.0/m};
    /* A_c, B_c documented for explicit continuous-time state-space reference */
    (void)A_c; (void)B_c;

    /* Simple Euler integration for small Ts */
    double x = 0.0, v = 0.0;
    for (int k = 0; k < N; k++) {
        /* Measurement: position x with noise */
        y[k] = x + 0.001 * ((double)rand() / RAND_MAX - 0.5); /* 1mm noise */

        /* Euler integration */
        double a = (-k * x - c * v + u[k]) / m;
        x = x + v * Ts;
        v = v + a * Ts;
    }

    printf("Data generated via Euler simulation.\n");

    /* Estimate ARX(2,2) model */
    int na = 2, nb = 2, nk = 1;
    int nparam = na + nb;

    int max_nrows = N - nk - nb + 1;
    double *Phi = (double *)calloc(max_nrows * nparam, sizeof(double));
    double *Y_vec = (double *)calloc(max_nrows, sizeof(double));
    int nrows;
    sysid_build_regression_arx(u, y, N, na, nb, nk, Phi, Y_vec, &nrows);

    double theta[4], sigma2;
    if (sysid_ls_solve(Phi, Y_vec, nrows, nparam, theta, NULL, &sigma2) != 0) {
        printf("LS estimation failed.\n");
        free(u); free(y); free(Phi); free(Y_vec);
        return 1;
    }
    free(Phi); free(Y_vec);

    /* Build ARX model */
    sysid_arx_t model;
    sysid_arx_alloc(&model, na, nb, nk, Ts);
    model.A[1] = theta[0];
    model.A[2] = theta[1];
    model.B[0] = theta[2];
    model.B[1] = theta[3];

    printf("\n--- Identified ARX Model ---\n");
    printf("A(q) = 1 + %.4f q^-1 + %.4f q^-2\n", theta[0], theta[1]);
    printf("B(q) = %.6f q^-1 + %.6f q^-2\n", theta[2], theta[3]);

    /* Compute poles of the identified model */
    double complex poles[2];
    sysid_arx_poles(&model, poles, 2);

    printf("\n--- Discrete Poles ---\n");
    for (int i = 0; i < 2; i++) {
        printf("z%d = %.6f %+.6f j  (|z| = %.6f)\n",
               i+1, creal(poles[i]), cimag(poles[i]), cabs(poles[i]));
    }

    /* Convert to continuous-time poles: s = ln(z)/Ts */
    printf("\n--- Continuous-Time Poles ---\n");
    double complex s_poles[2];
    for (int i = 0; i < 2; i++) {
        s_poles[i] = clog(poles[i]) / Ts;
        printf("s%d = %.3f %+.3f j rad/s\n",
               i+1, creal(s_poles[i]), cimag(s_poles[i]));
    }

    /* Modal parameters from complex conjugate poles:
     * s = -ζ ω_n ± j ω_n √(1-ζ²)
     * ω_n = |s|, ζ = -Re(s) / ω_n
     */
    double omega_n_est = cabs(s_poles[0]);
    double zeta_est = -creal(s_poles[0]) / omega_n_est;

    printf("\n--- Modal Parameter Recovery ---\n");
    printf("Natural frequency ω_n: %.3f rad/s (true: %.3f, error: %.2f%%)\n",
           omega_n_est, omega_n, fabs(omega_n_est - omega_n) / omega_n * 100.0);
    printf("Damping ratio ζ:     %.4f       (true: %.4f, error: %.2f%%)\n",
           zeta_est, zeta, fabs(zeta_est - zeta) / (zeta + 1e-10) * 100.0);
    printf("Stiffness k_est:     %.2f N/m   (true: %.1f)\n",
           omega_n_est * omega_n_est * m, k);
    printf("Damping c_est:       %.4f Ns/m  (true: %.2f)\n",
           2.0 * zeta_est * omega_n_est * m, c);

    /* Validation */
    double *ypred = (double *)calloc(N, sizeof(double));
    sysid_arx_predict(&model, u, y, ypred, N);

    sysid_fit_result_t result;
    result.N_data = N;
    result.nparam = nparam;
    sysid_compute_fit_metrics(y, ypred, N, &result);
    sysid_fill_ics(&result, sigma2);

    printf("\n--- Validation ---\n");
    printf("NRMSE Fit: %.2f %%\n", result.NRMSE_fit);
    printf("R-squared: %.4f\n", result.R_squared);
    printf("AIC:       %.2f\n", result.AIC);
    printf("AICc:      %.2f\n", result.AICc);

    /* Frequency response comparison */
    printf("\n--- Frequency Response (at 5 frequencies) ---\n");
    double omega_test[5] = {1.0, 5.0, omega_n * 0.7, omega_n, omega_n * 1.5};
    double complex G_est[5];
    sysid_arx_freqresp(&model, omega_test, 5, G_est);

    printf("ω [rad/s]   |G_est|     ∠G_est[°]   |G_true|  \n");
    printf("----------  ---------  ----------  ----------\n");
    for (int i = 0; i < 5; i++) {
        double complex s = I * omega_test[i];
        double complex G_true = 1.0 / (m * s * s + c * s + k);
        printf("%9.2f  %9.6f  %9.2f  %9.6f\n",
               omega_test[i], cabs(G_est[i]),
               atan2(cimag(G_est[i]), creal(G_est[i])) * 180.0 / M_PI,
               cabs(G_true));
    }

    /* Cross-validation for robustness check */
    printf("\n--- 5-fold Cross-Validation ---\n");
    double avg_fit, std_fit;
    sysid_cross_validate_arx(u, y, N, 2, 2, 1, 5, &avg_fit, &std_fit);
    printf("Average NRMSE fit across folds: %.2f ± %.2f %%\n", avg_fit, std_fit);

    /* Cleanup */
    free(u); free(y); free(ypred);
    sysid_arx_free(&model);

    printf("\nDone.\n");
    return 0;
}
