/**
 * @file example_dc_motor.c
 * @brief L6/L7: DC Motor System Identification
 *
 * Identifies the transfer function of a DC motor from I/O data.
 *
 * Physics: A DC motor is modeled as a second-order system:
 *   G(s) = K / [s (τ_m s + 1)]
 * where K = speed constant [rad/s/V], τ_m = mechanical time constant [s].
 *
 * Discretized (Euler, Ts=0.01s):
 *   ω(k) = a1 ω(k-1) + a2 ω(k-2) + b1 V(k-1) + b2 V(k-2)
 *
 * This example:
 *   1. Generates synthetic DC motor data (PRBS excitation)
 *   2. Estimates ARX model via Least Squares
 *   3. Validates with NRMSE fit, AICc, residual analysis
 *   4. Computes physical parameters K, τ_m from discrete estimates
 *
 * Reference: Ljung (1999) §1.3; Franklin, Powell, Emami-Naeini (2019)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_validation.h"
#include "sysid_signals.h"

int main(void) {
    printf("==============================\n");
    printf("  DC Motor System ID Example\n");
    printf("==============================\n\n");

    /* Motor parameters (ground truth) */
    double K_true = 30.0;     /* rad/s/V */
    double tau_m_true = 0.05; /* 50 ms mechanical time constant */
    double Ts = 0.001;        /* 1 kHz sampling */
    int N = 2000;

    printf("True motor: K=%.1f rad/s/V, tau=%.3f s\n", K_true, tau_m_true);
    printf("Sampling: Ts=%.4f s, N=%d samples\n\n", Ts, N);

    /* Generate synthetic data */
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    /* PRBS excitation (±5V) */
    sysid_prbs_generate_std(10, u, N, 5, 5.0);

    /* Discretize motor model: ω(k+1) = Φ ω(k) + Γ V(k)
     * Continuous: ẋ = Ax + Bu, x=[ω,α]ᵀ, A=[0 1; 0 -1/τ], B=[0; K/τ]
     * Discrete: Φ = exp(A*Ts), Γ = ∫exp(Aτ)B dτ
     *
     * Simplified: use bilinear approximation for ARX
     * The motor is type-1 (has integrator): A(q) has a root at q=1.
     *   (1 - q^{-1})(1 - a q^{-1}) ω = b₀ q^{-1} (1 + b₁ q^{-1}) V
     */
    double omega = 0.0;      /* current speed */
    double alpha = 0.0;      /* acceleration state */
    double A_disc[4];        /* 2x2 discrete state matrix */

    /* Discretize: A_cont = [[0, 1], [0, -1/tau_m]] */
    double lam = exp(-Ts / tau_m_true);
    A_disc[0] = 1.0;
    A_disc[1] = tau_m_true * (1.0 - lam);  /* ∫₀ᵀˢ e^{Aτ} dτ * B[0] */
    A_disc[2] = 0.0;
    A_disc[3] = lam;
    /* B_disc = [K*Ts - K*tau_m*(1-lam); K*(1-lam)] approximately */
    double b_alpha = K_true * (1.0 - lam);
    double b_omega = K_true * Ts - K_true * tau_m_true * (1.0 - lam);

    for (int k = 0; k < N; k++) {
        /* Output: speed measurement */
        y[k] = omega + 0.05 * ((double)rand() / RAND_MAX - 0.5); /* noise std=0.05 */

        /* State update */
        double omega_new = A_disc[0] * omega + A_disc[1] * alpha + b_omega * u[k];
        double alpha_new = A_disc[2] * omega + A_disc[3] * alpha + b_alpha * u[k];
        omega = omega_new;
        alpha = alpha_new;
    }

    printf("Data generated.\n");

    /* Estimate ARX model of order na=2, nb=2, nk=1 */
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

    /* Build identified ARX model */
    sysid_arx_t model;
    sysid_arx_alloc(&model, na, nb, nk, Ts);
    model.A[1] = theta[0];
    model.A[2] = theta[1];
    model.B[0] = theta[2];
    model.B[1] = theta[3];

    printf("\n--- Identified ARX Model ---\n");
    printf("A(q) = 1 + %.4f q^-1 + %.4f q^-2\n", theta[0], theta[1]);
    printf("B(q) = %.4f q^-1 + %.4f q^-2\n", theta[2], theta[3]);

    /* Validate model */
    double *ypred = (double *)calloc(N, sizeof(double));
    sysid_arx_predict(&model, u, y, ypred, N);

    sysid_fit_result_t result;
    result.N_data = N;
    result.nparam = nparam;
    sysid_compute_fit_metrics(y, ypred, N, &result);
    sysid_fill_ics(&result, sigma2);

    printf("\n--- Validation Results ---\n");
    printf("NRMSE Fit:      %.2f %%\n", result.NRMSE_fit);
    printf("R-squared:      %.4f\n", result.R_squared);
    printf("RMSE:           %.4f\n", result.RMSE);
    printf("Noise variance: %.6f\n", sigma2);
    printf("AIC:            %.2f\n", result.AIC);
    printf("AICc:           %.2f\n", result.AICc);
    printf("BIC:             %.2f\n", result.BIC);

    /* Residual whiteness test */
    /* Compute residuals for whiteness test */
    double *residuals = (double *)calloc(N, sizeof(double));
    sysid_compute_residuals(y, ypred, N, residuals);
    double Q_res = sysid_ljung_box_q(residuals, N, 20);
    printf("Ljung-Box Q(20): %.2f\n", Q_res);
    printf("Residual mean:   %.6f\n", result.residual_mean);

    /* Recover physical parameters */
    /* A(q) = (1 - q^{-1})(1 - a₁ q^{-1}) = 1 - (1+a₁)q^{-1} + a₁ q^{-2}
     * So: a₁ = A[2] = theta[1], and A[1] = -(1+a₁) ≈ theta[0]
     * τ_m ≈ -Ts / ln(a₁), K ≈ (b₀+b₁) / (1-a₁)
     */
    double a1_pole = theta[1];
    double K_est, tau_est;
    if (fabs(a1_pole) < 0.999) {
        tau_est = -Ts / log(fabs(a1_pole));
        K_est = (theta[2] + theta[3]) / (1.0 - a1_pole) / Ts;
        printf("\n--- Recovered Physical Parameters ---\n");
        printf("K (speed constant):  %.2f rad/s/V  (true: %.2f)\n", K_est, K_true);
        printf("tau_m (mech. time):  %.4f s        (true: %.4f)\n", tau_est, tau_m_true);
        printf("K estimate error:    %.1f %%\n", fabs(K_est - K_true) / K_true * 100.0);
        printf("tau_m estimate error: %.1f %%\n", fabs(tau_est - tau_m_true) / tau_m_true * 100.0);
    }

    /* DC gain: should be infinite for type-1 system (integrator)
     * The steady-state speed per volt is K. Check at k=1000 */
    printf("\nSteady-state speed check (at 5V): expected %.1f rad/s\n", K_true * 5.0);

    /* Model order selection via AICc scan */
    int best_na, best_nb;
    double best_aicc;
    sysid_order_selection_aicc(u, y, N, nk, 4, 4, &best_na, &best_nb, &best_aicc);
    printf("\n--- Order Selection (AICc) ---\n");
    printf("Best ARX order: na=%d, nb=%d, AICc=%.2f\n", best_na, best_nb, best_aicc);

    /* Cleanup */
    free(u); free(y); free(ypred); free(residuals);
    sysid_arx_free(&model);

    printf("\nDone.\n");
    return 0;
}
