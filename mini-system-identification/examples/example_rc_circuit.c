/**
 * @file example_rc_circuit.c
 * @brief L6: RC Circuit System Identification
 *
 * Identifies an RC low-pass filter from step response data.
 *
 * Physics: V_out(t) = V_in * (1 - e^{-t/τ}), τ = R·C
 * Discrete: y(k) = a y(k-1) + b u(k-1), where a=e^{-Ts/τ}, b=1-a
 *
 * This example:
 *   1. Generates RC circuit step response with measurement noise
 *   2. Estimates first-order ARX model via LS
 *   3. Recovers R·C time constant
 *   4. Compares with nonlinear least squares (PEM for OE model)
 *
 * Reference: Ljung (1999) Example 4.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_validation.h"

int main(void) {
    printf("================================\n");
    printf("  RC Circuit System ID Example\n");
    printf("================================\n\n");

    /* Circuit parameters */
    double R = 10000.0;    /* 10 kΩ */
    double C = 1e-6;       /* 1 µF */
    double tau_true = R * C; /* 10 ms */
    double Ts = 0.001;     /* 1 kHz sampling */
    int N = 500;

    printf("True RC circuit: R=%.0f ohm, C=%.1e F, tau=%.4f s\n", R, C, tau_true);
    printf("Sampling: Ts=%.4f s, N=%d\n\n", Ts, N);

    /* Generate step response */
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    /* Step input at k=20 */
    for (int k = 20; k < N; k++) u[k] = 5.0; /* 5V step */

    /* Generate true response with noise */
    double y_true = 0.0;
    double a_true = exp(-Ts / tau_true);
    double b_true = 1.0 - a_true;

    for (int k = 0; k < N; k++) {
        y_true = a_true * y_true + b_true * u[k];
        y[k] = y_true + 0.005 * ((double)rand() / RAND_MAX - 0.5); /* 5mV noise */
    }

    printf("Step response generated.\n\n");

    /* Method 1: ARX LS estimation */
    int na = 1, nb = 1, nk = 1;
    int nparam = na + nb;

    int max_nrows = N - nk - nb + 1;
    double *Phi = (double *)calloc(max_nrows * nparam, sizeof(double));
    double *Y_vec = (double *)calloc(max_nrows, sizeof(double));
    int nrows;
    sysid_build_regression_arx(u, y, N, na, nb, nk, Phi, Y_vec, &nrows);

    double theta[2], sigma2_ls;
    sysid_ls_solve(Phi, Y_vec, nrows, nparam, theta, NULL, &sigma2_ls);

    double a_est = -theta[0]; /* ARX model: theta[0] = -a1 */
    double b_est = theta[1];
    double tau_est_ls = -Ts / log(a_est);

    printf("--- Method 1: Least Squares (ARX) ---\n");
    printf("Estimated ARX: y(k) = %.4f y(k-1) + %.6f u(k-1)\n", a_est, b_est);
    printf("Estimated tau: %.6f s (true: %.6f s, error: %.2f%%)\n",
           tau_est_ls, tau_true, fabs(tau_est_ls - tau_true) / tau_true * 100.0);

    /* Validate LS model */
    sysid_arx_t model_ls;
    sysid_arx_alloc(&model_ls, na, nb, nk, Ts);
    model_ls.A[1] = theta[0]; /* -a1 */
    model_ls.B[0] = theta[1];

    double *ypred = (double *)calloc(N, sizeof(double));
    sysid_arx_predict(&model_ls, u, y, ypred, N);

    sysid_fit_result_t result;
    result.N_data = N;
    result.nparam = nparam;
    sysid_compute_fit_metrics(y, ypred, N, &result);
    sysid_fill_ics(&result, sigma2_ls);

    printf("NRMSE Fit: %.2f%%, R²: %.4f, AICc: %.2f\n\n",
           result.NRMSE_fit, result.R_squared, result.AICc);

    /* Method 2: Output Error PEM (Gauss-Newton) */
    printf("--- Method 2: PEM (Output Error) ---\n");

    sysid_oe_t model_oe;
    sysid_oe_alloc(&model_oe, nb, 1, nk, Ts);
    model_oe.B[0] = b_est;
    model_oe.F[1] = -a_est; /* F[1] = -a1 (F polynomial: 1 + f1 q^{-1}) */

    sysid_data_t data;
    sysid_data_alloc(&data, N, 1, 1, Ts);
    memcpy(data.u, u, N * sizeof(double));
    memcpy(data.y, y, N * sizeof(double));

    sysid_param_vector_t pv;
    sysid_param_alloc(&pv, nparam);
    pv.theta[0] = b_est; /* B coefficient */
    pv.theta[1] = -a_est; /* F[1] coefficient */

    double V_final;
    int pem_iters = sysid_pem_oe(&model_oe, &data, &pv, NULL, &V_final, 20, 1e-6);

    double b_pem = model_oe.B[0];
    double a_pem = -model_oe.F[1];
    double tau_pem = -Ts / log(a_pem);

    printf("PEM iterations: %d\n", pem_iters);
    printf("Refined OE: y(k)*F(q) = B(q) u(k-nk)\n");
    printf("B = %.6f, F = 1 + %.4f q^-1\n", b_pem, model_oe.F[1]);
    printf("Refined tau: %.6f s (error: %.2f%%)\n",
           tau_pem, fabs(tau_pem - tau_true) / tau_true * 100.0);
    printf("Loss V_N: %.2e\n\n", V_final);

    /* Compare methods using F-test */
    printf("--- Model Comparison ---\n");
    double F_stat = sysid_ftest(sigma2_ls, V_final, nparam, nparam, nrows);
    double p_value = sysid_ftest_pvalue(sigma2_ls, V_final, nparam, nparam, nrows);
    printf("F-statistic: %.4f, p-value: %.4f\n", F_stat, p_value);
    if (p_value < 0.05) {
        printf("Conclusion: PEM significantly better than LS (p < 0.05)\n");
    } else {
        printf("Conclusion: No significant difference (p >= 0.05)\n");
    }

    /* Model order selection: scan na=1..3, nb=1..3 */
    printf("\n--- AICc Order Selection ---\n");
    int best_na, best_nb;
    double best_aicc;
    sysid_order_selection_aicc(u, y, N, nk, 3, 3, &best_na, &best_nb, &best_aicc);
    printf("Best order: na=%d, nb=%d, AICc=%.2f (true: na=%d, nb=%d)\n",
           best_na, best_nb, best_aicc, na, nb);

    /* Cleanup */
    free(u); free(y); free(Phi); free(Y_vec); free(ypred);
    sysid_arx_free(&model_ls);
    sysid_oe_free(&model_oe);
    sysid_data_free(&data);
    sysid_param_free(&pv);

    printf("\nDone.\n");
    return 0;
}
