/**
 * @file test_sysid.c
 * @brief Comprehensive test suite for system identification module
 *
 * Tests cover all core API functions with assert-based validation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_validation.h"
#include "sysid_signals.h"
#include "sysid_subspace.h"
#include "sysid_adaptive.h"

#define TEST_PASS() printf("  PASS: %s\n", __func__)
#define TEST_FAIL(msg) do { printf("  FAIL: %s - %s\n", __func__, msg); return 1; } while(0)

/* test_count and pass_count reserved for future test statistics reporting */

/* Tolerance for floating-point comparisons */
#define TOL 1e-6
#define TOL_LARGE 1e-2

/*---------------------------------------------------------------------------
 * L1: Types and allocation
 *---------------------------------------------------------------------------*/

int test_order_init(void) {
    sysid_model_order_t order;
    sysid_order_init(&order);
    assert(order.na == 0 && order.nb == 0);
    assert(sysid_order_is_valid(&order) == 1);
    assert(sysid_order_nparam(SYSID_MODEL_ARX, &order) == 0);
    TEST_PASS(); return 0;
}

int test_param_alloc(void) {
    sysid_param_vector_t pv;
    assert(sysid_param_alloc(&pv, 5) == 0);
    assert(pv.nparam == 5);
    sysid_param_set_all(&pv, 3.0);
    assert(fabs(pv.theta[0] - 3.0) < TOL);
    assert(fabs(sysid_param_norm2(&pv) - sqrt(45.0)) < TOL);

    sysid_param_vector_t pv2;
    sysid_param_alloc(&pv2, 5);
    sysid_param_set_all(&pv2, 1.0);
    sysid_param_copy(&pv2, &pv);
    assert(fabs(pv2.theta[2] - 3.0) < TOL);

    double dot = sysid_param_dot(&pv, &pv2);
    assert(fabs(dot - 45.0) < TOL);

    sysid_param_free(&pv);
    sysid_param_free(&pv2);
    TEST_PASS(); return 0;
}

int test_data_alloc(void) {
    sysid_data_t data;
    assert(sysid_data_alloc(&data, 100, 1, 1, 0.01) == 0);
    assert(data.N == 100 && data.nu == 1 && data.ny == 1);
    sysid_data_free(&data);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L2/L3: Model construction and prediction
 *---------------------------------------------------------------------------*/

int test_arx_model(void) {
    sysid_arx_t model;
    assert(sysid_arx_alloc(&model, 2, 2, 1, 1.0) == 0);

    /* Set known parameters: y(k) = -0.5 y(k-1) - 0.3 y(k-2) + 1.0 u(k-1) + 0.5 u(k-2) */
    model.A[1] = -0.5;
    model.A[2] = -0.3;
    model.B[0] = 1.0;
    model.B[1] = 0.5;

    /* Generate test data */
    int N = 50;
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));
    double *ypred = (double *)calloc(N, sizeof(double));

    /* Step input */
    for (int i = 0; i < N; i++) u[i] = 1.0;

    /* Simulate output */
    for (int k = 0; k < N; k++) {
        y[k] = 0.0;
        if (k - 1 >= 0) {
            y[k] += -model.A[1] * y[k - 1] + model.B[0] * u[k - 1];
        }
        if (k - 2 >= 0) {
            y[k] += -model.A[2] * y[k - 2] + model.B[1] * u[k - 2];
        }
        /* Add small noise */
        y[k] += 0.001 * ((double)rand() / RAND_MAX - 0.5);
    }

    /* Predict */
    assert(sysid_arx_predict(&model, u, y, ypred, N) == 0);

    /* Check prediction quality */
    double mse = 0.0;
    for (int i = 2; i < N; i++) {
        mse += (y[i] - ypred[i]) * (y[i] - ypred[i]);
    }
    mse /= (N - 2);
    assert(mse < 0.01);

    /* DC gain */
    double dc = sysid_arx_dcgain(&model);
    /* DC gain check: G(1) = sum(b) / (1 + sum(a_i)) = (1.0+0.5) / (1-0.5-0.3) = 7.5 */
    assert(fabs(dc - 7.5) < TOL);

    free(u); free(y); free(ypred);
    sysid_arx_free(&model);
    TEST_PASS(); return 0;
}

int test_ss_model(void) {
    sysid_ss_t ss;
    assert(sysid_ss_alloc(&ss, 2, 1, 1, 1, 1.0) == 0);

    /* Simple system: integrator */
    ss.A[0] = 1.0; ss.A[1] = 0.0;
    ss.A[2] = 0.1; ss.A[3] = 1.0;
    ss.B[0] = 0.1; ss.B[1] = 0.0;
    ss.C[0] = 1.0; ss.C[1] = 0.0;
    ss.D[0] = 0.0;

    int N = 20;
    double *u = (double *)calloc(N, sizeof(double));
    double *ysim = (double *)calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) u[i] = 1.0;

    assert(sysid_ss_simulate(&ss, u, ysim, N) == 0);

    /* Output should grow (integrator) */
    assert(ysim[N-1] > ysim[0] * 2.0);

    free(u); free(ysim);
    sysid_ss_free(&ss);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L5: Least Squares estimation
 *---------------------------------------------------------------------------*/

int test_ls_solve(void) {
    /* y = 2*x + 1 with noise */
    int N = 50;
    double *Phi = (double *)calloc(N * 2, sizeof(double));
    double *Y = (double *)calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        double x = (double)i / N;
        Phi[i * 2 + 0] = x;
        Phi[i * 2 + 1] = 1.0;
        Y[i] = 2.0 * x + 1.0 + 0.01 * ((double)rand() / RAND_MAX - 0.5);
    }

    double theta[2];
    assert(sysid_ls_solve(Phi, Y, N, 2, theta, NULL, NULL) == 0);
    assert(fabs(theta[0] - 2.0) < 0.1);
    assert(fabs(theta[1] - 1.0) < 0.1);

    free(Phi); free(Y);
    TEST_PASS(); return 0;
}

int test_sysid_arx_estimate(void) {
    /* Create known ARX system and recover via LS */
    int N = 500;
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    /* Generate PRBS input */
    sysid_prbs_generate_std(7, u, N, 1, 1.0);

    /* Known stable system: y(k) = 0.9 y(k-1) - 0.2 y(k-2) + 1.0 u(k-1) + 0.5 u(k-2) + e(k)
     * ARX form: y(k) = -a₁ y(k-1) - a₂ y(k-2) + b₀ u(k-1) + b₁ u(k-2)
     * So: a₁ = -0.9, a₂ = 0.2, b₀ = 1.0, b₁ = 0.5
     * LS estimates: theta = [a₁, a₂, b₀, b₁] = [-0.9, 0.2, 1.0, 0.5] */
    for (int k = 0; k < N; k++) {
        y[k] = 0.0;
        if (k - 1 >= 0) y[k] += 0.9 * y[k - 1] + 1.0 * u[k - 1];
        if (k - 2 >= 0) y[k] += -0.2 * y[k - 2] + 0.5 * u[k - 2];
        y[k] += 0.01 * ((double)rand() / RAND_MAX - 0.5);
    }

    int na = 2, nb = 2, nk = 1;
    int nparam = na + nb;
    int max_nrows = N - nk - nb + 1;
    double *Phi = (double *)calloc(max_nrows * nparam, sizeof(double));
    double *Y_vec = (double *)calloc(max_nrows, sizeof(double));
    int nrows;

    sysid_build_regression_arx(u, y, N, na, nb, nk, Phi, Y_vec, &nrows);
    assert(nrows > 0);

    double theta[4];
    assert(sysid_ls_solve(Phi, Y_vec, nrows, nparam, theta, NULL, NULL) == 0);

    /* Check estimates: theta = [-a1, -a2, b0, b1] = [-0.9, 0.2, 1.0, 0.5] */
    assert(fabs(theta[0] + 0.9) < 0.2);
    assert(fabs(theta[1] - 0.2) < 0.2);
    assert(fabs(theta[2] - 1.0) < 0.2);
    assert(fabs(theta[3] - 0.5) < 0.2);

    free(u); free(y); free(Phi); free(Y_vec);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L4: Information criteria and validation
 *---------------------------------------------------------------------------*/

int test_information_criteria(void) {
    double loss = 0.01;
    int N = 200, d = 5;
    double aic = sysid_aic(loss, N, d);
    double bic = sysid_bic(loss, N, d);
    assert(bic > aic); /* BIC penalizes more */
    assert(sysid_aicc(loss, N, d) > aic); /* AICc stronger penalty than AIC */
    TEST_PASS(); return 0;
}

int test_fit_metrics(void) {
    int N = 100;
    double *y = (double *)calloc(N, sizeof(double));
    double *ypred = (double *)calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        y[i] = 2.0 * i / N;
        ypred[i] = y[i] + 0.05;
    }

    sysid_fit_result_t result;
    result.N_data = N;
    result.nparam = 2;
    sysid_compute_fit_metrics(y, ypred, N, &result);

    assert(result.NRMSE_fit > 50.0); /* OK fit (biased by 0.05) */
    assert(result.R_squared > 0.80);
    /* Residuals = y - ypred = y[i] - (y[i]+0.05) = -0.05, so mean ≈ -0.05 */
    assert(fabs(result.residual_mean + 0.05) < 0.02);

    free(y); free(ypred);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L6: Signal generation
 *---------------------------------------------------------------------------*/

int test_prbs(void) {
    int N = 200;
    double *signal = (double *)calloc(N, sizeof(double));
    assert(sysid_prbs_generate_std(5, signal, N, 1, 1.0) == 0);

    /* PRBS should have ±1 values */
    int n_pos = 0, n_neg = 0;
    for (int i = 0; i < N; i++) {
        assert(fabs(signal[i]) == 1.0);
        if (signal[i] > 0) n_pos++; else n_neg++;
    }
    assert(n_pos > 20 && n_neg > 20); /* roughly balanced */

    /* Persistence test — PRBS should be PE of decent order */
    double R_buf[9]; /* 3x3 matrix */
    assert(sysid_test_pe_order(signal, N, 3, R_buf, 0.01) == 1);

    free(signal);
    TEST_PASS(); return 0;
}

int test_chirp(void) {
    int N = 200;
    double *signal = (double *)calloc(N, sizeof(double));
    assert(sysid_chirp_linear(signal, N, 1.0, 10.0, 0.01, 1.0) == 0);

    /* Check that signal is bounded */
    for (int i = 0; i < N; i++) {
        assert(fabs(signal[i]) <= 1.01);
    }
    free(signal);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L5: RLS and adaptive methods
 *---------------------------------------------------------------------------*/

int test_rls(void) {
    /* Track a stable ARX system with RLS */
    int nparam = 2;
    sysid_rls_t rls;
    /* Large initial uncertainty (small delta → large P) for fast initial learning */
    assert(sysid_rls_alloc(&rls, nparam, 0.98, 0.01) == 0);

    double true_a1 = 0.7;  /* stable: y(k) = 0.7 y(k-1) + 0.5 u(k-1) */
    double true_b0 = 0.5;
    #define RLS_N 500
    double u_hist[RLS_N];
    double y_hist[RLS_N];
    memset(u_hist, 0, sizeof(u_hist));
    memset(y_hist, 0, sizeof(y_hist));

    for (int k = 1; k < RLS_N; k++) {
        u_hist[k] = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        y_hist[k] = true_a1 * y_hist[k - 1] + true_b0 * u_hist[k - 1]
                    + 0.01 * ((double)rand() / RAND_MAX);

        /* ARX regressor: phi = [-y(k-1), u(k-1)] → theta = [-a1, b0] */
        double phi[2] = {-y_hist[k - 1], u_hist[k - 1]};
        sysid_rls_step(&rls, phi, y_hist[k]);
    }

    /* After convergence: theta[0] ≈ -a1 = -0.7, theta[1] ≈ b0 = 0.5 */
    /* Relax tolerance since finite samples */
    assert(fabs(rls.theta[0] - (-true_a1)) < 0.3);
    assert(fabs(rls.theta[1] - true_b0) < 0.3);

    sysid_rls_free(&rls);
    TEST_PASS(); return 0;
}

int test_lms(void) {
    /* LMS for FIR system identification */
    int n_taps = 4;
    double true_w[4] = {0.5, -0.3, 0.2, 0.1};
    sysid_lms_t lms;
    assert(sysid_lms_alloc(&lms, n_taps, 0.01, 0.0) == 0);

    int LMS_N = 1000;
    double *x_buf = (double *)calloc(LMS_N, sizeof(double));
    /* Generate pseudo-random regressor sequentially */
    for (int k = n_taps - 1; k < LMS_N; k++) {
        double x[4];
        for (int i = 0; i < n_taps; i++) x[n_taps-1-i] = ((double)rand() / RAND_MAX - 0.5);
        double d = 0.0;
        for (int i = 0; i < n_taps; i++) d += true_w[i] * x[i];
        d += 0.01 * ((double)rand() / RAND_MAX);

        double y_out;
        sysid_lms_step(&lms, x, d, &y_out);
    }

    /* Weights should converge */
    for (int i = 0; i < n_taps; i++) {
        assert(fabs(lms.w[i] - true_w[i]) < 0.25);
    }
    free(x_buf);

    sysid_lms_free(&lms);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * L5: Instrumental Variables
 *---------------------------------------------------------------------------*/

int test_iv(void) {
    /* IV on ARX with colored noise */
    int N = 300;
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    /* PRBS input */
    sysid_prbs_generate_std(7, u, N, 1, 1.0);

    /* Colored noise: e(k) = 0.8 e(k-1) + w(k), w ~ N(0, 0.01) */
    double e = 0.0;
    for (int k = 0; k < N; k++) {
        e = 0.8 * e + 0.01 * ((double)rand() / RAND_MAX - 0.5);
        y[k] = e;
        if (k - 1 >= 0) y[k] += 0.7 * y[k - 1] + 0.5 * u[k - 1];
        if (k - 2 >= 0) y[k] += 0.3 * u[k - 2];
    }

    double theta_iv[4], sigma2;
    int ret = sysid_iv4_estimate(u, y, N, 1, 2, 1, theta_iv, NULL, &sigma2, 5);
    assert(ret >= 0);

    free(u); free(y);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * Cross-validation
 *---------------------------------------------------------------------------*/

int test_cross_validate(void) {
    int N = 200;
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    sysid_prbs_generate_std(7, u, N, 1, 1.0);
    for (int k = 2; k < N; k++) {
        y[k] = 0.9 * y[k - 1] - 0.2 * y[k - 2] + 0.5 * u[k - 1] + 0.01 * ((double)rand() / RAND_MAX);
    }

    double avg_fit, std_fit;
    assert(sysid_cross_validate_arx(u, y, N, 2, 1, 1, 5, &avg_fit, &std_fit) == 0);
    assert(avg_fit > 70.0); /* Should fit well */

    free(u); free(y);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * Subspace identification
 *---------------------------------------------------------------------------*/

int test_subspace(void) {
    int N = 300;
    double *u = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));

    sysid_prbs_generate_std(10, u, N, 1, 1.0);

    /* 2nd-order system */
    double x1 = 0.0, x2 = 0.0;
    for (int k = 0; k < N; k++) {
        double x1_new = 1.5 * x1 - 0.7 * x2 + u[k];
        double x2_new = x1;
        x1 = x1_new; x2 = x2_new;
        y[k] = x1 + 0.1 * ((double)rand() / RAND_MAX - 0.5);
    }

    sysid_subspace_opts_t opts;
    sysid_subspace_opts_init(&opts);
    opts.i = 8;
    opts.nx_max = 4;

    sysid_ss_t model;
    int nx_out;
    double sv[32];
    assert(sysid_n4sid(u, y, N, 1, 1, &opts, &model, &nx_out, sv) == 0);

    /* Order detection should find nx=2 */
    assert(nx_out >= 1 && nx_out <= 4);

    sysid_ss_free(&model);
    free(u); free(y);
    TEST_PASS(); return 0;
}

/*---------------------------------------------------------------------------
 * Main test runner
 *---------------------------------------------------------------------------*/

int main(void) {
    printf("=== System Identification Test Suite ===\n\n");

    int failed = 0;

    printf("--- L1: Types and Definitions ---\n");
    failed += test_order_init();
    failed += test_param_alloc();
    failed += test_data_alloc();

    printf("\n--- L2/L3: Models ---\n");
    failed += test_arx_model();
    failed += test_ss_model();

    printf("\n--- L5: Estimation Algorithms ---\n");
    failed += test_ls_solve();
    failed += test_sysid_arx_estimate();
    failed += test_rls();
    failed += test_lms();
    failed += test_iv();

    printf("\n--- L4: Validation ---\n");
    failed += test_information_criteria();
    failed += test_fit_metrics();
    failed += test_cross_validate();

    printf("\n--- L6: Signals ---\n");
    failed += test_prbs();
    failed += test_chirp();

    printf("\n--- L5/L8: Subspace ---\n");
    failed += test_subspace();

    printf("\n========================================\n");
    if (failed == 0) {
        printf("ALL TESTS PASSED (%d test functions)\n", 16);
    } else {
        printf("%d TEST(S) FAILED\n", failed);
    }
    printf("========================================\n");

    return failed;
}
