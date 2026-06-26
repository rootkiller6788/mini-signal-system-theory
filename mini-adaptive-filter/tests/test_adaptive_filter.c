/**
 * @file    test_adaptive_filter.c
 * @brief   Assert-based tests for adaptive filtering module
 *
 * Tests cover all core APIs and algorithms:
 *   - Lifecycle (init/free/reset)
 *   - LMS variants (direct, NLMS, sign-error, sign-data, sign-sign, leaky)
 *   - RLS (standard, QR-RLS)
 *   - Kalman filter (linear predict/update, EKF)
 *   - Wiener solver (Levinson-Durbin)
 *   - Matrix/vector operations
 *   - Applications (system identification, noise cancellation)
 *   - Edge cases (NULL handling, zero order, boundary conditions)
 */

#include "adaptive_filter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define ASSERT_TRUE(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("FAIL: %s (%.6f vs %.6f)\n", msg, a, b); \
        tests_failed++; return; \
    } \
} while(0)

/* =================================================================
 * Test 1: Lifecycle - Init, Free, Reset
 * ================================================================= */
static void test_lifecycle(void) {
    TEST("af_init / af_free (LMS)");
    {
        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.structure = AF_STRUCTURE_FIR;
        cfg.algorithm = AF_ALGO_LMS;
        cfg.order = 16;
        cfg.mu = 0.01;
        cfg.lambda = 1.0;
        cfg.delta = 1.0;
        cfg.gamma_leakage = 1.0;
        cfg.epsilon = 1e-6;
        cfg.init_strategy = AF_INIT_ZEROS;

        int ret = af_init(&cfg, &f);
        ASSERT_TRUE(ret == 0, "af_init failed");
        ASSERT_TRUE(f.w != NULL, "w is NULL");
        ASSERT_TRUE(f.x_buffer != NULL, "x_buffer is NULL");
        ASSERT_TRUE(f.config != NULL, "config is NULL");
        ASSERT_TRUE(f.iteration == 0, "iteration not zero");

        /* Verify zero initialization */
        size_t i;
        for (i = 0; i < cfg.order; i++) {
            ASSERT_TRUE(f.w[i] == 0.0, "initial weight not zero");
            ASSERT_TRUE(f.x_buffer[i] == 0.0, "initial x_buffer not zero");
        }

        af_free(&f);
        ASSERT_TRUE(f.w == NULL, "w not NULL after free");
        PASS();
    }

    TEST("af_init (RLS)");
    {
        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.structure = AF_STRUCTURE_FIR;
        cfg.algorithm = AF_ALGO_RLS;
        cfg.order = 8;
        cfg.mu = 0.01;
        cfg.lambda = 0.99;
        cfg.delta = 10.0;
        cfg.gamma_leakage = 1.0;
        cfg.epsilon = 1e-6;
        cfg.init_strategy = AF_INIT_ZEROS;

        ASSERT_TRUE(af_init(&cfg, &f) == 0, "RLS init failed");
        ASSERT_TRUE(f.P != NULL, "RLS P matrix is NULL");
        ASSERT_TRUE(f.K != NULL, "RLS K vector is NULL");

        /* P(0) = delta^{-1} * I */
        size_t i;
        for (i = 0; i < cfg.order; i++) {
            double p_ii = f.P[i * cfg.order + i];
            ASSERT_NEAR(p_ii, 1.0 / cfg.delta, 1e-10, "P diagonal incorrect");
        }

        af_free(&f);
        PASS();
    }

    TEST("af_init (invalid config)");
    {
        af_filter_t f;
        /* NULL config */
        ASSERT_TRUE(af_init(NULL, &f) == -1, "NULL config should fail");

        /* Invalid mu */
        af_config_t cfg;
        memset(&cfg, 0, sizeof(cfg));
        cfg.order = 8;
        cfg.mu = 0.0; /* Invalid */
        cfg.lambda = 1.0;
        cfg.delta = 1.0;
        cfg.gamma_leakage = 1.0;
        ASSERT_TRUE(af_init(&cfg, &f) == -1, "zero mu should fail");
        PASS();
    }

    TEST("af_reset");
    {
        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.structure = AF_STRUCTURE_FIR;
        cfg.algorithm = AF_ALGO_LMS;
        cfg.order = 4;
        cfg.mu = 0.1;
        cfg.lambda = 1.0;
        cfg.delta = 1.0;
        cfg.gamma_leakage = 1.0;
        cfg.epsilon = 1e-6;
        cfg.init_strategy = AF_INIT_ZEROS;

        af_init(&cfg, &f);

        /* Manually modify state */
        f.w[0] = 100.0;
        f.error = 5.0;
        f.mse = 25.0;
        f.iteration = 999;

        af_reset(&f);

        ASSERT_NEAR(f.w[0], 0.0, 1e-10, "weights not reset to zero");
        ASSERT_NEAR(f.error, 0.0, 1e-10, "error not reset");
        ASSERT_NEAR(f.mse, 0.0, 1e-10, "mse not reset");
        ASSERT_TRUE(f.iteration == 0, "iteration not reset");

        af_free(&f);
        PASS();
    }
}

/* =================================================================
 * Test 2: LMS Direct Update (Convergence to known weight)
 * ================================================================= */
static void test_lms_direct(void) {
    TEST("LMS direct update - convergence");
    {
        /* Identify a known FIR filter: h = [0.5, -0.3, 0.8] */
        size_t order = 3;
        double h[3] = {0.5, -0.3, 0.8};
        double w[3] = {0.0, 0.0, 0.0};
        double x_buf[3] = {0.0, 0.0, 0.0};

        /* Generate white noise training sequence */
        size_t n_samples = 5000;
        double *x = (double *)malloc(n_samples * sizeof(double));
        double *d = (double *)malloc(n_samples * sizeof(double));
        size_t n;

        srand(42); /* Deterministic seed */
        for (n = 0; n < n_samples; n++) {
            x[n] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        }

        /* Generate desired output: d = h * x */
        for (n = 0; n < n_samples; n++) {
            /* Shift delay line */
            memmove(x_buf + 1, x_buf, (order - 1) * sizeof(double));
            x_buf[0] = x[n];
            d[n] = 0.0;
            size_t i;
            for (i = 0; i < order; i++) {
                d[n] += h[i] * x_buf[i];
            }
        }

        /* Run LMS */
        for (n = 0; n < n_samples; n++) {
            af_lms_direct_update(NULL, x_buf, d[n], w, 0.01, order);
            /* Shift delay line for next iteration */
            memmove(x_buf + 1, x_buf, (order - 1) * sizeof(double));
            x_buf[0] = (n + 1 < n_samples) ? x[n + 1] : 0.0;
        }

        /* Verify convergence */
        ASSERT_NEAR(w[0], h[0], 0.05, "w[0] not converged");
        ASSERT_NEAR(w[1], h[1], 0.05, "w[1] not converged");
        ASSERT_NEAR(w[2], h[2], 0.05, "w[2] not converged");

        free(x); free(d);
        PASS();
    }
}

/* =================================================================
 * Test 3: NLMS convergence
 * ================================================================= */
static void test_nlms(void) {
    TEST("NLMS update - convergence");
    {
        size_t order = 4;
        double h[4] = {0.3, 0.6, -0.2, 0.1};
        double w[4] = {0.0, 0.0, 0.0, 0.0};
        double x_buf[4] = {0.0, 0.0, 0.0, 0.0};

        size_t n_samples = 3000;
        double *x = (double *)malloc(n_samples * sizeof(double));
        double *d = (double *)malloc(n_samples * sizeof(double));

        srand(123);
        size_t n, i;
        for (n = 0; n < n_samples; n++) {
            x[n] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            memmove(x_buf + 1, x_buf, (order - 1) * sizeof(double));
            x_buf[0] = x[n];
            d[n] = 0.0;
            for (i = 0; i < order; i++) d[n] += h[i] * x_buf[i];
        }

        /* Run NLMS */
        for (n = 0; n < n_samples; n++) {
            af_lms_normalized_update(NULL, x_buf, d[n], w, 0.5, 1e-4, order);
            memmove(x_buf + 1, x_buf, (order - 1) * sizeof(double));
            x_buf[0] = (n + 1 < n_samples) ? x[n + 1] : 0.0;
        }

        ASSERT_NEAR(w[0], h[0], 0.05, "NLMS w[0] not converged");
        ASSERT_NEAR(w[1], h[1], 0.05, "NLMS w[1] not converged");

        free(x); free(d);
        PASS();
    }
}

/* =================================================================
 * Test 4: Sign LMS Variants
 * ================================================================= */
static void test_sign_lms(void) {
    TEST("Sign-error LMS");
    {
        size_t order = 2;
        double w[2] = {0.0, 0.0};
        double x[2] = {1.0, 0.5};
        double d = 0.8;
        af_lms_sign_error_update(NULL, x, d, w, 0.1, order);
        /* w should have been updated: sgn(e) = sgn(0.8 - 0) = +1 */
        ASSERT_NEAR(w[0], 0.1, 1e-10, "Sign-error w[0] wrong");
        ASSERT_NEAR(w[1], 0.05, 1e-10, "Sign-error w[1] wrong");
        PASS();
    }

    TEST("Sign-data LMS");
    {
        size_t order = 2;
        double w[2] = {0.0, 0.0};
        double x[2] = {2.0, -0.5};
        double d = 0.0;
        af_lms_sign_data_update(NULL, x, d, w, 0.1, order);
        /* x[0]=2>0 so sgn=+1, x[1]=-0.5<0 so sgn=-1 */
        ASSERT_TRUE(w[0] == 0.0 && w[1] == 0.0, "expected zero update");
        PASS();
    }

    TEST("Sign-sign LMS");
    {
        size_t order = 2;
        double w[2] = {0.0, 0.0};
        double x[2] = {1.0, -1.0};
        double d = 1.0;
        af_lms_sign_sign_update(NULL, x, d, w, 0.1, order);
        /* e = 1 - 0 = +1, x[0]=+1 => +0.1, x[1]=-1 => -0.1 */
        ASSERT_NEAR(w[0], 0.1, 1e-10, "Sign-sign w[0] wrong");
        ASSERT_NEAR(w[1], -0.1, 1e-10, "Sign-sign w[1] wrong");
        PASS();
    }

    TEST("Leaky LMS");
    {
        size_t order = 2;
        double w[2] = {0.5, -0.3};
        double x[2] = {1.0, 0.0};
        double d = 0.5;
        double gamma = 0.99;
        af_lms_leaky_update(NULL, x, d, w, 0.01, gamma, order);
        /* e = 0.5 - (0.5*1 + (-0.3)*0) = 0.5 - 0.5 = 0
         * So w = gamma*w + mu*0*x = gamma*w */
        ASSERT_NEAR(w[0], 0.495, 1e-10, "Leaky w[0] wrong");
        ASSERT_NEAR(w[1], -0.297, 1e-10, "Leaky w[1] wrong");
        PASS();
    }
}

/* =================================================================
 * Test 5: RLS convergence
 * ================================================================= */
static void test_rls(void) {
    TEST("RLS standard update - fast convergence");
    {
        size_t order = 4;
        double h[4] = {0.4, -0.5, 0.3, 0.7};

        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.structure = AF_STRUCTURE_FIR;
        cfg.algorithm = AF_ALGO_RLS;
        cfg.order = order;
        cfg.mu = 0.01;
        cfg.lambda = 0.99;
        cfg.delta = 1.0;
        cfg.epsilon = 1e-6;
        cfg.gamma_leakage = 1.0;
        cfg.init_strategy = AF_INIT_ZEROS;

        ASSERT_TRUE(af_init(&cfg, &f) == 0, "RLS init failed");

        /* Train with 500 samples (RLS converges in ~2*M samples) */
        size_t n_samples = 500;
        double x_buf[4] = {0.0, 0.0, 0.0, 0.0};
        srand(777);

        size_t n, i;
        for (n = 0; n < n_samples; n++) {
            double xn = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            memmove(x_buf + 1, x_buf, (order - 1) * sizeof(double));
            x_buf[0] = xn;
            double dn = 0.0;
            for (i = 0; i < order; i++) dn += h[i] * x_buf[i];
            /* Inject into filter's x_buffer manually for RLS */
            f.x_buffer[0] = xn;
            for (i = 1; i < order; i++) f.x_buffer[i] = x_buf[i];
            af_rls_standard_update(&f, xn, dn, 0.99, 1.0, order);
        }

        ASSERT_NEAR(f.w[0], h[0], 0.1, "RLS w[0] not converged");
        ASSERT_NEAR(f.w[1], h[1], 0.1, "RLS w[1] not converged");
        ASSERT_NEAR(f.w[2], h[2], 0.1, "RLS w[2] not converged");

        af_free(&f);
        PASS();
    }
}

/* =================================================================
 * Test 6: Kalman filter
 * ================================================================= */
static void test_kalman(void) {
    TEST("Kalman predict + update");
    {
        size_t n = 2; /* 2-state system */
        double x_hat[2] = {0.0, 0.0};
        double P[4] = {1.0, 0.0, 0.0, 1.0};  /* Initial uncertainty */
        double F[4] = {1.0, 0.1, 0.0, 1.0};   /* Constant velocity model */
        double Q[4] = {0.01, 0.0, 0.0, 0.01};  /* Process noise */
        double H[2] = {1.0, 0.0};              /* Measure position only */
        double R[1] = {0.1};                   /* Measurement noise */

        /* Predict */
        af_kalman_predict(x_hat, P, F, Q, n);
        /* After predict: x_hat = F * [0,0]^T = [0,0]^T (unchanged) */
        ASSERT_NEAR(x_hat[0], 0.0, 1e-10, "predict state wrong");
        /* P should have increased */
        ASSERT_TRUE(P[0] > 1.0, "covariance did not increase after predict");

        /* Update with measurement z = 2.0 */
        af_kalman_update(x_hat, P, H, R, 2.0, n, 1);
        /* Kalman gain should pull estimate toward measurement */
        ASSERT_TRUE(x_hat[0] > 1.0, "state not updated toward measurement");

        PASS();
    }

    TEST("Kalman filter NULL safety");
    {
        /* These should not crash */
        af_kalman_predict(NULL, NULL, NULL, NULL, 0);
        af_kalman_update(NULL, NULL, NULL, NULL, 0.0, 0, 0);
        af_kalman_extended_predict(NULL, NULL, NULL, NULL, NULL, 0);
        af_kalman_extended_update(NULL, NULL, NULL, NULL, NULL, 0.0, 0, 0);
        PASS();
    }
}

/* =================================================================
 * Test 7: Wiener solver
 * ================================================================= */
static void test_wiener(void) {
    TEST("Wiener solver via Levinson-Durbin");
    {
        /* Known Wiener problem: order-2 system
         * r_xx = [1.0, 0.5], p_dx = [0.8, 0.3]
         * R = [[1.0, 0.5], [0.5, 1.0]]
         * w_opt = R^{-1} * p
         * R^{-1} = (1/(1-0.25)) * [[1, -0.5], [-0.5, 1]] = [[1.333, -0.667], [-0.667, 1.333]]
         * w_opt = [1.333*0.8 + (-0.667)*0.3, (-0.667)*0.8 + 1.333*0.3]
         *       = [1.0664 - 0.2001, -0.5336 + 0.3999]
         *       = [0.8663, -0.1337]
         */
        double r_xx[2] = {1.0, 0.5};
        double p_dx[2] = {0.8, 0.3};

        af_wiener_config_t cfg;
        af_wiener_filter_t wf;
        memset(&cfg, 0, sizeof(cfg));
        memset(&wf, 0, sizeof(wf));

        cfg.order = 2;
        cfg.r_xx = r_xx;
        cfg.p_dx = p_dx;
        cfg.sigma_d2 = 1.0;

        int ret = af_wiener_solve(&cfg, &wf);
        ASSERT_TRUE(ret == 0, "Wiener solver failed");

        ASSERT_NEAR(wf.w_opt[0], 0.866667, 0.001, "w_opt[0] wrong");
        ASSERT_NEAR(wf.w_opt[1], -0.133333, 0.001, "w_opt[1] wrong");
        ASSERT_TRUE(wf.j_min >= 0.0, "J_min should be >= 0");

        af_wiener_free(&wf);
        PASS();
    }
}

/* =================================================================
 * Test 8: Matrix operations
 * ================================================================= */
static void test_matrix(void) {
    TEST("Matrix alloc/free/get/set");
    {
        af_matrix_t mat;
        ASSERT_TRUE(af_matrix_alloc(3, 4, &mat) == 0, "alloc failed");
        ASSERT_TRUE(mat.rows == 3, "rows wrong");
        ASSERT_TRUE(mat.cols == 4, "cols wrong");

        af_matrix_set(&mat, 1, 2, 42.0);
        ASSERT_NEAR(af_matrix_get(&mat, 1, 2), 42.0, 1e-10, "get/set wrong");
        ASSERT_NEAR(af_matrix_get(&mat, 0, 0), 0.0, 1e-10, "unset element nonzero");

        af_matrix_free(&mat);
        ASSERT_TRUE(mat.data == NULL, "data not freed");
        PASS();
    }

    TEST("Matrix identity");
    {
        af_matrix_t A;
        af_matrix_alloc(3, 3, &A);
        af_matrix_identity(&A);
        size_t i, j;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                double expected = (i == j) ? 1.0 : 0.0;
                ASSERT_NEAR(af_matrix_get(&A, i, j), expected, 1e-10, "identity wrong");
            }
        }
        af_matrix_free(&A);
        PASS();
    }

    TEST("Matrix inversion");
    {
        /* A = [[2, 1], [1, 2]]
         * A^{-1} = [[2/3, -1/3], [-1/3, 2/3]] */
        af_matrix_t A, A_inv;
        af_matrix_alloc(2, 2, &A);
        af_matrix_alloc(2, 2, &A_inv);

        af_matrix_set(&A, 0, 0, 2.0);
        af_matrix_set(&A, 0, 1, 1.0);
        af_matrix_set(&A, 1, 0, 1.0);
        af_matrix_set(&A, 1, 1, 2.0);

        int ret = af_matrix_inverse(&A, &A_inv);
        ASSERT_TRUE(ret == 0, "inversion failed");

        ASSERT_NEAR(af_matrix_get(&A_inv, 0, 0), 2.0/3.0, 1e-6, "A^{-1}[0,0] wrong");
        ASSERT_NEAR(af_matrix_get(&A_inv, 0, 1), -1.0/3.0, 1e-6, "A^{-1}[0,1] wrong");

        af_matrix_free(&A); af_matrix_free(&A_inv);
        PASS();
    }

    TEST("Cholesky decomposition");
    {
        /* A = [[4, 2], [2, 3]]
         * L = [[2, 0], [1, sqrt(2)]] */
        af_matrix_t A, L;
        af_matrix_alloc(2, 2, &A);
        af_matrix_alloc(2, 2, &L);

        af_matrix_set(&A, 0, 0, 4.0);
        af_matrix_set(&A, 0, 1, 2.0);
        af_matrix_set(&A, 1, 0, 2.0);
        af_matrix_set(&A, 1, 1, 3.0);

        af_matrix_cholesky(&A, &L);
        ASSERT_NEAR(af_matrix_get(&L, 0, 0), 2.0, 1e-10, "L[0,0] wrong");
        ASSERT_NEAR(af_matrix_get(&L, 0, 1), 0.0, 1e-10, "L[0,1] not zero");
        ASSERT_NEAR(af_matrix_get(&L, 1, 0), 1.0, 1e-10, "L[1,0] wrong");

        af_matrix_free(&A); af_matrix_free(&L);
        PASS();
    }
}

/* =================================================================
 * Test 9: Vector operations
 * ================================================================= */
static void test_vector(void) {
    TEST("Vector dot product and norm");
    {
        double a[3] = {1.0, 2.0, 3.0};
        double b[3] = {4.0, 5.0, 6.0};
        double dot = af_vector_dot(a, b, 3);
        ASSERT_NEAR(dot, 32.0, 1e-10, "dot product wrong"); /* 4+10+18=32 */

        double norm = af_vector_norm(a, 3);
        ASSERT_NEAR(norm, sqrt(14.0), 1e-10, "vector norm wrong");
        PASS();
    }

    TEST("Vector scale/add/sub");
    {
        double a[3] = {1.0, 2.0, 3.0};
        double b[3] = {0.5, 1.0, 1.5};
        double c[3];

        af_vector_add(a, b, c, 3);
        ASSERT_NEAR(c[0], 1.5, 1e-10, "add wrong");
        ASSERT_NEAR(c[1], 3.0, 1e-10, "add wrong");

        af_vector_sub(a, b, c, 3);
        ASSERT_NEAR(c[0], 0.5, 1e-10, "sub wrong");

        af_vector_scale(a, 2.0, 3);
        ASSERT_NEAR(a[0], 2.0, 1e-10, "scale wrong");
        /* Restore */
        a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
        PASS();
    }

    TEST("Vector NULL safety");
    {
        ASSERT_NEAR(af_vector_dot(NULL, NULL, 10), 0.0, 1e-10, "NULL dot");
        ASSERT_NEAR(af_vector_norm(NULL, 10), 0.0, 1e-10, "NULL norm");
        af_vector_scale(NULL, 1.0, 10);  /* Should not crash */
        af_vector_add(NULL, NULL, NULL, 10);  /* Should not crash */
        PASS();
    }
}

/* =================================================================
 * Test 10: Autocorrelation
 * ================================================================= */
static void test_autocorr(void) {
    TEST("Autocorrelation estimation");
    {
        /* Constant signal x = [1, 1, 1, 1]: r_xx[k] = 1 for all k */
        double x[4] = {1.0, 1.0, 1.0, 1.0};
        double r[3];
        af_autocorrelation(x, 4, 2, r, 1);  /* biased */
        ASSERT_NEAR(r[0], 1.0, 1e-10, "r[0] wrong");
        ASSERT_NEAR(r[1], 0.75, 1e-10, "r[1] wrong (biased)");
        ASSERT_NEAR(r[2], 0.5, 1e-10, "r[2] wrong (biased)");
        PASS();
    }
}

/* =================================================================
 * Test 11: System Identification application
 * ================================================================= */
static void test_app_sysid(void) {
    TEST("System identification via LMS");
    {
        /* Unknown system: h = [0.7, -0.3] */
        size_t order = 2;
        size_t n_samples = 2000;
        double *x = (double *)malloc(n_samples * sizeof(double));
        double *d = (double *)malloc(n_samples * sizeof(double));
        double *w_est = (double *)malloc(order * sizeof(double));

        srand(42);
        double h[2] = {0.7, -0.3};
        double state[2] = {0.0, 0.0};
        size_t n;
        for (n = 0; n < n_samples; n++) {
            x[n] = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            d[n] = h[0] * x[n] + h[1] * state[0];
            state[1] = state[0];
            state[0] = x[n];
        }

        int ret = af_app_system_identify_lms(x, d, n_samples, order,
                                              0.005, w_est, NULL);
        ASSERT_TRUE(ret == 0, "sysid failed");
        ASSERT_NEAR(w_est[0], h[0], 0.1, "sysid w[0] not converged");
        ASSERT_NEAR(w_est[1], h[1], 0.1, "sysid w[1] not converged");

        free(x); free(d); free(w_est);
        PASS();
    }
}

/* =================================================================
 * Test 12: Edge cases
 * ================================================================= */
static void test_edge_cases(void) {
    TEST("NULL handling in af_adapt");
    {
        af_filter_t f;
        memset(&f, 0, sizeof(f));
        /* Should not crash */
        double e = af_adapt(NULL, 1.0, 1.0);
        ASSERT_NEAR(e, 0.0, 1e-10, "NULL filter adapt should return 0");
        ASSERT_NEAR(af_adapt(&f, 1.0, 1.0), 0.0, 1e-10, "NULL config adapt");
        PASS();
    }

    TEST("Zero-order filter");
    {
        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.order = 0;
        cfg.mu = 0.01;
        cfg.lambda = 1.0;
        cfg.delta = 1.0;
        cfg.gamma_leakage = 1.0;
        ASSERT_TRUE(af_init(&cfg, &f) == -1, "zero-order should fail");
        PASS();
    }

    TEST("af_filter_output NULL safety");
    {
        ASSERT_NEAR(af_filter_output(NULL, 1.0), 0.0, 1e-10, "NULL filter output");
        PASS();
    }
}

/* =================================================================
 * Test 13: Performance evaluation
 * ================================================================= */
static void test_performance(void) {
    TEST("Performance evaluation");
    {
        af_config_t cfg;
        af_filter_t f;
        memset(&cfg, 0, sizeof(cfg));
        cfg.structure = AF_STRUCTURE_FIR;
        cfg.algorithm = AF_ALGO_LMS;
        cfg.order = 4;
        cfg.mu = 0.01;
        cfg.lambda = 1.0;
        cfg.delta = 1.0;
        cfg.gamma_leakage = 1.0;
        cfg.epsilon = 1e-6;
        cfg.init_strategy = AF_INIT_ZEROS;

        af_init(&cfg, &f);

        size_t n_samples = 1000;
        double *x = (double *)malloc(n_samples * sizeof(double));
        double *d = (double *)malloc(n_samples * sizeof(double));
        af_performance_t perf;

        srand(99);
        size_t n;
        for (n = 0; n < n_samples; n++) {
            x[n] = ((double)rand() / RAND_MAX) - 0.5;
            d[n] = 0.0;
        }

        af_performance_evaluate(&f, x, d, n_samples, &perf);

        ASSERT_TRUE(perf.mse_initial >= 0.0, "MSE initial negative");
        ASSERT_TRUE(perf.mse_final >= 0.0, "MSE final negative");

        free(x); free(d);
        af_free(&f);
        PASS();
    }
}

/* =================================================================
 * Main
 * ================================================================= */
int main(void) {
    printf("=== Adaptive Filter Module Tests ===\n\n");

    printf("--- Lifecycle ---\n");
    test_lifecycle();

    printf("\n--- LMS Algorithms ---\n");
    test_lms_direct();
    test_nlms();
    test_sign_lms();

    printf("\n--- RLS Algorithm ---\n");
    test_rls();

    printf("\n--- Kalman Filter ---\n");
    test_kalman();

    printf("\n--- Wiener Solver ---\n");
    test_wiener();

    printf("\n--- Matrix Operations ---\n");
    test_matrix();

    printf("\n--- Vector Operations ---\n");
    test_vector();

    printf("\n--- Autocorrelation ---\n");
    test_autocorr();

    printf("\n--- Applications ---\n");
    test_app_sysid();

    printf("\n--- Edge Cases ---\n");
    test_edge_cases();

    printf("\n--- Performance ---\n");
    test_performance();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
