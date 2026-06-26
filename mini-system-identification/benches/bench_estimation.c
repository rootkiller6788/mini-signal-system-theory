#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_signals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double get_time_sec(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static void bench_ls_scale(void) {
    printf("--- B1: LS Throughput vs Problem Size ---\n");
    printf("  %6s  %6s  %10s\n", "N", "d", "Time[ms]");
    int sizes_d[] = {4, 8, 16, 32};
    for (int si = 0; si < 4; si++) {
        int d = sizes_d[si];
        int N = 1000;
        double *Phi = (double *)malloc((size_t)(N * d) * sizeof(double));
        double *Y   = (double *)malloc((size_t)N * sizeof(double));
        double *theta = (double *)malloc((size_t)d * sizeof(double));
        for (int i = 0; i < N * d; i++)
            Phi[i] = (double)rand() / RAND_MAX - 0.5;
        for (int i = 0; i < N; i++)
            Y[i] = (double)rand() / RAND_MAX - 0.5;
        int n_trials = 50;
        double t0 = get_time_sec();
        for (int t = 0; t < n_trials; t++)
            sysid_ls_solve(Phi, Y, N, d, theta, NULL, NULL);
        double t_ms = (get_time_sec() - t0) / n_trials * 1000.0;
        printf("  %6d  %6d  %10.3f\n", N, d, t_ms);
        free(Phi); free(Y); free(theta);
    }
}

static void bench_rls_vs_batch(void) {
    printf("\n--- B2: Recursive LS vs Batch LS ---\n");
    int d = 8, N_steps = 2000;
    double *theta_batch = (double *)malloc((size_t)d * sizeof(double));
    double *theta_rls   = (double *)malloc((size_t)d * sizeof(double));
    double *phi = (double *)malloc((size_t)d * sizeof(double));
    double *P   = (double *)calloc((size_t)(d * d), sizeof(double));
    double *Phi_batch = (double *)malloc((size_t)(N_steps * d) * sizeof(double));
    double *Y_batch   = (double *)malloc((size_t)N_steps * sizeof(double));
    for (int k = 0; k < N_steps; k++) {
        for (int j = 0; j < d; j++)
            Phi_batch[k * d + j] = (double)rand() / RAND_MAX - 0.5;
        Y_batch[k] = 0.5 * Phi_batch[k * d + 0]
                     + 0.3 * Phi_batch[k * d + 1]
                     + 0.1 * ((double)rand() / RAND_MAX - 0.5);
    }
    double t0 = get_time_sec();
    for (int rep = 0; rep < 20; rep++)
        sysid_ls_solve(Phi_batch, Y_batch, N_steps, d, theta_batch, NULL, NULL);
    double t_batch = (get_time_sec() - t0) / 20.0 * 1000.0;

    sysid_rls_init(d, P, 1000.0);
    memset(theta_rls, 0, (size_t)d * sizeof(double));
    t0 = get_time_sec();
    for (int rep = 0; rep < 10; rep++) {
        memset(theta_rls, 0, (size_t)d * sizeof(double));
        sysid_rls_init(d, P, 1000.0);
        for (int k = 0; k < N_steps; k++) {
            for (int j = 0; j < d; j++)
                phi[j] = Phi_batch[k * d + j];
            sysid_rls_update(theta_rls, P, phi, Y_batch[k], 0.995, d);
        }
    }
    double t_rls = (get_time_sec() - t0) / 10.0 * 1000.0;

    printf("  Batch LS (%d x %d): %.3f ms/solve\n", N_steps, d, t_batch);
    printf("  RLS (%d steps):     %.3f ms (%.1f us/step)\n",
           N_steps, t_rls, t_rls / N_steps * 1000.0);
    double err = 0.0;
    for (int j = 0; j < d; j++)
        err += fabs(theta_batch[j] - theta_rls[j]);
    printf("  ||theta_batch - theta_rls||_1 = %.2e\n", err);

    free(phi); free(theta_batch); free(theta_rls); free(P);
    free(Phi_batch); free(Y_batch);
}

static void bench_arx_predict(void) {
    printf("\n--- B3: ARX Prediction Throughput ---\n");
    int N = 10000;
    double *u = (double *)malloc((size_t)N * sizeof(double));
    double *y = (double *)malloc((size_t)N * sizeof(double));
    double *ypred = (double *)malloc((size_t)N * sizeof(double));
    for (int k = 0; k < N; k++) {
        u[k] = (double)rand() / RAND_MAX - 0.5;
        y[k] = u[k] + 0.1 * ((double)rand() / RAND_MAX - 0.5);
    }
    int orders[] = {2, 4, 8, 16, 32};
    printf("  %6s  %10s\n", "na=nb", "Time[us]");
    for (int oi = 0; oi < 5; oi++) {
        int na = orders[oi], nb = orders[oi];
        sysid_arx_t model;
        sysid_arx_alloc(&model, na, nb, 1, 0.001);
        model.A[0] = 1.0;
        for (int i = 1; i <= na; i++) model.A[i] = -0.5 / na;
        for (int i = 0; i < nb; i++) model.B[i] = 1.0 / nb;
        int n_trials = 200;
        double t0 = get_time_sec();
        for (int t = 0; t < n_trials; t++)
            sysid_arx_predict(&model, u, y, ypred, N);
        double t_us = (get_time_sec() - t0) / n_trials * 1e6;
        printf("  %6d  %10.1f\n", na, t_us);
        sysid_arx_free(&model);
    }
    free(u); free(y); free(ypred);
}

static void bench_prbs(void) {
    printf("\n--- B4: PRBS Generation Throughput ---\n");
    int N = 100000;
    int taps_10[] = {10, 7};
    sysid_prbs_t prbs;
    sysid_prbs_init(&prbs, 10, taps_10, 2, 1);
    int n_trials = 100;
    double t0 = get_time_sec();
    for (int rep = 0; rep < n_trials; rep++) {
        for (int k = 0; k < N; k++)
            sysid_prbs_next(&prbs);
    }
    double t_us = (get_time_sec() - t0) / n_trials * 1e6;
    double msps = N / (t_us + 1e-30);
    printf("  Order %d: %.1f us  (%.1f Msamples/s)\n", 10, t_us, msps);
}

static void bench_model_conversion(void) {
    printf("\n--- B5: ARX <-> State-Space Conversion ---\n");
    sysid_arx_t arx;
    sysid_arx_alloc(&arx, 4, 3, 1, 0.001);
    arx.A[1] = -1.5; arx.A[2] = 0.9; arx.A[3] = -0.3; arx.A[4] = 0.1;
    arx.B[0] = 0.05; arx.B[1] = 0.10; arx.B[2] = 0.02;
    sysid_ss_t ss;
    int n_trials = 5000;
    double t0 = get_time_sec();
    for (int i = 0; i < n_trials; i++) {
        sysid_arx_to_ss(&arx, &ss);
        if (i < n_trials - 1) {
            free(ss.A); free(ss.B); free(ss.C);
            free(ss.D); free(ss.K); free(ss.x0);
        }
    }
    double t_arx2ss = (get_time_sec() - t0) / n_trials * 1e6;
    printf("  ARX -> SS (n=4): %.2f us/conversion\n", t_arx2ss);
    sysid_arx_free(&arx);
    free(ss.A); free(ss.B); free(ss.C);
    free(ss.D); free(ss.K); free(ss.x0);
    printf("  Conversions verified.\n");
}

int main(void) {
    printf("========================================\n");
    printf("  System Identification Benchmarks\n");
    printf("========================================\n\n");
    srand(12345);
    bench_ls_scale();
    bench_rls_vs_batch();
    bench_arx_predict();
    bench_prbs();
    bench_model_conversion();
    printf("\n===== All Benchmarks Complete =====\n");
    return 0;
}
