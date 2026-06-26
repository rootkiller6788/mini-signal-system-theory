/**
 * @file sysid_validation.c
 * @brief Model validation implementations (L4 Fundamental Laws, L2 Core Concepts)
 */

#include "sysid_validation.h"
#include "sysid_estimation.h"
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>

/*---------------------------------------------------------------------------
 * Information criteria (L4)
 *---------------------------------------------------------------------------*/

double sysid_aic(double loss_value, int N_data, int nparam) {
    if (loss_value <= 1e-15 || N_data <= nparam) return 1e30;
    return (double)N_data * log(loss_value) + 2.0 * (double)nparam;
}

double sysid_aicc(double loss_value, int N_data, int nparam) {
    if (N_data <= nparam + 1) return 1e30;
    double aic = sysid_aic(loss_value, N_data, nparam);
    return aic + 2.0 * (double)(nparam) * (double)(nparam + 1) /
           (double)(N_data - nparam - 1);
}

double sysid_bic(double loss_value, int N_data, int nparam) {
    if (loss_value <= 1e-15 || N_data <= nparam) return 1e30;
    return (double)N_data * log(loss_value) + (double)nparam * log((double)N_data);
}

double sysid_mdl(double loss_value, int N_data, int nparam) {
    if (loss_value <= 1e-15 || N_data <= nparam) return 1e30;
    /* MDL ≈ BIC for large N, with additional terms */
    double bic = sysid_bic(loss_value, N_data, nparam);
    /* Add Rissanen's universal prior term */
    return bic + 0.5 * (double)nparam * log(2.0 * M_PI);
}

double sysid_fpe(double loss_value, int N_data, int nparam) {
    if (N_data <= nparam) return 1e30;
    double d = (double)nparam;
    double N = (double)N_data;
    return loss_value * (1.0 + d / N) / (1.0 - d / N);
}

void sysid_fill_ics(sysid_fit_result_t *result, double loss_value) {
    if (!result) return;
    result->AIC  = sysid_aic(loss_value, result->N_data, result->nparam);
    result->AICc = sysid_aicc(loss_value, result->N_data, result->nparam);
    result->BIC  = sysid_bic(loss_value, result->N_data, result->nparam);
    result->MDL  = sysid_mdl(loss_value, result->N_data, result->nparam);
    result->FPE  = sysid_fpe(loss_value, result->N_data, result->nparam);
}

/*---------------------------------------------------------------------------
 * Residual analysis (L4)
 *---------------------------------------------------------------------------*/

int sysid_compute_residuals(const double *y, const double *ypred, int N,
                             double *residuals) {
    if (!y || !ypred || !residuals || N <= 0) return -1;
    for (int i = 0; i < N; i++) {
        residuals[i] = y[i] - ypred[i];
    }
    return 0;
}

int sysid_residual_autocorr(const double *residuals, int N,
                             double *r, int max_lag) {
    if (!residuals || !r || N <= 0 || max_lag <= 0) return -1;
    /* Compute sample variance first */
    double mean = 0.0;
    for (int i = 0; i < N; i++) mean += residuals[i];
    mean /= (double)N;

    double var = 0.0;
    for (int i = 0; i < N; i++) {
        double d = residuals[i] - mean;
        var += d * d;
    }
    var /= (double)N;
    if (var < 1e-15) {
        memset(r, 0, (size_t)(max_lag + 1) * sizeof(double));
        r[0] = 0.0;
        return 0;
    }

    for (int lag = 0; lag <= max_lag; lag++) {
        double cov = 0.0;
        for (int i = lag; i < N; i++) {
            cov += (residuals[i] - mean) * (residuals[i - lag] - mean);
        }
        r[lag] = cov / ((double)(N - lag) * var);
    }
    return 0;
}

double sysid_ljung_box_q(const double *residuals, int N, int h) {
    if (!residuals || N <= 0 || h <= 0) return -1.0;
    double *r = (double *)calloc((size_t)(h + 1), sizeof(double));
    if (!r) return -1.0;

    sysid_residual_autocorr(residuals, N, r, h);

    double Q = 0.0;
    for (int k = 1; k <= h; k++) {
        double rk2 = r[k] * r[k];
        Q += rk2 / (double)(N - k);
    }
    Q *= (double)(N * (N + 2));
    free(r);
    return Q;
}

int sysid_crosscorr_test(const double *residuals, const double *u,
                          int N, int max_lag, double confidence) {
    if (!residuals || !u || N <= 0 || max_lag <= 0) return -1;

    /* Compute cross-correlation r_{εu}(τ) */
    double mean_e = 0.0, mean_u = 0.0;
    for (int i = 0; i < N; i++) {
        mean_e += residuals[i];
        mean_u += u[i];
    }
    mean_e /= (double)N; mean_u /= (double)N;

    double var_e = 0.0, var_u = 0.0;
    for (int i = 0; i < N; i++) {
        var_e += (residuals[i] - mean_e) * (residuals[i] - mean_e);
        var_u += (u[i] - mean_u) * (u[i] - mean_u);
    }
    var_e /= (double)N; var_u /= (double)N;

    double std_prod = sqrt(var_e * var_u);
    if (std_prod < 1e-15) return 0;

    /* 95% confidence band: ±1.96 / √N */
    double z_alpha;
    if (confidence >= 0.99) z_alpha = 2.576;
    else if (confidence >= 0.95) z_alpha = 1.96;
    else z_alpha = 1.645; /* 90% */

    double threshold = z_alpha / sqrt((double)N);

    int pass = 1;
    for (int tau = -max_lag; tau <= max_lag; tau++) {
        if (tau == 0) continue; /* zero-lag cross-corr not a test */
        double cov = 0.0;
        int t_start = (tau > 0) ? tau : 0;
        int t_end = (tau > 0) ? N : (N + tau);
        int cnt = t_end - t_start;
        if (cnt <= 0) continue;
        for (int k = t_start; k < t_end; k++) {
            cov += (residuals[k] - mean_e) * (u[k - tau] - mean_u);
        }
        double r_cross = cov / ((double)cnt * std_prod);
        if (fabs(r_cross) > threshold) {
            pass = 0;
            break;
        }
    }
    return pass;
}

/*---------------------------------------------------------------------------
 * Fit metrics (L2)
 *---------------------------------------------------------------------------*/

double sysid_mean(const double *x, int N) {
    if (!x || N <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < N; i++) sum += x[i];
    return sum / (double)N;
}

double sysid_var(const double *x, int N) {
    if (!x || N <= 1) return 0.0;
    double m = sysid_mean(x, N);
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double d = x[i] - m;
        sum += d * d;
    }
    return sum / (double)(N - 1);
}

double sysid_rmse(const double *y, const double *ypred, int N) {
    if (!y || !ypred || N <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double e = y[i] - ypred[i];
        sum += e * e;
    }
    return sqrt(sum / (double)N);
}

double sysid_mae(const double *y, const double *ypred, int N) {
    if (!y || !ypred || N <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += fabs(y[i] - ypred[i]);
    }
    return sum / (double)N;
}

double sysid_nrmse_fit(const double *y, const double *ypred, int N) {
    if (!y || !ypred || N <= 0) return -1e6;
    double y_mean = sysid_mean(y, N);
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < N; i++) {
        double e = y[i] - ypred[i];
        double d = y[i] - y_mean;
        ss_res += e * e;
        ss_tot += d * d;
    }
    if (ss_tot < 1e-15) return (ss_res < 1e-15) ? 100.0 : -1e6;
    return 100.0 * (1.0 - sqrt(ss_res / ss_tot));
}

double sysid_r_squared(const double *y, const double *ypred, int N) {
    if (!y || !ypred || N <= 0) return 0.0;
    double y_mean = sysid_mean(y, N);
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < N; i++) {
        double e = y[i] - ypred[i];
        double d = y[i] - y_mean;
        ss_res += e * e;
        ss_tot += d * d;
    }
    if (ss_tot < 1e-15) return 1.0;
    return 1.0 - ss_res / ss_tot;
}

void sysid_compute_fit_metrics(const double *y, const double *ypred,
                                int N, sysid_fit_result_t *result) {
    if (!result || !y || !ypred || N <= 0) return;
    result->NRMSE_fit = sysid_nrmse_fit(y, ypred, N);
    result->RMSE      = sysid_rmse(y, ypred, N);
    result->MAE       = sysid_mae(y, ypred, N);
    result->R_squared = sysid_r_squared(y, ypred, N);

    /* Residual stats */
    double *res = (double *)calloc((size_t)N, sizeof(double));
    if (res) {
        for (int i = 0; i < N; i++) res[i] = y[i] - ypred[i];
        result->residual_mean = sysid_mean(res, N);
        result->residual_var  = sysid_var(res, N);
        result->residual_whiteness = sysid_ljung_box_q(res, N, (N < 20) ? N / 4 : 20);
        free(res);
    }
}

/*---------------------------------------------------------------------------
 * Cross-validation (L2)
 *---------------------------------------------------------------------------*/

int sysid_cross_validate_arx(const double *u, const double *y, int N,
                              int na, int nb, int nk, int K,
                              double *avg_fit, double *std_fit) {
    if (!u || !y || !avg_fit || !std_fit || N <= 0 || K <= 1) return -1;

    double *fits = (double *)calloc((size_t)K, sizeof(double));
    if (!fits) return -1;

    int fold_size = N / K;

    for (int fold = 0; fold < K; fold++) {
        int val_start = fold * fold_size;
        int val_end = (fold == K - 1) ? N : ((fold + 1) * fold_size);

        /* Train on all data except validation fold.
         * For simplicity, concatenate training segments. */
        int N_train = N - (val_end - val_start);
        double *u_train = (double *)calloc((size_t)N_train, sizeof(double));
        double *y_train = (double *)calloc((size_t)N_train, sizeof(double));
        double *u_val = (double *)calloc((size_t)(val_end - val_start), sizeof(double));
        double *y_val = (double *)calloc((size_t)(val_end - val_start), sizeof(double));
        if (!u_train || !y_train || !u_val || !y_val) {
            free(u_train); free(y_train); free(u_val); free(y_val);
            free(fits); return -1;
        }

        int train_idx = 0, val_idx = 0;
        for (int i = 0; i < N; i++) {
            if (i >= val_start && i < val_end) {
                u_val[val_idx] = u[i];
                y_val[val_idx] = y[i];
                val_idx++;
            } else {
                u_train[train_idx] = u[i];
                y_train[train_idx] = y[i];
                train_idx++;
            }
        }

        /* Train ARX model on training data */
        int nparam = na + nb;
        int max_nrows = N_train - nk - nb + 1;
        double *Phi = (double *)calloc((size_t)(max_nrows * nparam), sizeof(double));
        double *Y_vec = (double *)calloc((size_t)max_nrows, sizeof(double));
        int nrows;
        sysid_build_regression_arx(u_train, y_train, N_train, na, nb, nk,
                                    Phi, Y_vec, &nrows);

        double *theta = (double *)calloc((size_t)nparam, sizeof(double));
        if (sysid_ls_solve(Phi, Y_vec, nrows, nparam, theta, NULL, NULL) == 0) {
            /* Validate */
            sysid_arx_t model;
            sysid_arx_alloc(&model, na, nb, nk, 1.0);
            sysid_param_vector_t pv = {theta, nparam, nparam, NULL};
            sysid_arx_set_params(&model, &pv);

            double *ypred = (double *)calloc((size_t)(val_end - val_start), sizeof(double));
            sysid_arx_predict(&model, u_val, y_val, ypred, val_end - val_start);
            fits[fold] = sysid_nrmse_fit(y_val, ypred, val_end - val_start);
            free(ypred);
            sysid_arx_free(&model);
        } else {
            fits[fold] = -1e6;
        }

        free(theta); free(Phi); free(Y_vec);
        free(u_train); free(y_train); free(u_val); free(y_val);
    }

    /* Compute mean and std of fits */
    double sum = 0.0, sum2 = 0.0;
    int valid = 0;
    for (int i = 0; i < K; i++) {
        if (fits[i] > -1e5) {
            sum += fits[i];
            sum2 += fits[i] * fits[i];
            valid++;
        }
    }
    *avg_fit = (valid > 0) ? (sum / (double)valid) : -1e6;
    *std_fit = (valid > 1) ? sqrt((sum2 - sum * sum / (double)valid) / (double)(valid - 1)) : 0.0;

    free(fits);
    return 0;
}

int sysid_order_selection_aicc(const double *u, const double *y, int N,
                                 int nk, int na_max, int nb_max,
                                 int *best_na, int *best_nb, double *best_aicc) {
    if (!u || !y || !best_na || !best_nb || !best_aicc) return -1;
    *best_aicc = 1e30;
    *best_na = 1; *best_nb = 1;

    for (int na = 1; na <= na_max; na++) {
        for (int nb = 1; nb <= nb_max; nb++) {
            int nparam = na + nb;
            int max_nrows = N - nk - nb + 1;
            if (max_nrows <= nparam) continue;

            double *Phi = (double *)calloc((size_t)(max_nrows * nparam), sizeof(double));
            double *Y_vec = (double *)calloc((size_t)max_nrows, sizeof(double));
            int nrows;
            sysid_build_regression_arx(u, y, N, na, nb, nk, Phi, Y_vec, &nrows);

            double *theta = (double *)calloc((size_t)nparam, sizeof(double));
            double sigma2;
            if (sysid_ls_solve(Phi, Y_vec, nrows, nparam, theta, NULL, &sigma2) == 0) {
                double aicc = sysid_aicc(sigma2, nrows, nparam);
                if (aicc < *best_aicc) {
                    *best_aicc = aicc;
                    *best_na = na;
                    *best_nb = nb;
                }
            }
            free(theta); free(Phi); free(Y_vec);
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * F-test for nested model comparison (L4)
 *---------------------------------------------------------------------------*/

double sysid_ftest(double loss1, double loss2, int nparam1, int nparam2, int N) {
    if (nparam2 <= nparam1 || N <= nparam2) return 0.0;
    if (loss2 <= 0.0) return 1e10;
    double V1 = loss1, V2 = loss2;
    double df1 = (double)(nparam2 - nparam1);
    double df2 = (double)(N - nparam2);
    if (df2 <= 0.0) return 0.0;
    /* F = [(V1 - V2)/df1] / [V2/df2] */
    return ((V1 - V2) / df1) / (V2 / df2);
}

double sysid_ftest_pvalue(double loss1, double loss2, int nparam1, int nparam2, int N) {
    double F = sysid_ftest(loss1, loss2, nparam1, nparam2, N);
    if (F <= 0.0) return 1.0;
    /* Approximate p-value using chi-square limit:
     * For large N, (nparam2-nparam1)*F ~ χ²(d1) approximately.
     * Use Wilson-Hilferty approximation for p-value.
     */
    double d1 = (double)(nparam2 - nparam1);
    (void)(N - nparam2); /* d2 reserved for full F-distribution */
    /* Use F-distribution approximation via regularized incomplete beta
     * Abramowitz & Stegun 26.6.2:
     * P(F ≤ x) = I_{d1*x/(d1*x+d2)}(d1/2, d2/2)
     *
     * Simplified: return approximate p-value using chi-square.
     */
    double chi2 = d1 * F;
    /* Chi-square survival function approximation for moderate df */
    double z = (pow(chi2 / d1, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * d1))) /
               sqrt(2.0 / (9.0 * d1));
    /* Normal CDF approximation */
    double p = 0.5 * (1.0 + erf(-z / sqrt(2.0)));
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    return p;
}
