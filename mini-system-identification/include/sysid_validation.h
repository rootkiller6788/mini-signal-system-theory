/**
 * @file sysid_validation.h
 * @brief Model validation methods (L4 Fundamental Laws)
 *
 * Knowledge coverage:
 *   L4 - Akaike Information Criterion (AIC): AIC = N ln(V_N) + 2d
 *   L4 - Corrected AIC (AICc): AICc = AIC + 2d(d+1)/(N-d-1)
 *   L4 - Bayesian Information Criterion: BIC = N ln(V_N) + d ln(N)
 *   L4 - Minimum Description Length principle (Rissanen, 1978)
 *   L4 - Akaike's Final Prediction Error: FPE = V_N (1 + d/N) / (1 - d/N)
 *   L4 - Residual whiteness test (Ljung-Box Q-test, Portmanteau test)
 *   L2 - Cross-validation principle: train/validation split
 *   L2 - Bias-variance tradeoff in model order selection
 *   L2 - Overfitting detection
 *   L4 - Likelihood ratio test (F-test) for nested model comparison
 *
 * The model selection problem (Ljung 1999, §16.4):
 *   Choose model order d that minimizes the expected prediction error.
 *   AIC asymptotically equals the Kullback-Leibler divergence.
 *   BIC is consistent: chooses true order with probability → 1 as N → ∞.
 *
 * References:
 *   Akaike, H. (1974) "A new look at the statistical model identification"
 *   Schwarz, G. (1978) "Estimating the dimension of a model"
 *   Rissanen, J. (1978) "Modeling by shortest data description"
 *   Ljung, G. & Box, G. (1978) "On a measure of lack of fit in time series models"
 */

#ifndef SYSID_VALIDATION_H
#define SYSID_VALIDATION_H

#include "sysid_types.h"
#include "sysid_models.h"

/*---------------------------------------------------------------------------
 * L4: Information criteria computation
 *---------------------------------------------------------------------------*/

/**
 * Compute AIC (Akaike Information Criterion)
 *   AIC = N * ln(V_N) + 2 * d
 * where V_N = loss function = (1/N) Σ ε²(k), d = number of parameters.
 */
double sysid_aic(double loss_value, int N_data, int nparam);

/**
 * Compute AICc (corrected AIC for small samples)
 *   AICc = AIC + 2d(d+1)/(N-d-1)
 * Use when N/d < 40 (Burnham & Anderson, 2002).
 */
double sysid_aicc(double loss_value, int N_data, int nparam);

/**
 * Compute BIC (Bayesian/Schwarz Information Criterion)
 *   BIC = N * ln(V_N) + d * ln(N)
 * Consistent model order estimator; penalty grows with log(N).
 */
double sysid_bic(double loss_value, int N_data, int nparam);

/**
 * Compute MDL (Minimum Description Length, Rissanen 1978)
 *   MDL = N * ln(V_N) + d * ln(N) + (d/2) * ln(2π) + ...
 * Simplified form commonly used in practice.
 */
double sysid_mdl(double loss_value, int N_data, int nparam);

/**
 * Compute FPE (Final Prediction Error, Akaike 1969)
 *   FPE = V_N * (1 + d/N) / (1 - d/N)
 */
double sysid_fpe(double loss_value, int N_data, int nparam);

/** Fill all ICs in a fit result structure */
void sysid_fill_ics(sysid_fit_result_t *result, double loss_value);

/*---------------------------------------------------------------------------
 * L4: Residual analysis (whiteness test, cross-correlation test)
 *---------------------------------------------------------------------------*/

/**
 * Compute residuals ε(k) = y(k) - ŷ(k|θ)
 */
int  sysid_compute_residuals(const double *y, const double *ypred, int N,
                             double *residuals);

/**
 * Residual autocorrelation: r̂_εε(τ) = (1/N) Σ ε(k) ε(k-τ) for τ ≥ 0
 * Returns r[0..max_lag]. r[0] is the variance.
 */
int  sysid_residual_autocorr(const double *residuals, int N,
                             double *r, int max_lag);

/**
 * Ljung-Box Q-statistic (Portmanteau lack-of-fit test):
 *   Q = N(N+2) Σ_{k=1}^{h} r̂²(k) / (N-k)
 * Under H₀ (white noise), Q ~ χ²(h) asymptotically.
 * Returns the Q statistic.
 */
double sysid_ljung_box_q(const double *residuals, int N, int h);

/**
 * Cross-correlation test between residuals and input:
 *   r̂_{εu}(τ) = (1/N) Σ ε(k) u(k-τ)
 * Independence implies r̂_{εu}(τ) ≈ 0 for all τ.
 * Returns 1 if test passes at 95% confidence, 0 otherwise.
 */
int  sysid_crosscorr_test(const double *residuals, const double *u,
                          int N, int max_lag, double confidence);

/*---------------------------------------------------------------------------
 * L4/L2: Fit metrics and model comparison
 *---------------------------------------------------------------------------*/

/**
 * Normalized Root Mean Square Error (NRMSE) / FIT percentage
 *   NRMSE = 100 * (1 - ||y - ŷ|| / ||y - mean(y)||)
 * Returns negative value if model is worse than constant mean.
 */
double sysid_nrmse_fit(const double *y, const double *ypred, int N);

/** Root Mean Square Error */
double sysid_rmse(const double *y, const double *ypred, int N);

/** Mean Absolute Error */
double sysid_mae(const double *y, const double *ypred, int N);

/** Coefficient of determination R² = 1 - SS_res / SS_tot */
double sysid_r_squared(const double *y, const double *ypred, int N);

/** Sample mean */
double sysid_mean(const double *x, int N);

/** Sample variance (unbiased estimator: 1/(N-1)) */
double sysid_var(const double *x, int N);

/** Fill all fit metrics into result struct */
void sysid_compute_fit_metrics(const double *y, const double *ypred,
                               int N, sysid_fit_result_t *result);

/*---------------------------------------------------------------------------
 * L2: Cross-validation
 *---------------------------------------------------------------------------*/

/**
 * K-fold cross-validation:
 * Split data into K folds, train on K-1, validate on 1, average.
 *
 * Uses Monte-Carlo style: repeatable random split with fixed seed for reproducibility.
 *
 * @param u          Input data [N × nu], row-major
 * @param y          Output data [N × 1]
 * @param N          Total samples
 * @param na,nb,nk   ARX model orders to test
 * @param K          Number of folds (typically 5 or 10)
 * @param avg_fit    Output: average NRMSE fit across folds
 * @param std_fit    Output: std dev of NRMSE across folds
 * @return           0 on success
 */
int  sysid_cross_validate_arx(const double *u, const double *y, int N,
                               int na, int nb, int nk, int K,
                               double *avg_fit, double *std_fit);

/**
 * Model order selection by scanning na, nb and choosing minimum AICc
 *
 * @param u,y,N      Data
 * @param nk         Fixed input delay
 * @param na_max     Maximum AR order to scan
 * @param nb_max     Maximum input order to scan
 * @param best_na    Output: best AR order
 * @param best_nb    Output: best input order
 * @param best_aicc  Output: AICc of best model
 * @return           0 on success
 */
int  sysid_order_selection_aicc(const double *u, const double *y, int N,
                                 int nk, int na_max, int nb_max,
                                 int *best_na, int *best_nb, double *best_aicc);

/**
 * F-test for nested model comparison:
 *   F = [(V₁ - V₂) / (d₂ - d₁)] / [V₂ / (N - d₂)]
 * where model 2 nests model 1. Under H₀: model 1 is sufficient, F ~ F(d₂-d₁, N-d₂).
 * Returns the F-statistic.
 */
double sysid_ftest(double loss1, double loss2, int nparam1, int nparam2, int N);

/** F-test p-value approximation using asymptotic chi-square relationship */
double sysid_ftest_pvalue(double loss1, double loss2, int nparam1, int nparam2, int N);

#endif /* SYSID_VALIDATION_H */
