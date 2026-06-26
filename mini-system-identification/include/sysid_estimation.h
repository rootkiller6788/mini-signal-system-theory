/**
 * @file sysid_estimation.h
 * @brief Parameter estimation algorithms for system identification (L5 Algorithms/Methods)
 *
 * Knowledge coverage:
 *   L5 - Least Squares (LS) estimation (Gauss, 1795; Aitken, 1935)
 *   L5 - Weighted Least Squares (WLS)
 *   L5 - Prediction Error Method (PEM) with Gauss-Newton/LM optimization
 *   L5 - Instrumental Variable (IV) method (Reiersøl, 1941)
 *   L5 - Maximum Likelihood estimation (Fisher, 1922)
 *   L5 - Correlation analysis (impulse response estimation)
 *   L4 - Cramér-Rao lower bound on parameter covariance
 *   L4 - Persistence of excitation condition for LS uniqueness
 *
 * The central estimation problem (Ljung 1999, §7.2):
 *   Minimize V_N(θ) = (1/N) Σₖ₌₁ᴺ || ε(k,θ) ||²
 *   where ε(k,θ) = y(k) - ŷ(k|θ) is the prediction error.
 *
 * References:
 *   Ljung (1999) Ch. 7 (Parameter Estimation Methods)
 *   Söderström & Stoica (1989) Ch. 4 (Prediction Error Methods)
 */

#ifndef SYSID_ESTIMATION_H
#define SYSID_ESTIMATION_H

#include "sysid_types.h"
#include "sysid_models.h"

/*---------------------------------------------------------------------------
 * L5: Regression matrix construction for linear-in-parameter models
 *---------------------------------------------------------------------------*/

/**
 * Build ARX regression matrix Φ and regressand vector Y
 *
 * The ARX predictor is linear in its parameters:
 *   ŷ(k|θ) = φ(k)ᵀ θ
 * where φ(k) = [-y(k-1), ..., -y(k-na), u(k-nk), ..., u(k-nk-nb+1)]ᵀ
 *
 * This yields the overdetermined system: Y ≈ Φ θ
 *
 * @param u        Input signal [N × nu]
 * @param y        Output signal [N × ny] (SISO assumed, ny=1)
 * @param N        Number of samples
 * @param na       AR order
 * @param nb       Input order
 * @param nk       Input delay
 * @param Phi      Output regression matrix [(N-nk-nb+1) × (na+nb)], row-major
 * @param Y        Output vector [N-nk-nb+1]
 * @param nrows    Output: number of rows in Phi/Y
 * @return         0 on success, -1 on error
 */
int  sysid_build_regression_arx(const double *u, const double *y, int N,
                                int na, int nb, int nk,
                                double *Phi, double *Y, int *nrows);

/** Build FIR regression matrix (input-only) */
int  sysid_build_regression_fir(const double *u, int N, int nb, int nk,
                                double *Phi, double *Y, int *nrows);

/** Compute information matrix R_N = (1/N) Σ φ(k) φ(k)ᵀ */
int  sysid_information_matrix(const double *Phi, int nrows, int ncols,
                              double *R);

/*---------------------------------------------------------------------------
 * L5: Least Squares family (explicit normal equations + QR variants)
 *---------------------------------------------------------------------------*/

/**
 * Ordinary Least Squares: θ̂ = (ΦᵀΦ)⁻¹ Φᵀ Y
 *
 * Implements normal equations with Cholesky decomposition of ΦᵀΦ.
 * Complexity: O(N·d² + d³) for N data points and d parameters.
 *
 * @param Phi      Regression matrix [nrows × ncols], row-major
 * @param Y        Observation vector [nrows]
 * @param nrows    Number of data points
 * @param ncols    Number of parameters
 * @param theta    Output: estimated parameters [ncols]
 * @param cov      Output: parameter covariance matrix [ncols × ncols], row-major
 *                  (NULL = don't compute)
 * @param sigma2   Output: estimated noise variance σ̂² = V_N(θ̂)
 * @return         0 on success, -1 on singularity
 */
int  sysid_ls_solve(const double *Phi, const double *Y, int nrows, int ncols,
                    double *theta, double *cov, double *sigma2);

/**
 * Weighted Least Squares: θ̂ = (ΦᵀWΦ)⁻¹ ΦᵀW Y
 * W is a diagonal weight matrix (stored as vector).
 */
int  sysid_wls_solve(const double *Phi, const double *Y, const double *W,
                     int nrows, int ncols,
                     double *theta, double *cov, double *sigma2);

/**
 * Recursive Least Squares (RLS) - single-step update (Ljung 1999, §11.2)
 *
 * Given current estimate θ(k-1) and covariance P(k-1), update with new data:
 *   K(k) = P(k-1) φ(k) / [λ + φ(k)ᵀ P(k-1) φ(k)]
 *   θ(k) = θ(k-1) + K(k) [y(k) - φ(k)ᵀ θ(k-1)]
 *   P(k) = [I - K(k) φ(k)ᵀ] P(k-1) / λ
 *
 * @param theta      In/out: parameter estimate [nparam]
 * @param P          In/out: covariance matrix [nparam × nparam], row-major
 * @param phi        Regressor vector [nparam]
 * @param y          New measurement
 * @param lambda     Forgetting factor (0 < λ ≤ 1); λ=1 = no forgetting
 * @param nparam     Number of parameters
 * @return           0 on success
 */
int  sysid_rls_update(double *theta, double *P, const double *phi, double y,
                      double lambda, int nparam);

/** Initialize RLS covariance to P0 = (1/delta) * I (large initial uncertainty) */
void sysid_rls_init(int nparam, double *P, double delta);

/**
 * Total Least Squares (TLS) via SVD
 * Solves [Φ | Y] ≈ [Φ̂ | Ŷ] minimizing ||[ΔΦ | ΔY]||_F subject to Ŷ = Φ̂ θ
 */
int  sysid_tls_solve(const double *Phi, const double *Y, int nrows, int ncols,
                     double *theta, double *sigma2);

/*---------------------------------------------------------------------------
 * L5: Instrumental Variable (IV) Methods
 *
 * IV solves the inconsistency of LS with colored/correlated noise:
 *   θ̂_IV = (ZᵀΦ)⁻¹ Zᵀ Y
 * where Z is the instrument matrix, correlated with Φ but uncorrelated with noise.
 *---------------------------------------------------------------------------*/

/**
 * Generate instrument matrix from delayed inputs (IV4, Söderström & Stoica 1983)
 *   z(k) = [u(k-nk), ..., u(k-nk-nz+1)]ᵀ
 */
int  sysid_iv_instruments_from_input(const double *u, int N, int nz,
                                     double *Z, int *nrows);

/** IV estimation with pre-computed instrument matrix */
int  sysid_iv_estimate(const double *Phi, const double *Z, const double *Y,
                       int nrows, int ncols, int nz,
                       double *theta, double *cov, double *sigma2);

/** Two-stage Least Squares (2SLS) - instruments from auxiliary model */
int  sysid_2sls_estimate(const double *u, const double *y, int N,
                         int na, int nb, int nk,
                         double *theta, double *sigma2);

/** Extended IV (IV4): four-step optimal IV (Stoica & Söderström, 1983) */
int  sysid_iv4_estimate(const double *u, const double *y, int N,
                        int na, int nb, int nk,
                        double *theta, double *cov, double *sigma2,
                        int max_iter);

/*---------------------------------------------------------------------------
 * L5: Prediction Error Method (PEM) - Gauss-Newton and Levenberg-Marquardt
 *
 * PEM minimizes V_N(θ) = (1/N) Σ ε²(k,θ) via iterative optimization.
 * The Gauss-Newton search direction uses the Jacobian Ψ of the predictor.
 *---------------------------------------------------------------------------*/

/**
 * PEM Gauss-Newton iteration for ARMAX model
 *
 * Returns the loss function value V_N. Updates theta and computes covariance.
 *
 * @param model     ARMAX model (in/out: updated parameter estimates)
 * @param data      Input-output data
 * @param theta     In/out: parameter vector
 * @param cov       Output: parameter covariance (NULL = skip)
 * @param V         Output: final loss function value
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance on relative parameter change
 * @return          Number of iterations used, -1 on error
 */
int  sysid_pem_armax(sysid_armax_t *model, const sysid_data_t *data,
                     sysid_param_vector_t *theta, double *cov,
                     double *V, int max_iter, double tol);

/**
 * PEM Gauss-Newton for OE model (simpler: predictor independent of y)
 */
int  sysid_pem_oe(sysid_oe_t *model, const sysid_data_t *data,
                  sysid_param_vector_t *theta, double *cov,
                  double *V, int max_iter, double tol);

/**
 * PEM for Box-Jenkins model (full generality)
 */
int  sysid_pem_bj(sysid_bj_t *model, const sysid_data_t *data,
                  sysid_param_vector_t *theta, double *cov,
                  double *V, int max_iter, double tol);

/**
 * Levenberg-Marquardt damped Gauss-Newton for ill-conditioned problems
 * (d² + μI) Δθ = -g, where μ is adaptively tuned.
 */
int  sysid_pem_lm_armax(sysid_armax_t *model, const sysid_data_t *data,
                        sysid_param_vector_t *theta, double *cov,
                        double *V, int max_iter, double tol,
                        double lambda0, double lambda_factor);

/*---------------------------------------------------------------------------
 * L5: Spectral / Correlation analysis (non-parametric methods)
 *---------------------------------------------------------------------------*/

/** Estimate impulse response via cross-correlation (ETF: Empirical Transfer Function)
 *  ĝ(k) = R_{yu}(k) / R_{uu}(0)  for white-ish input
 */
int  sysid_corr_impulse(const double *u, const double *y, int N,
                        double *impulse, int n_impulse);

/** Estimate frequency response via ETFE (Empirical Transfer Function Estimate)
 *  Ĝ(ω) = Y(ω) / U(ω)
 */
int  sysid_etfe(const double *u, const double *y, int N,
                const double *omega, int Nf, double complex *G);

/** Smoothed spectral estimate (Blackman-Tukey)
 *  Ĝ(e^{jω}) = Φ_{yu}(ω) / Φ_{uu}(ω) with windowed correlation
 */
int  sysid_blackman_tukey(const double *u, const double *y, int N,
                          const double *omega, int Nf, double complex *G,
                          int window_len);

/**
 * Frequency-domain Least Squares fitting of transfer function to FRF data
 * (Levy's method, 1959; Sanathanan-Koerner iteration)
 */
int  sysid_freq_ls_fit(const sysid_freq_data_t *fd,
                       int num_order, int den_order,
                       double *num, double *den, double *sigma2);

/*---------------------------------------------------------------------------
 * L4: Cramér-Rao Lower Bound (CRLB) — Fisher Information Matrix
 *
 * For an unbiased estimator θ̂ of parameter θ:
 *   cov(θ̂) ≥ I(θ)⁻¹   (matrix inequality, positive semidefinite sense)
 *
 * The Fisher Information Matrix for the prediction error framework:
 *   I(θ) = (1/σ²) Σₖ ψ(k,θ) ψ(k,θ)ᵀ
 * where ψ(k,θ) = ∂ŷ(k|θ)/∂θ is the predictor sensitivity (gradient).
 *
 * Reference: Ljung (1999) §9.3, Appendix 9A; Kay (1993) Ch. 3
 *---------------------------------------------------------------------------*/

/**
 * Compute Cramér-Rao Lower Bound from sensitivity matrix and noise variance
 *
 * Given the sensitivity (Jacobian) matrix Ψ = [ψ(1) ... ψ(N)]ᵀ [N × d]
 * and noise variance σ², computes:
 *   I(θ) = (1/σ²) Ψᵀ Ψ           (Fisher Information Matrix)
 *   crlb[i] = (I⁻¹)_{ii}          (diagonal entries = minimum variance bound)
 *   efficiency[i] = crlb[i] / cov[ii]   (0 ≤ η ≤ 1; η=1 → efficient estimator)
 *
 * @param Psi       Sensitivity matrix [N × d], row-major (∂ŷ/∂θ)
 * @param N         Number of data points
 * @param d         Number of parameters
 * @param sigma2    Noise variance σ² (> 0)
 * @param crlb      Output: CRLB diagonal entries [d] (minimum achievable variance)
 * @param efficiency Output: estimator efficiency per parameter [d] (NULL = skip)
 * @return          0 on success, -1 on error (singular FI matrix, sigma2 ≤ 0)
 */
int  sysid_cramer_rao_bound(const double *Psi, int N, int d,
                             double sigma2,
                             double *crlb, double *efficiency);

#endif /* SYSID_ESTIMATION_H */
