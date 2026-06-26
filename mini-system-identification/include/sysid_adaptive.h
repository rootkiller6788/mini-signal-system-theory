/**
 * @file sysid_adaptive.h
 * @brief Adaptive and recursive identification methods (L5 Algorithms, L8 Advanced)
 *
 * Knowledge coverage:
 *   L5 - Recursive Least Squares (RLS) with exponential forgetting
 *   L5 - Recursive Instrumental Variable (RIV)
 *   L5 - Recursive Prediction Error Method (RPEM)
 *   L5 - LMS (Least Mean Squares) for FIR identification (Widrow & Hoff, 1960)
 *   L5 - Normalized LMS (NLMS)
 *   L5 - Kalman filter for parameter estimation (time-varying parameters)
 *   L5 - Extended Kalman Filter (EKF) for joint state-parameter estimation
 *   L8 - Variable forgetting factor methods (Fortescue, 1981)
 *   L8 - Directional forgetting (Kulhavý, 1987)
 *
 * Recursive identification updates model parameters online as new data arrives,
 * avoiding batch recomputation. The RLS algorithm (Plackett, 1950) minimizes:
 *   V_k(θ) = Σ_{i=1}^{k} λ^{k-i} (y(i) - φ(i)ᵀθ)²
 * with forgetting factor 0 < λ ≤ 1.
 *
 * References:
 *   Ljung (1999) Ch. 11 (Recursive Identification)
 *   Ljung & Söderström (1983) "Theory and Practice of Recursive Identification"
 *   Fortescue, Kershenbaum & Ydstie (1981) "Variable forgetting factor"
 *   Widrow & Stearns (1985) "Adaptive Signal Processing"
 */

#ifndef SYSID_ADAPTIVE_H
#define SYSID_ADAPTIVE_H

#include "sysid_types.h"
#include "sysid_models.h"

/*---------------------------------------------------------------------------
 * L5: Recursive Least Squares (RLS) with forgetting factor
 *---------------------------------------------------------------------------*/

/** RLS estimator state */
typedef struct {
    int     nparam;       /**< Number of parameters */
    double *theta;        /**< Parameter estimate [nparam] */
    double *P;            /**< Covariance matrix [nparam × nparam], row-major */
    double  lambda;       /**< Forgetting factor (0 < λ ≤ 1) */
    double  delta;        /**< Initial covariance scaling */
} sysid_rls_t;

/** Initialize RLS estimator with given dimension */
int  sysid_rls_alloc(sysid_rls_t *rls, int nparam, double lambda, double delta);
void sysid_rls_free(sysid_rls_t *rls);

/** Single RLS update step (Ljung 1999, Eq 11.10-11.12) */
int  sysid_rls_step(sysid_rls_t *rls, const double *phi, double y);

/** Batch RLS: process N data points sequentially */
int  sysid_rls_batch(sysid_rls_t *rls, const double *Phi, const double *Y,
                     int nrows, int ncols);

/*---------------------------------------------------------------------------
 * L8: Variable forgetting factor (Fortescue et al., 1981)
 *     λ(k) = 1 - [1 - φᵀP(k-1)φ] ε²(k) / Σ₀
 *     Automatically reduces λ when parameter changes are detected.
 *---------------------------------------------------------------------------*/

/** RLS with variable forgetting factor */
int  sysid_rls_vff_step(sysid_rls_t *rls, const double *phi, double y,
                        double sigma0, double lambda_min, double lambda_max);

/*---------------------------------------------------------------------------
 * L8: Directional forgetting (Kulhavý, 1987; Hägglund, 1983)
 *     Forgets only in directions with new information, avoids covariance windup.
 *---------------------------------------------------------------------------*/

/** RLS step with directional forgetting */
int  sysid_rls_df_step(sysid_rls_t *rls, const double *phi, double y,
                       double lambda, double epsilon);

/*---------------------------------------------------------------------------
 * L5: Least Mean Squares (LMS) for FIR filter identification
 *     θ(k+1) = θ(k) + μ · e(k) · φ(k)
 *     Simple gradient descent with fixed step size μ.
 *     Complexity: O(n) per sample vs O(n²) for RLS.
 *---------------------------------------------------------------------------*/

/** LMS configuration */
typedef struct {
    int     n_taps;       /**< FIR filter length */
    double *w;            /**< Filter weights [n_taps] */
    double  mu;           /**< Step size (0 < μ < 2/tr(R)) */
    double  leak;         /**< Leakage factor (0 for standard LMS) */
} sysid_lms_t;

/** Initialize LMS filter */
int  sysid_lms_alloc(sysid_lms_t *lms, int n_taps, double mu, double leak);
void sysid_lms_free(sysid_lms_t *lms);

/** Single LMS update: w(k+1) = w(k) + μ · e(k) · x(k) */
int  sysid_lms_step(sysid_lms_t *lms, const double *x, double d, double *y_out);

/** Normalized LMS: μ(k) = μ₀ / (ε + ||x(k)||²)
 *  More robust to input power variations.
 */
int  sysid_nlms_step(sysid_lms_t *lms, const double *x, double d,
                     double mu0, double epsilon, double *y_out);

/** LMS batch training over N samples */
int  sysid_lms_batch(sysid_lms_t *lms, const double *X, const double *d,
                     int N, int n_taps, int epochs);

/** Sign-error LMS: w(k+1) = w(k) + μ · sign(e(k)) · x(k) (reduced complexity) */
int  sysid_sign_lms_step(sysid_lms_t *lms, const double *x, double d, double *y_out);

/*---------------------------------------------------------------------------
 * L5: Kalman filter for parameter estimation
 *     θ(k+1) = θ(k) + w(k)    (random walk parameter model)
 *     y(k)    = φ(k)ᵀ θ(k) + v(k)
 *     Kalman gain gives optimal tradeoff between prior and data.
 *---------------------------------------------------------------------------*/

/** Kalman filter for parameter tracking */
typedef struct {
    int     nparam;
    double *theta;        /**< Parameter estimate [nparam] */
    double *P;            /**< Error covariance [nparam × nparam] */
    double  R1;           /**< Parameter drift covariance (process noise) */
    double  R2;           /**< Measurement noise variance */
} sysid_kalman_t;

/** Initialize Kalman parameter estimator */
int  sysid_kalman_alloc(sysid_kalman_t *kf, int nparam, double R1, double R2);
void sysid_kalman_free(sysid_kalman_t *kf);

/** Kalman filter update for parameter estimation */
int  sysid_kalman_step(sysid_kalman_t *kf, const double *phi, double y);

/*---------------------------------------------------------------------------
 * L5/L8: Recursive PEM (RPEM) — Recursive Prediction Error Method
 *     Implements recursive Gauss-Newton for general model structures.
 *---------------------------------------------------------------------------*/

/** RPEM state for ARMAX model */
typedef struct {
    sysid_armax_t model;    /**< Current model */
    sysid_rls_t   rls;      /**< RLS for parameter update */
    double       *psi;      /**< Gradient ψ(k) = dŷ/dθ [nparam] */
    double       *phi_y;    /**< Past outputs buffer */
    double       *phi_u;    /**< Past inputs buffer */
    double       *phi_e;    /**< Past prediction errors */
    int           nparam;
} sysid_rpem_armax_t;

/** Initialize RPEM for ARMAX */
int  sysid_rpem_armax_alloc(sysid_rpem_armax_t *rpem, int na, int nb, int nc,
                             int nk, double Ts, double lambda);
void sysid_rpem_armax_free(sysid_rpem_armax_t *rpem);

/** Single RPEM update step with new I/O sample */
int  sysid_rpem_armax_step(sysid_rpem_armax_t *rpem, double u, double y);

/*---------------------------------------------------------------------------
 * L8: Extended Kalman Filter for joint state-parameter estimation
 *     Augments state vector with parameters: [x; θ]
 *     Nonlinear dynamics require linearization at each step.
 *---------------------------------------------------------------------------*/

/** Joint state-parameter estimation for first-order system identification
 *  Estimates τ (time constant) and K (gain) online from I/O data.
 */
int  sysid_ekf_1st_order(double u, double y, double Ts,
                         double *K_est, double *tau_est,
                         double *P, double Q, double R);

#endif /* SYSID_ADAPTIVE_H */
