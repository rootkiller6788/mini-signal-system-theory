/**
 * @file kalman_filter.h
 * @brief L7/L8 — Kalman Filter: Optimal State Estimation for Linear
 *        and Nonlinear Dynamic Systems
 *
 * The Kalman filter (Kalman, 1960) is the optimal minimum-variance
 * estimator for linear Gaussian state-space models. It recursively
 * fuses noisy measurements with a process model to produce state
 * estimates with quantified uncertainty.
 *
 * Knowledge Coverage:
 *   L7-A04: GPS position/velocity tracking via Kalman filter
 *   L7-A05: Sensor fusion (IMU + GPS)
 *   L8-T03: Stochastic state estimation (Bayesian filtering)
 *   L8-T04: Extended Kalman filter (nonlinear systems)
 *   L8-T05: Information filter (canonical form)
 *   L8-T06: Monte Carlo covariance analysis
 *
 * Core Theorem (Kalman, 1960):
 *   Given x_{k+1} = F x_k + w_k,  w_k ~ N(0, Q)
 *            z_k   = H x_k + v_k,  v_k ~ N(0, R)
 *   The Kalman filter produces the MMSE estimate x̂_{k|k} with
 *   covariance P_{k|k} via predict-update recursions.
 *
 * Course Mapping:
 *   Stanford EE363 (Linear Dynamical Systems), MIT 6.241 (Dynamic Sys)
 *   Berkeley EE221A (Linear Systems), CMU 16-745 (Optimal Control)
 *   ETH 227-0216 (Estimation and Control)
 *
 * References:
 *   Kalman, R.E. "A New Approach to Linear Filtering and Prediction
 *      Problems," Trans. ASME–J. Basic Eng., 82:35–45, 1960.
 *   Simon, D. "Optimal State Estimation," Wiley, 2006.
 *   Bar-Shalom, Li, Kirubarajan. "Estimation with Applications to
 *      Tracking and Navigation," Wiley, 2001.
 */

#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include "system_defs.h"
#include "state_space.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Kalman Filter State Structure ─────────────────────────────────────── */

/** Complete Kalman filter state.
 *    x_{k+1} = F x_k + B u_k + w_k,  w_k ~ N(0, Q)
 *    z_k     = H x_k + v_k,          v_k ~ N(0, R)
 *  All matrices in row-major order.
 */
typedef struct {
    double *x;       /* n x 1   state estimate x̂_{k|k} */
    double *P;       /* n x n   error covariance P_{k|k} */
    double *F;       /* n x n   state transition matrix */
    double *H;       /* l x n   measurement matrix */
    double *Q;       /* n x n   process noise covariance */
    double *R;       /* l x l   measurement noise covariance */
    double *B;       /* n x m   input matrix (optional, NULL if no control) */
    double *K;       /* n x l   Kalman gain (computed internally) */
    double *S;       /* l x l   innovation covariance (computed internally) */

    /* Scratch buffers */
    double *temp_nn;  /* n x n */
    double *temp_nl;  /* n x l */
    double *temp_ll;  /* l x l */
    double *temp_ln;  /* l x n */
    double *temp_n1;  /* n x 1 */
    double *temp_l1;  /* l x 1 */

    int     n;        /* State dimension */
    int     m;        /* Input dimension */
    int     l;        /* Measurement dimension */
    int     owns_data;
} kalman_filter_t;

/* ── Lifecycle ─────────────────────────────────────────────────────────── */

kalman_filter_t kf_alloc(int n, int m, int l);
void            kf_free(kalman_filter_t *kf);

int kf_initialize(kalman_filter_t *kf,
                  const double *x0, const double *P0,
                  const double *F, const double *H,
                  const double *Q, const double *R,
                  const double *B);

/* ── Core Predict-Update Cycle ─────────────────────────────────────────── */

int kf_predict(kalman_filter_t *kf, const double *u);
int kf_update(kalman_filter_t *kf, const double *z);
int kf_step(kalman_filter_t *kf, const double *u, const double *z);

/* ── Innovation and Diagnostics ────────────────────────────────────────── */

int    kf_innovation(const kalman_filter_t *kf, const double *z, double *y);
double kf_nis(const kalman_filter_t *kf, const double *z);
double kf_log_likelihood(const kalman_filter_t *kf, const double *z);

/* ── Steady-State Kalman Filter ────────────────────────────────────────── */

int kf_steady_state(kalman_filter_t *kf, double tol, int max_iter);

/* ── Extended Kalman Filter (EKF) ──────────────────────────────────────── */

int ekf_predict(kalman_filter_t *kf, const double *u,
                void (*f)(const double *x, const double *u, double *x_next,
                          void *user_data),
                void (*jac_f)(const double *x, const double *u, double *F_out,
                              void *user_data),
                void *user_data);

int ekf_update(kalman_filter_t *kf, const double *z,
               void (*h)(const double *x, double *z_pred, void *user_data),
               void (*jac_h)(const double *x, double *H_out, void *user_data),
               void *user_data);

/* ── Rauch-Tung-Striebel Smoother ──────────────────────────────────────── */

int kf_rts_smooth(const double *x_filt, const double *P_filt,
                  const double *x_pred, const double *P_pred,
                  const double *F_arr, int N, int n,
                  double *x_smooth, double *P_smooth);

/* ── Information Filter Form ───────────────────────────────────────────── */

int kf_to_information_form(const kalman_filter_t *kf,
                            double *Y_out, double *y_out);

/* ── Performance Metrics ───────────────────────────────────────────────── */

double kf_rmse(const double *x_est, const double *x_true, int M, int n);
int    kf_nees(const double *x_est, const double *x_true,
               const double *P, int M, int n, double *nees_out);

/* ── GPS Position Tracking (L7 Application) ────────────────────────────── */

int kf_configure_gps_tracker(kalman_filter_t *kf, double dt,
                              double pos_std, double vel_std,
                              double proc_noise);
int kf_gps_update(kalman_filter_t *kf, const double *z, int z_len);

/* ── Specialized Utilities (L8 Advanced) ───────────────────────────────── */

/** Explicit 2x2 Cholesky decomposition for PD matrix A = [[a,b],[b,c]].
 *  L = [[sqrt(a), 0], [b/sqrt(a), sqrt(c - b^2/a)]]
 *  ~3x faster than general n x n Cholesky for GPS tracking loops.
 *  @return 0 on success, -1 if A is not positive definite. */
int kf_cholesky_2x2(const double *A, double *L);

/** Generate DWNA (Discrete White Noise Acceleration) process noise Q_d.
 *  For constant-velocity model with dt and acceleration PSD q.
 *  Q_d = q * [[dt^3/3*I, dt^2/2*I], [dt^2/2*I, dt*I]]
 *  @param dt    Sampling interval (seconds).
 *  @param q     Acceleration PSD (m^2/s^3 for position tracking).
 *  @param dim   Spatial dimension (1, 2, or 3).
 *  @param Q_out Output matrix 2*dim x 2*dim, row-major. */
int kf_process_noise_discrete(double dt, double q, int dim, double *Q_out);

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_FILTER_H */
