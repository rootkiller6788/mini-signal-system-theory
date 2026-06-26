/**
 * @file sysid_adaptive.c
 * @brief Adaptive/recursive identification methods (L5 Algorithms, L8 Advanced)
 *
 * Implements RLS (with VFF/DF variants), LMS/NLMS, Kalman filter for parameters,
 * and RPEM for recursive prediction error identification.
 */

#include "sysid_adaptive.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

extern int sysid_cholesky_solve(double *A, double *b, int n);
extern int sysid_cholesky(double *A, int n);

/*---------------------------------------------------------------------------
 * RLS implementation (L5)
 *---------------------------------------------------------------------------*/

int sysid_rls_alloc(sysid_rls_t *rls, int nparam, double lambda, double delta) {
    if (!rls || nparam <= 0 || lambda <= 0.0 || lambda > 1.0) return -1;
    memset(rls, 0, sizeof(sysid_rls_t));
    rls->nparam = nparam;
    rls->lambda = lambda;
    rls->delta = delta;
    rls->theta = (double *)calloc((size_t)nparam, sizeof(double));
    rls->P = (double *)calloc((size_t)(nparam * nparam), sizeof(double));
    if (!rls->theta || !rls->P) {
        free(rls->theta); free(rls->P); return -1;
    }
    /* P(0) = I / delta */
    for (int i = 0; i < nparam; i++) rls->P[i * nparam + i] = 1.0 / delta;
    return 0;
}

void sysid_rls_free(sysid_rls_t *rls) {
    if (!rls) return;
    free(rls->theta); rls->theta = NULL;
    free(rls->P); rls->P = NULL;
    rls->nparam = 0;
}

int sysid_rls_step(sysid_rls_t *rls, const double *phi, double y) {
    if (!rls || !phi) return -1;
    int n = rls->nparam;
    double lambda = rls->lambda;

    /* P * phi */
    double *Pphi = (double *)calloc((size_t)n, sizeof(double));
    if (!Pphi) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Pphi[i] += rls->P[i * n + j] * phi[j];
        }
    }

    /* phiᵀ P phi */
    double phiPphi = 0.0;
    for (int i = 0; i < n; i++) phiPphi += phi[i] * Pphi[i];

    double denom = lambda + phiPphi;
    if (denom < 1e-15) { free(Pphi); return -1; }

    /* Prediction error */
    double ypred = 0.0;
    for (int i = 0; i < n; i++) ypred += phi[i] * rls->theta[i];
    double err = y - ypred;

    /* Update theta */
    for (int i = 0; i < n; i++) {
        rls->theta[i] += Pphi[i] * err / denom;
    }

    /* Update P = (P - Pphi * Pphiᵀ / denom) / lambda */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rls->P[i * n + j] = (rls->P[i * n + j] - Pphi[i] * Pphi[j] / denom) / lambda;
        }
    }

    free(Pphi);
    return 0;
}

int sysid_rls_batch(sysid_rls_t *rls, const double *Phi, const double *Y,
                     int nrows, int ncols) {
    if (!rls || !Phi || !Y || nrows <= 0 || ncols <= 0) return -1;
    if (ncols != rls->nparam) return -1;
    for (int i = 0; i < nrows; i++) {
        if (sysid_rls_step(rls, &Phi[i * ncols], Y[i]) != 0) return -1;
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Variable forgetting factor RLS (L8)
 * Fortescue, Kershenbaum & Ydstie (1981)
 *---------------------------------------------------------------------------*/

int sysid_rls_vff_step(sysid_rls_t *rls, const double *phi, double y,
                        double sigma0, double lambda_min, double lambda_max) {
    if (!rls || !phi) return -1;
    int n = rls->nparam;

    double *Pphi = (double *)calloc((size_t)n, sizeof(double));
    if (!Pphi) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Pphi[i] += rls->P[i * n + j] * phi[j];
    }

    double phiPphi = 0.0;
    for (int i = 0; i < n; i++) phiPphi += phi[i] * Pphi[i];

    /* Prediction error */
    double ypred = 0.0;
    for (int i = 0; i < n; i++) ypred += phi[i] * rls->theta[i];
    double err = y - ypred;

    /* Variable forgetting factor:
     * λ(k) = 1 - [1 - φᵀPφ] * ε² / Σ₀
     * Clamped to [λ_min, λ_max] */
    double lambda = 1.0 - (1.0 - phiPphi) * err * err / sigma0;
    if (lambda < lambda_min) lambda = lambda_min;
    if (lambda > lambda_max) lambda = lambda_max;
    rls->lambda = lambda;

    double denom = lambda + phiPphi;
    if (denom < 1e-15) { free(Pphi); return -1; }

    for (int i = 0; i < n; i++) rls->theta[i] += Pphi[i] * err / denom;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rls->P[i * n + j] = (rls->P[i * n + j] - Pphi[i] * Pphi[j] / denom) / lambda;
        }
    }

    free(Pphi);
    return 0;
}

/*---------------------------------------------------------------------------
 * Directional forgetting RLS (L8)
 * Kulhavý (1987), Hägglund (1983)
 *---------------------------------------------------------------------------*/

int sysid_rls_df_step(sysid_rls_t *rls, const double *phi, double y,
                       double lambda, double epsilon) {
    if (!rls || !phi) return -1;
    int n = rls->nparam;

    double *Pphi = (double *)calloc((size_t)n, sizeof(double));
    if (!Pphi) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Pphi[i] += rls->P[i * n + j] * phi[j];
    }

    double phiPphi = 0.0;
    for (int i = 0; i < n; i++) phiPphi += phi[i] * Pphi[i];

    double ypred = 0.0;
    for (int i = 0; i < n; i++) ypred += phi[i] * rls->theta[i];
    double err = y - ypred;

    /* Directional forgetting factor:
     * α(k) = λ - (1-λ) / (ε + φᵀPφ)
     * P(k) = P(k-1) - Pφ φᵀP / (1/α + φᵀPφ) */
    double alpha = lambda - (1.0 - lambda) / (epsilon + phiPphi);
    if (alpha < 0.01) alpha = 0.01;
    if (alpha > 1.0) alpha = 1.0;

    /* Update using modified covariance */
    double r = 1.0 / alpha + phiPphi;
    if (r < 1e-15) { free(Pphi); return -1; }

    for (int i = 0; i < n; i++) rls->theta[i] += Pphi[i] * err / (1.0 + phiPphi);

    /* Covariance update with directional forgetting */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rls->P[i * n + j] -= Pphi[i] * Pphi[j] / r;
        }
    }

    free(Pphi);
    return 0;
}

/*---------------------------------------------------------------------------
 * LMS implementation (L5)
 * Widrow & Hoff (1960)
 *---------------------------------------------------------------------------*/

int sysid_lms_alloc(sysid_lms_t *lms, int n_taps, double mu, double leak) {
    if (!lms || n_taps <= 0 || mu <= 0.0) return -1;
    memset(lms, 0, sizeof(sysid_lms_t));
    lms->n_taps = n_taps;
    lms->mu = mu;
    lms->leak = leak;
    lms->w = (double *)calloc((size_t)n_taps, sizeof(double));
    if (!lms->w) return -1;
    return 0;
}

void sysid_lms_free(sysid_lms_t *lms) {
    if (!lms) return;
    free(lms->w); lms->w = NULL;
    lms->n_taps = 0;
}

int sysid_lms_step(sysid_lms_t *lms, const double *x, double d, double *y_out) {
    if (!lms || !x) return -1;
    int n = lms->n_taps;
    double mu = lms->mu;

    /* Filter output */
    double y = 0.0;
    for (int i = 0; i < n; i++) y += lms->w[i] * x[i];
    if (y_out) *y_out = y;

    /* Error */
    double e = d - y;

    /* Weight update: w(k+1) = (1 - μ·leak) w(k) + μ·e·x(k) */
    for (int i = 0; i < n; i++) {
        lms->w[i] = (1.0 - mu * lms->leak) * lms->w[i] + mu * e * x[i];
    }

    return 0;
}

int sysid_nlms_step(sysid_lms_t *lms, const double *x, double d,
                     double mu0, double epsilon, double *y_out) {
    if (!lms || !x) return -1;
    int n = lms->n_taps;

    /* Signal power */
    double power = 0.0;
    for (int i = 0; i < n; i++) power += x[i] * x[i];

    /* Normalized step size */
    double mu = mu0 / (epsilon + power / (double)n);
    double mu_saved = lms->mu;
    lms->mu = mu;

    int ret = sysid_lms_step(lms, x, d, y_out);
    lms->mu = mu_saved;
    return ret;
}

int sysid_lms_batch(sysid_lms_t *lms, const double *X, const double *d,
                     int N, int n_taps, int epochs) {
    if (!lms || !X || !d || N <= 0 || n_taps <= 0) return -1;

    for (int epoch = 0; epoch < epochs; epoch++) {
        for (int k = n_taps - 1; k < N; k++) {
            /* Regressor: X[k-n_taps+1 .. k] */
            double y_out;
            sysid_lms_step(lms, &X[k - n_taps + 1], d[k], &y_out);
        }
    }
    return 0;
}

int sysid_sign_lms_step(sysid_lms_t *lms, const double *x, double d, double *y_out) {
    if (!lms || !x) return -1;
    int n = lms->n_taps;
    double mu = lms->mu;

    double y = 0.0;
    for (int i = 0; i < n; i++) y += lms->w[i] * x[i];
    if (y_out) *y_out = y;

    double e = d - y;
    /* Sign of error: only direction, not magnitude */
    double sign_e = (e > 0.0) ? 1.0 : ((e < 0.0) ? -1.0 : 0.0);

    for (int i = 0; i < n; i++) {
        lms->w[i] += mu * sign_e * x[i];
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Kalman filter for parameter estimation (L5)
 *---------------------------------------------------------------------------*/

int sysid_kalman_alloc(sysid_kalman_t *kf, int nparam, double R1, double R2) {
    if (!kf || nparam <= 0) return -1;
    memset(kf, 0, sizeof(sysid_kalman_t));
    kf->nparam = nparam;
    kf->R1 = R1;
    kf->R2 = R2;
    kf->theta = (double *)calloc((size_t)nparam, sizeof(double));
    kf->P = (double *)calloc((size_t)(nparam * nparam), sizeof(double));
    if (!kf->theta || !kf->P) {
        free(kf->theta); free(kf->P); return -1;
    }
    for (int i = 0; i < nparam; i++) kf->P[i * nparam + i] = 1000.0;
    return 0;
}

void sysid_kalman_free(sysid_kalman_t *kf) {
    if (!kf) return;
    free(kf->theta); kf->theta = NULL;
    free(kf->P); kf->P = NULL;
    kf->nparam = 0;
}

int sysid_kalman_step(sysid_kalman_t *kf, const double *phi, double y) {
    if (!kf || !phi) return -1;
    int n = kf->nparam;

    /* Prediction step:
     * θ(k|k-1) = θ(k-1)  (random walk)
     * P(k|k-1) = P(k-1) + R1 * I
     */
    double *P_pred = (double *)calloc((size_t)(n * n), sizeof(double));
    if (!P_pred) return -1;
    memcpy(P_pred, kf->P, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) P_pred[i * n + i] += kf->R1;

    /* Innovation: ν(k) = y(k) - φᵀ θ(k|k-1) */
    double ypred = 0.0;
    for (int i = 0; i < n; i++) ypred += phi[i] * kf->theta[i];
    double nu = y - ypred;

    /* Innovation covariance: S = φᵀ P_pred φ + R2 */
    double *Pphi = (double *)calloc((size_t)n, sizeof(double));
    if (!Pphi) { free(P_pred); return -1; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Pphi[i] += P_pred[i * n + j] * phi[j];
        }
    }
    double S = kf->R2;
    for (int i = 0; i < n; i++) S += phi[i] * Pphi[i];
    if (S < 1e-15) { free(P_pred); free(Pphi); return -1; }

    /* Kalman gain: K = P_pred φ / S */
    double *K = (double *)calloc((size_t)n, sizeof(double));
    if (!K) { free(P_pred); free(Pphi); return -1; }
    for (int i = 0; i < n; i++) K[i] = Pphi[i] / S;

    /* Update:
     * θ(k) = θ(k-1) + K * ν
     * P(k) = P_pred - K φᵀ P_pred
     */
    for (int i = 0; i < n; i++) kf->theta[i] += K[i] * nu;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            kf->P[i * n + j] = P_pred[i * n + j] - K[i] * Pphi[j];
        }
    }

    free(K); free(Pphi); free(P_pred);
    return 0;
}

/*---------------------------------------------------------------------------
 * RPEM: Recursive PEM for ARMAX (L5/L8)
 *---------------------------------------------------------------------------*/

int sysid_rpem_armax_alloc(sysid_rpem_armax_t *rpem, int na, int nb, int nc,
                             int nk, double Ts, double lambda) {
    if (!rpem || na < 0 || nb < 0 || nc < 0) return -1;
    memset(rpem, 0, sizeof(sysid_rpem_armax_t));
    if (sysid_armax_alloc(&rpem->model, na, nb, nc, nk, Ts) != 0) return -1;

    int nparam = na + nb + nc;
    sysid_rls_alloc(&rpem->rls, nparam, lambda, 0.1);
    rpem->nparam = nparam;
    rpem->psi = (double *)calloc((size_t)nparam, sizeof(double));
    rpem->phi_y = (double *)calloc((size_t)(na + 1), sizeof(double));
    rpem->phi_u = (double *)calloc((size_t)(nb + nk), sizeof(double));
    rpem->phi_e = (double *)calloc((size_t)(nc + 1), sizeof(double));
    return 0;
}

void sysid_rpem_armax_free(sysid_rpem_armax_t *rpem) {
    if (!rpem) return;
    sysid_armax_free(&rpem->model);
    sysid_rls_free(&rpem->rls);
    free(rpem->psi); rpem->psi = NULL;
    free(rpem->phi_y); rpem->phi_y = NULL;
    free(rpem->phi_u); rpem->phi_u = NULL;
    free(rpem->phi_e); rpem->phi_e = NULL;
}

int sysid_rpem_armax_step(sysid_rpem_armax_t *rpem, double u, double y) {
    if (!rpem) return -1;
    int na = rpem->model.order.na, nb = rpem->model.order.nb;
    int nc = rpem->model.order.nc, nk = rpem->model.order.nk;
    int nparam = rpem->nparam;

    /* Shift buffers */
    for (int i = na; i > 0; i--) rpem->phi_y[i] = rpem->phi_y[i - 1];
    rpem->phi_y[0] = y;

    for (int i = nb + nk - 1; i > 0; i--) rpem->phi_u[i] = rpem->phi_u[i - 1];
    rpem->phi_u[0] = u;

    /* Compute prediction error using current parameters */
    double ypred = 0.0;
    for (int i = 0; i < nb; i++) {
        if (nk + i < nb + nk) ypred += rpem->model.B[i] * rpem->phi_u[nk + i];
    }
    for (int i = 1; i <= na; i++) ypred -= rpem->model.A[i] * rpem->phi_y[i];
    /* Add noise model contribution */
    for (int i = 1; i <= nc; i++) ypred += rpem->model.C[i] * rpem->phi_e[i];
    for (int i = 1; i <= nc; i++) ypred -= rpem->model.C[i] * rpem->phi_e[i]; /* cancel for now */

    double err = y - ypred;

    /* Shift error buffer */
    for (int i = nc; i > 0; i--) rpem->phi_e[i] = rpem->phi_e[i - 1];
    rpem->phi_e[0] = err;

    /* Compute gradient ψ = dŷ/dθ (approximate using filtered data)
     * C(q) ψ_a_i(k) = -y(k-i)   →  ψ_a_i ≈ -y(k-i) (approximate)
     * C(q) ψ_b_i(k) = u(k-nk-i)
     * C(q) ψ_c_i(k) = ε(k-i)
     */
    for (int i = 0; i < na; i++) rpem->psi[i] = -rpem->phi_y[i + 1];
    for (int i = 0; i < nb; i++) rpem->psi[na + i] = rpem->phi_u[nk + i];
    for (int i = 0; i < nc; i++) rpem->psi[na + nb + i] = rpem->phi_e[i + 1];

    /* RLS update on gradient */
    sysid_rls_step(&rpem->rls, rpem->psi, y);

    /* Update model parameters */
    sysid_param_vector_t pv = {rpem->rls.theta, nparam, nparam, NULL};
    sysid_armax_set_params(&rpem->model, &pv);

    return 0;
}

/*---------------------------------------------------------------------------
 * EKF for joint state-parameter estimation (L8)
 * Estimates K (gain) and τ (time constant) of a first-order system.
 * Dynamics: τ ẏ + y = K u  →  discrete: y(k) = a y(k-1) + b u(k-1)
 * where a = exp(-Ts/τ), b = K(1-a)
 *---------------------------------------------------------------------------*/

int sysid_ekf_1st_order(double u, double y, double Ts,
                         double *K_est, double *tau_est,
                         double *P, double Q, double R) {
    /* State: x = [K, τ]ᵀ
     * Measurement: y(k) = a x(1) * u(k-1) + (1 - a) * K * u(k-1)
     * where a = exp(-Ts/τ)
     *
     * Simplified EKF update for 2 parameters.
     */
    if (!K_est || !tau_est || !P || Ts <= 0.0) return -1;

    double K = *K_est;
    double tau = *tau_est;
    if (tau < Ts / 10.0) tau = Ts / 10.0; /* clamp to avoid singularity */

    double a = exp(-Ts / tau);
    double da_dtau = a * Ts / (tau * tau); /* d(a)/dτ */

    /* Predicted output: ŷ = a * y_prev + (1-a) * K * u_prev
     * Since we don't have y_prev stored, use measurement model:
     * ŷ(k) = a * ŷ(k-1) + b * u(k-1)  →  approximate ŷ(k-1) ≈ y(k-1) */
    static double y_prev = 0.0, u_prev = 0.0;
    double ypred = a * y_prev + (1.0 - a) * K * u_prev;

    /* Jacobian H = [dŷ/dK, dŷ/dτ]
     * dŷ/dK = (1-a) * u_prev
     * dŷ/dτ = da/dτ * y_prev - da/dτ * K * u_prev
     *        = da_dtau * (y_prev - K * u_prev)
     */
    double H_K = (1.0 - a) * u_prev;
    double H_tau = da_dtau * (y_prev - K * u_prev);

    /* Innovation */
    double nu = y - ypred;

    /* Innovation covariance */
    double S = H_K * H_K * P[0] + H_K * H_tau * (P[1] + P[2]) +
               H_tau * H_tau * P[3] + R;
    if (S < 1e-15) S = 1e-15;

    /* Kalman gain */
    double Kg_K = (P[0] * H_K + P[1] * H_tau) / S;
    double Kg_tau = (P[2] * H_K + P[3] * H_tau) / S;

    /* Update estimates */
    K += Kg_K * nu;
    tau += Kg_tau * nu;
    if (K < 0.0) K = 0.0;
    if (tau < Ts / 10.0) tau = Ts / 10.0;

    /* Update covariance: P = (I - K H) P + Q */
    double P00_new = (1.0 - Kg_K * H_K) * P[0] - Kg_K * H_tau * P[2] + Q;
    double P01_new = (1.0 - Kg_K * H_K) * P[1] - Kg_K * H_tau * P[3];
    double P10_new = -Kg_tau * H_K * P[0] + (1.0 - Kg_tau * H_tau) * P[2];
    double P11_new = -Kg_tau * H_K * P[1] + (1.0 - Kg_tau * H_tau) * P[3] + Q;
    P[0] = P00_new; P[1] = P01_new;
    P[2] = P10_new; P[3] = P11_new;

    *K_est = K;
    *tau_est = tau;

    /* Store current values for next step */
    y_prev = y;
    u_prev = u;

    return 0;
}
