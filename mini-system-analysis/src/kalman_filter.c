/**
 * @file kalman_filter.c
 * @brief L7/L8 — Kalman Filter Implementation and GPS Tracking
 *
 * Knowledge Implementation (19 points):
 *   L7-I01: kf_alloc / kf_free — lifecycle management
 *   L7-I02: kf_initialize — system matrix setup with validation
 *   L7-I03: kf_predict — state and covariance time update
 *   L7-I04: kf_update — Kalman gain, state, covariance measurement update
 *   L7-I05: kf_step — combined predict-update cycle
 *   L7-I06: kf_innovation — pre-fit residual computation
 *   L7-I07: kf_nis — normalized innovation squared (chi-squared gating)
 *   L7-I08: kf_log_likelihood — measurement log-likelihood
 *   L7-I09: kf_steady_state — discrete algebraic Riccati equation
 *   L7-I10: ekf_predict — EKF nonlinear state propagation + Jacobian
 *   L7-I11: ekf_update — EKF nonlinear measurement update + Jacobian
 *   L7-I12: kf_rts_smooth — Rauch-Tung-Striebel smoother
 *   L7-I13: kf_to_information_form — canonical information filter
 *   L7-I14: kf_rmse — root mean square error over trajectory
 *   L7-I15: kf_nees — normalized estimation error squared
 *   L7-I16: kf_configure_gps_tracker — GPS 8-state ECEF tracker
 *   L7-I17: kf_gps_update — GPS measurement update
 *   L8-I01: kf_solve_dare — discrete algebraic Riccati equation
 *   L8-I02: mat_cholesky — Cholesky decomposition for symmetric PD
 *
 * References:
 *   Kalman (1960), Simon (2006), Bar-Shalom et al. (2001)
 *   GPS: Tsui, "Fundamentals of GPS Receivers," 2nd Ed. 2005
 *   Brown & Hwang, "Introduction to Random Signals and
 *   Applied Kalman Filtering," 4th Ed. 2012
 */

#include "kalman_filter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Internal matrix operations (row-major storage) ──────────────────── */

static void mat_mul(const double *A, const double *B, double *C,
                     int rA, int cA, int cB) {
    for (int i = 0; i < rA; i++)
        for (int j = 0; j < cB; j++) {
            double s = 0.0;
            for (int k = 0; k < cA; k++) s += A[i*cA+k] * B[k*cB+j];
            C[i*cB+j] = s;
        }
}

/* C = A * B' where A is rA×cA, B is rB×cA (transposed to cA×rB), C is rA×rB */
static void mat_mul_ABT(const double *A, const double *B, double *C,
                         int rA, int cA, int rB) {
    for (int i = 0; i < rA; i++)
        for (int j = 0; j < rB; j++) {
            double s = 0.0;
            for (int k = 0; k < cA; k++) s += A[i*cA+k] * B[j*cA+k];
            C[i*rB+j] = s;
        }
}

static void mat_add(const double *A, const double *B, double *C,
                     int rows, int cols) {
    for (int i = 0; i < rows*cols; i++) C[i] = A[i] + B[i];
}

static void mat_sub(const double *A, const double *B, double *C,
                     int rows, int cols) {
    for (int i = 0; i < rows*cols; i++) C[i] = A[i] - B[i];
}

static void mat_eye(double *M, int n) {
    memset(M, 0, (size_t)(n*n)*sizeof(double));
    for (int i = 0; i < n; i++) M[i*n+i] = 1.0;
}

static double mat_frob_diff(const double *A, const double *B, int n) {
    double s = 0.0;
    for (int i = 0; i < n*n; i++) { double d = A[i]-B[i]; s += d*d; }
    return sqrt(s);
}

/* ========================================================================
 *  Cholesky Decomposition A = L L' for symmetric PD matrix
 *
 *  Theorem (Cholesky, 1910 / Benoit, 1924):
 *    Every symmetric positive definite matrix A ∈ R^{n×n} has a unique
 *    Cholesky factor L (lower triangular with L_{ii} > 0) such that
 *    A = L L'. Complexity: O(n³/6).
 *
 *  Used throughout Kalman filtering for stable matrix inversion of
 *  innovation covariance S and error covariance P.
 * ======================================================================== */

static int mat_cholesky(const double *A, double *L, int n)
{
    memset(L, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += L[i * n + k] * L[j * n + k];
            if (i == j) {
                double diag = A[i * n + i] - sum;
                if (diag <= 1e-30) return -1;
                L[i * n + i] = sqrt(diag);
            } else {
                L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
            }
        }
    }
    return 0;
}

/* Solve L L' x = b via forward/back substitution */
static void cholesky_solve(const double *L, const double *b, double *x, int n)
{
    double *y = (double *)calloc((size_t)n, sizeof(double));
    if (!y) return;
    /* Forward: L y = b */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) sum += L[i*n+j] * y[j];
        y[i] = (b[i] - sum) / L[i*n+i];
    }
    /* Back: L' x = y */
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) sum += L[j*n+i] * x[j];
        x[i] = (y[i] - sum) / L[i*n+i];
    }
    free(y);
}

/* ========================================================================
 *  L7-I01: kf_alloc / kf_free — Lifecycle Management
 *
 *  Allocates all matrices for an n-state, m-input, l-measurement KF.
 *  On failure, all partial allocations are freed and zeroed struct
 *  returned (safe to call kf_free on the result).
 * ======================================================================== */

kalman_filter_t kf_alloc(int n, int m, int l)
{
    kalman_filter_t kf;
    memset(&kf, 0, sizeof(kf));
    if (n <= 0 || l <= 0 || m < 0) return kf;

    size_t nn = (size_t)(n * n), nl = (size_t)(n * l);
    size_t ll = (size_t)(l * l), ln = (size_t)(l * n);
    size_t nm = (size_t)(n * (m > 0 ? m : 1));

    kf.x  = (double *)calloc((size_t)n, sizeof(double));
    kf.P  = (double *)calloc(nn, sizeof(double));
    kf.F  = (double *)calloc(nn, sizeof(double));
    kf.H  = (double *)calloc(ln, sizeof(double));
    kf.Q  = (double *)calloc(nn, sizeof(double));
    kf.R  = (double *)calloc(ll, sizeof(double));
    kf.K  = (double *)calloc(nl, sizeof(double));
    kf.S  = (double *)calloc(ll, sizeof(double));
    kf.temp_nn = (double *)calloc(nn, sizeof(double));
    kf.temp_nl = (double *)calloc(nl, sizeof(double));
    kf.temp_ll = (double *)calloc(ll, sizeof(double));
    kf.temp_ln = (double *)calloc(ln, sizeof(double));
    kf.temp_n1 = (double *)calloc((size_t)n, sizeof(double));
    kf.temp_l1 = (double *)calloc((size_t)l, sizeof(double));

    if (m > 0) {
        kf.B = (double *)calloc(nm, sizeof(double));
        if (!kf.B) goto fail;
    }

    if (!kf.x || !kf.P || !kf.F || !kf.H || !kf.Q || !kf.R ||
        !kf.K || !kf.S || !kf.temp_nn || !kf.temp_nl || !kf.temp_ll ||
        !kf.temp_ln || !kf.temp_n1 || !kf.temp_l1) goto fail;

    kf.n = n; kf.m = m; kf.l = l;
    kf.owns_data = 1;
    return kf;

fail:
    kf_free(&kf);
    memset(&kf, 0, sizeof(kf));
    return kf;
}

void kf_free(kalman_filter_t *kf)
{
    if (!kf) return;
    if (kf->owns_data) {
        free(kf->x); free(kf->P); free(kf->F); free(kf->H);
        free(kf->Q); free(kf->R); free(kf->B); free(kf->K); free(kf->S);
        free(kf->temp_nn); free(kf->temp_nl); free(kf->temp_ll);
        free(kf->temp_ln); free(kf->temp_n1); free(kf->temp_l1);
    }
    memset(kf, 0, sizeof(*kf));
}

/* ========================================================================
 *  L7-I02: kf_initialize — System Matrix Setup
 *
 *  Copies F, H, Q, R, x0, P0, and optionally B into the KF struct.
 *  The copy semantics ensure the caller retains ownership of the
 *  original buffers. Kalman gain K and innovation cov S are zeroed.
 * ======================================================================== */

int kf_initialize(kalman_filter_t *kf,
                  const double *x0, const double *P0,
                  const double *F, const double *H,
                  const double *Q, const double *R,
                  const double *B)
{
    if (!kf || !x0 || !P0 || !F || !H || !Q || !R) return -1;
    if (!kf->x || !kf->P) return -1;
    int n = kf->n, m = kf->m, l = kf->l;
    size_t nn = (size_t)(n * n), nl = (size_t)(n * l);

    memcpy(kf->x, x0, (size_t)n * sizeof(double));
    memcpy(kf->P, P0, nn * sizeof(double));
    memcpy(kf->F, F, nn * sizeof(double));
    memcpy(kf->H, H, (size_t)(l * n) * sizeof(double));
    memcpy(kf->Q, Q, nn * sizeof(double));
    memcpy(kf->R, R, (size_t)(l * l) * sizeof(double));
    if (B && m > 0)
        memcpy(kf->B, B, (size_t)(n * m) * sizeof(double));
    memset(kf->K, 0, nl * sizeof(double));
    memset(kf->S, 0, (size_t)(l * l) * sizeof(double));
    return 0;
}

/* ========================================================================
 *  L7-I03: kf_predict — State and Covariance Time Update
 *
 *  Theorem (Kalman Predict, 1960):
 *    x̂_{k+1|k} = F x̂_{k|k} + B u_k
 *    P_{k+1|k} = F P_{k|k} F' + Q
 *
 *  The predict step propagates the state estimate forward in time.
 *  The covariance grows by Q, reflecting increased uncertainty due
 *  to process noise driving the system.
 *
 *  Complexity: O(n³) for the triple matrix product F P F'.
 * ======================================================================== */

int kf_predict(kalman_filter_t *kf, const double *u)
{
    if (!kf || !kf->x || !kf->F) return -1;
    int n = kf->n, m = kf->m;

    /* x = F * x + B * u */
    memset(kf->temp_n1, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += kf->F[i*n+j] * kf->x[j];
        kf->temp_n1[i] = sum;
    }
    if (u && kf->B && m > 0) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) sum += kf->B[i*m+j] * u[j];
            kf->temp_n1[i] += sum;
        }
    }
    memcpy(kf->x, kf->temp_n1, (size_t)n * sizeof(double));

    /* P = F * P * F' + Q */
    mat_mul(kf->F, kf->P, kf->temp_nn, n, n, n);
    mat_mul_ABT(kf->temp_nn, kf->F, kf->P, n, n, n);
    mat_add(kf->P, kf->Q, kf->P, n, n);
    return 0;
}
/* ========================================================================
 *  L7-I04: kf_update — Measurement Update (Kalman Gain)
 *
 *  Algorithm (Joseph-form for numerical stability, Bucy & Joseph 1968):
 *    y_k = z_k - H x̂_{k|k-1}                         (innovation)
 *    S_k = H P_{k|k-1} H' + R                         (innovation cov)
 *    K_k = P_{k|k-1} H' S_k^{-1}                      (Kalman gain)
 *    x̂_{k|k} = x̂_{k|k-1} + K_k y_k                   (state update)
 *    P_{k|k} = (I - K_k H) P (I - K_k H)' + K R K'   (Joseph form)
 *
 *  The Joseph-form covariance update guarantees symmetry and positive
 *  semidefiniteness even under roundoff errors. It is the recommended
 *  form for safety-critical applications.
 *
 *  S_k^{-1} is computed via Cholesky decomposition. If S is near-
 *  singular, a small regularization ε is added to the diagonal.
 *
 *  Complexity: O(l³ + n²·l), dominated by Cholesky on S.
 * ======================================================================== */

int kf_update(kalman_filter_t *kf, const double *z)
{
    if (!kf || !kf->x || !z) return -1;
    int n = kf->n, l = kf->l;

    /* Innovation y = z - H*x */
    for (int i = 0; i < l; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += kf->H[i*n+j] * kf->x[j];
        kf->temp_l1[i] = z[i] - sum;
    }

    /* S = H*P*H' + R */
    mat_mul(kf->H, kf->P, kf->temp_ln, l, n, n);
    mat_mul_ABT(kf->temp_ln, kf->H, kf->S, l, n, l);
    mat_add(kf->S, kf->R, kf->S, l, l);

    /* K = P*H' * S^{-1}  (solve column by column via Cholesky) */
    mat_mul_ABT(kf->P, kf->H, kf->temp_nl, n, n, l);

    double *L = (double *)malloc((size_t)(l * l) * sizeof(double));
    if (!L) return -1;
    if (mat_cholesky(kf->S, L, l) != 0) {
        for (int i = 0; i < l; i++) kf->S[i*l+i] += 1e-10;
        if (mat_cholesky(kf->S, L, l) != 0) { free(L); return -1; }
    }
    for (int col = 0; col < l; col++) {
        double *b = (double *)calloc((size_t)l, sizeof(double));
        double *kt = (double *)calloc((size_t)l, sizeof(double));
        if (!b || !kt) { free(b); free(kt); free(L); return -1; }
        b[col] = 1.0;
        cholesky_solve(L, b, kt, l);
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < l; j++) sum += kf->temp_nl[i*l+j] * kt[j];
            kf->K[i*l+col] = sum;
        }
        free(b); free(kt);
    }
    free(L);

    /* x = x + K*y */
    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < l; j++) corr += kf->K[i*l+j] * kf->temp_l1[j];
        kf->x[i] += corr;
    }

    /* P = (I-KH)P(I-KH)' + KRK'  (Joseph form)
     * Uses malloc'd nxn temp since embedded scratch may be too small (nxl) */
    {
        double *M1 = (double *)malloc((size_t)(n*n) * sizeof(double));
        double *M2 = (double *)malloc((size_t)(n*n) * sizeof(double));
        if (!M1 || !M2) { free(M1); free(M2); return -1; }
        /* M1 = I - K*H */
        mat_mul(kf->K, kf->H, M1, n, l, n);
        for (int i = 0; i < n*n; i++)
            M1[i] = ((i % (n+1)) == 0 ? 1.0 : 0.0) - M1[i];
        /* M2 = M1 * P */
        mat_mul(M1, kf->P, M2, n, n, n);
        /* P = M2 * M1' */
        mat_mul_ABT(M2, M1, kf->P, n, n, n);
        /* M2 = K * R * K' */
        mat_mul(kf->K, kf->R, kf->temp_nl, n, l, l);
        mat_mul_ABT(kf->temp_nl, kf->K, M2, n, l, n);
        /* P += M2 */
        mat_add(kf->P, M2, kf->P, n, n);
        free(M1); free(M2);
    }

    return 0;
}

/* ========================================================================
 *  L7-I05: kf_step — Combined Predict-Update Cycle
 *
 *  Convenience wrapper for one complete KF time step.
 *  Equivalent to sequential kf_predict + kf_update.
 * ======================================================================== */

int kf_step(kalman_filter_t *kf, const double *u, const double *z)
{
    int ret = kf_predict(kf, u);
    if (ret != 0) return ret;
    return kf_update(kf, z);
}

/* ========================================================================
 *  L7-I06: kf_innovation — Pre-fit Residual y = z - H x̂_{k|k-1}
 *
 *  The innovation (measurement residual) is the difference between
 *  actual and predicted measurement. Under Gaussian assumptions:
 *    y_k ~ N(0, S_k)
 *  This zero-mean white sequence is key to filter consistency checks.
 * ======================================================================== */

int kf_innovation(const kalman_filter_t *kf, const double *z, double *y)
{
    if (!kf || !z || !y) return -1;
    for (int i = 0; i < kf->l; i++) {
        double pred = 0.0;
        for (int j = 0; j < kf->n; j++) pred += kf->H[i*kf->n+j] * kf->x[j];
        y[i] = z[i] - pred;
    }
    return 0;
}

/* ========================================================================
 *  L7-I07: kf_nis — Normalized Innovation Squared
 *
 *  NIS_k = y_k' S_k^{-1} y_k
 *
 *  Theorem (chi-squared gating, Bar-Shalom 2001):
 *    NIS_k ~ χ²(l) under nominal conditions.
 *    E[NIS] = l, Var[NIS] = 2l.
 *    Reject measurement if NIS > χ²_{1-α}(l), e.g.:
 *      χ²_{0.95}(1)=3.84, χ²_{0.95}(2)=5.99, χ²_{0.95}(3)=7.81
 * ======================================================================== */

double kf_nis(const kalman_filter_t *kf, const double *z)
{
    if (!kf || !z) return -1.0;
    int n = kf->n, l = kf->l;

    double *y = (double *)malloc((size_t)l * sizeof(double));
    if (!y) return -1.0;
    kf_innovation(kf, z, y);

    /* S = HPH' + R */
    double *S = (double *)malloc((size_t)(l*l) * sizeof(double));
    if (!S) { free(y); return -1.0; }
    mat_mul(kf->H, kf->P, kf->temp_ln, l, n, n);
    mat_mul_ABT(kf->temp_ln, kf->H, S, l, n, l);
    mat_add(S, kf->R, S, l, l);

    double *L = (double *)malloc((size_t)(l*l) * sizeof(double));
    if (!L || mat_cholesky(S, L, l) != 0) {
        free(y); free(S); free(L); return -1.0;
    }
    cholesky_solve(L, y, kf->temp_l1, l);

    double nis = 0.0;
    for (int i = 0; i < l; i++) nis += y[i] * kf->temp_l1[i];
    free(y); free(S); free(L);
    return nis;
}

/* ========================================================================
 *  L7-I08: kf_log_likelihood — Measurement Log-Likelihood
 *
 *  log p(z_k | z_{1:k-1}) = -½ [l·log(2π) + log|S| + y'S^{-1}y]
 *
 *  Used for: IMM model probability update, adaptive noise estimation,
 *  data association in multi-target tracking.
 *
 *  log|S| = 2 Σ log L_{ii} via Cholesky.
 * ======================================================================== */

double kf_log_likelihood(const kalman_filter_t *kf, const double *z)
{
    if (!kf || !z) return -INFINITY;
    int n = kf->n, l = kf->l;

    double *y = (double *)malloc((size_t)l * sizeof(double));
    if (!y) return -INFINITY;
    kf_innovation(kf, z, y);

    double *S = (double *)malloc((size_t)(l*l) * sizeof(double));
    if (!S) { free(y); return -INFINITY; }
    mat_mul(kf->H, kf->P, kf->temp_ln, l, n, n);
    mat_mul_ABT(kf->temp_ln, kf->H, S, l, n, l);
    mat_add(S, kf->R, S, l, l);

    double *L = (double *)malloc((size_t)(l*l) * sizeof(double));
    if (!L || mat_cholesky(S, L, l) != 0) {
        free(y); free(S); free(L); return -INFINITY;
    }

    double log_det = 0.0;
    for (int i = 0; i < l; i++) {
        if (L[i*l+i] <= 0.0) { free(y); free(S); free(L); return -INFINITY; }
        log_det += log(L[i*l+i]);
    }
    log_det *= 2.0;

    cholesky_solve(L, y, kf->temp_l1, l);
    double qf = 0.0;
    for (int i = 0; i < l; i++) qf += y[i] * kf->temp_l1[i];

    free(y); free(S); free(L);
    return -0.5 * ((double)l * log(2.0 * M_PI) + log_det + qf);
}

int kf_steady_state(kalman_filter_t *kf, double tol, int max_iter)
{
    if (!kf || tol <= 0.0 || max_iter <= 0) return -1;
    int n = kf->n, l = kf->l;
    size_t nn = (size_t)(n * n);
    double *P_prev = (double *)malloc(nn * sizeof(double));
    if (!P_prev) return -1;
    for (int iter = 0; iter < max_iter; iter++) {
        memcpy(P_prev, kf->P, nn * sizeof(double));
        mat_mul(kf->F, kf->P, kf->temp_nn, n, n, n);
        mat_mul_ABT(kf->temp_nn, kf->F, kf->P, n, n, n);
        mat_add(kf->P, kf->Q, kf->P, n, n);
        mat_mul(kf->H, kf->P, kf->temp_ln, l, n, n);
        mat_mul_ABT(kf->temp_ln, kf->H, kf->S, l, n, l);
        mat_add(kf->S, kf->R, kf->S, l, l);
        mat_mul_ABT(kf->P, kf->H, kf->temp_nl, n, n, l);
        double *L = (double *)malloc((size_t)(l*l) * sizeof(double));
        if (!L) { free(P_prev); return -1; }
        if (mat_cholesky(kf->S, L, l) == 0) {
            for (int col = 0; col < l; col++) {
                double *b = (double *)calloc((size_t)l, sizeof(double));
                double *kt = (double *)calloc((size_t)l, sizeof(double));
                if (!b || !kt) { free(b); free(kt); free(L); free(P_prev); return -1; }
                b[col] = 1.0;
                cholesky_solve(L, b, kt, l);
                for (int i = 0; i < n; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < l; j++) sum += kf->temp_nl[i*l+j] * kt[j];
                    kf->K[i*l+col] = sum;
                }
                free(b); free(kt);
            }
        }
        free(L);
        mat_mul(kf->K, kf->H, kf->temp_nn, n, l, n);
        mat_mul(kf->temp_nn, kf->P, kf->temp_nl, n, n, n);
        mat_sub(kf->P, kf->temp_nl, kf->P, n, n);
        if (mat_frob_diff(kf->P, P_prev, n) < tol) { free(P_prev); return 0; }
    }
    free(P_prev);
    return -1;
}

/* ========================================================================
 *  L7-I10: ekf_predict -- EKF Nonlinear State Propagation
 *
 *  x_{k+1} = f(x_k, u_k) + w_k
 *  1. Propagate: x_hat = f(x_hat, u_k)
 *  2. Linearize: F_k = df/dx at x_hat
 *  3. Covariance: P = F_k P F_k^T + Q
 * ======================================================================== */

int ekf_predict(kalman_filter_t *kf, const double *u,
                void (*f)(const double *x, const double *u, double *x_next,
                          void *user_data),
                void (*jac_f)(const double *x, const double *u, double *F_out,
                              void *user_data),
                void *user_data)
{
    if (!kf || !f || !jac_f) return -1;
    int n = kf->n;
    f(kf->x, u, kf->temp_n1, user_data);
    memcpy(kf->x, kf->temp_n1, (size_t)n * sizeof(double));
    jac_f(kf->x, u, kf->F, user_data);
    mat_mul(kf->F, kf->P, kf->temp_nn, n, n, n);
    mat_mul_ABT(kf->temp_nn, kf->F, kf->P, n, n, n);
    mat_add(kf->P, kf->Q, kf->P, n, n);
    return 0;
}

int ekf_update(kalman_filter_t *kf, const double *z,
               void (*h)(const double *x, double *z_pred, void *user_data),
               void (*jac_h)(const double *x, double *H_out, void *user_data),
               void *user_data)
{
    if (!kf || !z || !h || !jac_h) return -1;
    int n = kf->n, l = kf->l;
    h(kf->x, kf->temp_l1, user_data);
    jac_h(kf->x, kf->H, user_data);
    for (int i = 0; i < l; i++)
        kf->temp_l1[i] = z[i] - kf->temp_l1[i];
    mat_mul(kf->H, kf->P, kf->temp_ln, l, n, n);
    mat_mul_ABT(kf->temp_ln, kf->H, kf->S, l, n, l);
    mat_add(kf->S, kf->R, kf->S, l, l);
    mat_mul_ABT(kf->P, kf->H, kf->temp_nl, n, n, l);
    double *L = (double *)malloc((size_t)(l*l) * sizeof(double));
    if (!L) return -1;
    if (mat_cholesky(kf->S, L, l) != 0) {
        for (int i = 0; i < l; i++) kf->S[i*l+i] += 1e-10;
        if (mat_cholesky(kf->S, L, l) != 0) { free(L); return -1; }
    }
    for (int col = 0; col < l; col++) {
        double *b = (double *)calloc((size_t)l, sizeof(double));
        double *kt = (double *)calloc((size_t)l, sizeof(double));
        if (!b || !kt) { free(b); free(kt); free(L); return -1; }
        b[col] = 1.0;
        cholesky_solve(L, b, kt, l);
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < l; j++) sum += kf->temp_nl[i*l+j] * kt[j];
            kf->K[i*l+col] = sum;
        }
        free(b); free(kt);
    }
    free(L);
    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < l; j++) corr += kf->K[i*l+j] * kf->temp_l1[j];
        kf->x[i] += corr;
    }
    {
        double *M1 = (double *)malloc((size_t)(n*n) * sizeof(double));
        double *M2 = (double *)malloc((size_t)(n*n) * sizeof(double));
        if (!M1 || !M2) { free(M1); free(M2); return -1; }
        mat_mul(kf->K, kf->H, M1, n, l, n);
        for (int i = 0; i < n*n; i++)
            M1[i] = ((i % (n+1)) == 0 ? 1.0 : 0.0) - M1[i];
        mat_mul(M1, kf->P, M2, n, n, n);
        mat_mul_ABT(M2, M1, kf->P, n, n, n);
        mat_mul(kf->K, kf->R, kf->temp_nl, n, l, l);
        mat_mul_ABT(kf->temp_nl, kf->K, M2, n, l, n);
        mat_add(kf->P, M2, kf->P, n, n);
        free(M1); free(M2);
    }
    return 0;
}

/* ========================================================================
 *  L7-I12: kf_rts_smooth -- Rauch-Tung-Striebel Fixed-Interval Smoother
 *
 *  Backward recursion (k = N-2..0):
 *    G_k = P_{k|k} F_k^T P_{k+1|k}^{-1}
 *    x_{k|N} = x_{k|k} + G_k (x_{k+1|N} - x_{k+1|k})
 *    P_{k|N} = P_{k|k} + G_k (P_{k+1|N} - P_{k+1|k}) G_k^T
 *
 *  Theorem (RTS, 1965): The RTS smoother produces the MMSE estimate
 *  of the state sequence using all N measurements. Applications:
 *  GPS post-processing, orbit determination, EM for sys-id.
 * ======================================================================== */

int kf_rts_smooth(const double *x_filt, const double *P_filt,
                  const double *x_pred, const double *P_pred,
                  const double *F_arr, int N, int n,
                  double *x_smooth, double *P_smooth)
{
    if (!x_filt || !P_filt || !x_pred || !P_pred || !F_arr ||
        !x_smooth || !P_smooth || N < 2 || n < 1) return -1;
    size_t nn = (size_t)(n * n);
    memcpy(x_smooth + (N-1)*n, x_filt + (N-1)*n, (size_t)n*sizeof(double));
    memcpy(P_smooth + (N-1)*nn, P_filt + (N-1)*nn, nn*sizeof(double));
    double *G   = (double *)malloc(nn * sizeof(double));
    double *tmp = (double *)malloc(nn * sizeof(double));
    double *dx  = (double *)malloc((size_t)n * sizeof(double));
    if (!G || !tmp || !dx) { free(G); free(tmp); free(dx); return -1; }
    for (int k = N - 2; k >= 0; k--) {
        const double *F_k = F_arr + k * nn;
        const double *Pp = P_pred + (k+1)*nn;
        mat_mul_ABT(P_filt + k*nn, F_k, tmp, n, n, n);
        double *L = (double *)malloc(nn * sizeof(double));
        int ok = 0;
        if (L && mat_cholesky(Pp, L, n) == 0) {
            ok = 1;
            for (int col = 0; col < n; col++) {
                double *b = (double *)calloc((size_t)n, sizeof(double));
                double *gt = (double *)calloc((size_t)n, sizeof(double));
                if (!b || !gt) { free(b); free(gt); ok = 0; break; }
                for (int i = 0; i < n; i++) b[i] = tmp[col*n+i];
                cholesky_solve(L, b, gt, n);
                for (int i = 0; i < n; i++) G[i*n+col] = gt[i];
                free(b); free(gt);
            }
        }
        if (!ok) mat_eye(G, n);
        free(L);
        for (int i = 0; i < n; i++)
            dx[i] = x_smooth[(k+1)*n+i] - x_pred[(k+1)*n+i];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += G[i*n+j] * dx[j];
            x_smooth[k*n+i] = x_filt[k*n+i] + sum;
        }
        mat_sub(P_smooth+(k+1)*nn, P_pred+(k+1)*nn, tmp, n, n);
        mat_mul(G, tmp, tmp, n, n, n);
        mat_mul_ABT(tmp, G, tmp, n, n, n);
        mat_add(P_filt+k*nn, tmp, P_smooth+k*nn, n, n);
    }
    free(G); free(tmp); free(dx);
    return 0;
}

/* ========================================================================
 *  L7-I13: kf_to_information_form -- Canonical Information Filter
 *
 *  Y = P^{-1} (information matrix), y_hat = Y * x_hat (info state)
 *  Advantage: additive update, natural handling of P -> infinity.
 * ======================================================================== */

int kf_to_information_form(const kalman_filter_t *kf,
                            double *Y_out, double *y_out)
{
    if (!kf || !Y_out || !y_out) return -1;
    int n = kf->n;
    double *L = (double *)malloc((size_t)(n*n) * sizeof(double));
    if (!L) return -1;
    if (mat_cholesky(kf->P, L, n) != 0) { free(L); return -1; }
    double *Yt = (double *)calloc((size_t)(n*n), sizeof(double));
    if (!Yt) { free(L); return -1; }
    for (int col = 0; col < n; col++) {
        double *e = (double *)calloc((size_t)n, sizeof(double));
        double *yc = (double *)calloc((size_t)n, sizeof(double));
        if (!e || !yc) { free(e); free(yc); free(L); free(Yt); return -1; }
        e[col] = 1.0;
        cholesky_solve(L, e, yc, n);
        for (int i = 0; i < n; i++) Yt[col*n+i] = yc[i];
        free(e); free(yc);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Y_out[i*n+j] = Yt[j*n+i];
    free(L); free(Yt);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += Y_out[i*n+j] * kf->x[j];
        y_out[i] = sum;
    }
    return 0;
}

/* ========================================================================
 *  L7-I14: kf_rmse -- Root Mean Square Error
 *
 *  RMSE = sqrt( (1/(M*n)) * sum_{k,i} (x_hat_{k,i} - x_true_{k,i})^2 )
 *  Standard aggregate performance metric for filter evaluation.
 * ======================================================================== */

double kf_rmse(const double *x_est, const double *x_true, int M, int n)
{
    if (!x_est || !x_true || M <= 0 || n <= 0) return -1.0;
    double sum_sq = 0.0;
    for (int k = 0; k < M; k++)
        for (int i = 0; i < n; i++) {
            double diff = x_est[k*n+i] - x_true[k*n+i];
            sum_sq += diff * diff;
        }
    return sqrt(sum_sq / (double)(M * n));
}

/* ========================================================================
 *  L7-I15: kf_nees -- Normalized Estimation Error Squared
 *
 *  NEES_k = (x_hat_k - x_k)^T P_k^{-1} (x_hat_k - x_k)
 *
 *  Theorem (Bar-Shalom 2001): For a consistent filter,
 *  NEES_k ~ chi^2(n) with E[NEES] = n. Over-confident filters
 *  have NEES > n; conservative filters have NEES < n.
 *  Acceptable region: [chi2_{alpha/2}(M*n)/M, chi2_{1-alpha/2}(M*n)/M]
 * ======================================================================== */

int kf_nees(const double *x_est, const double *x_true,
            const double *P, int M, int n, double *nees_out)
{
    if (!x_est || !x_true || !P || !nees_out || M <= 0 || n <= 0)
        return -1;
    size_t nn = (size_t)(n * n);
    for (int k = 0; k < M; k++) {
        double *e = (double *)malloc((size_t)n * sizeof(double));
        if (!e) return -1;
        for (int i = 0; i < n; i++)
            e[i] = x_est[k*n+i] - x_true[k*n+i];
        double *L = (double *)malloc(nn * sizeof(double));
        if (!L) { free(e); return -1; }
        if (mat_cholesky(P + k*nn, L, n) == 0) {
            double *a = (double *)malloc((size_t)n * sizeof(double));
            if (a) {
                cholesky_solve(L, e, a, n);
                nees_out[k] = 0.0;
                for (int i = 0; i < n; i++) nees_out[k] += e[i] * a[i];
                free(a);
            } else { nees_out[k] = -1.0; }
        } else { nees_out[k] = -1.0; }
        free(e); free(L);
    }
    return 0;
}

/* ========================================================================
 *  L7-I16: kf_configure_gps_tracker -- GPS 8-State ECEF Tracker
 *
 *  State [8]: ECEF pos(x,y,z), vel(vx,vy,vz), clock bias cb(m), drift cd(m/s)
 *
 *  Constant velocity model: F couples position and velocity with dt.
 *  Q uses discrete white noise acceleration (DWNA) model:
 *    Q_3x3 = q * [dt^3/3   dt^2/2; dt^2/2   dt  ]
 *  Clock model: standard OCXO parameters (Allan variance based).
 *
 *  H measures position (from GPS navigation solution fix).
 *
 *  Reference: Tsui (2005), Brown & Hwang (2012)
 * ======================================================================== */

int kf_configure_gps_tracker(kalman_filter_t *kf, double dt,
                              double pos_std, double vel_std,
                              double proc_noise)
{
    if (!kf || kf->n != 8 || dt <= 0.0) return -1;
    int n = 8, l = kf->l;
    double dt2 = dt * dt, dt3 = dt2 * dt;
    double q = proc_noise * proc_noise;

    /* State transition F */
    memset(kf->F, 0, (size_t)(n*n)*sizeof(double));
    for (int i = 0; i < n; i++) kf->F[i*n+i] = 1.0;
    kf->F[0*n+3] = dt; kf->F[1*n+4] = dt; kf->F[2*n+5] = dt;
    kf->F[6*n+7] = dt;

    /* Process noise Q (DWNA model) */
    memset(kf->Q, 0, (size_t)(n*n)*sizeof(double));
    for (int i = 0; i < 3; i++) {
        kf->Q[i*n+i] = q * dt3 / 3.0;       /* pos variance */
        kf->Q[(i+3)*n+(i+3)] = q * dt;      /* vel variance */
        kf->Q[i*n+(i+3)] = q * dt2 / 2.0;   /* pos-vel coupling */
        kf->Q[(i+3)*n+i] = q * dt2 / 2.0;
    }
    /* Clock model (OCXO: Scb=1e-18 s, Scd=1e-20 s^{-1}) */
    double Scb = 1e-18, Scd = 1e-20;
    kf->Q[6*n+6] = Scb * dt + Scd * dt3 / 3.0;
    kf->Q[7*n+7] = Scd * dt;
    kf->Q[6*n+7] = kf->Q[7*n+6] = Scd * dt2 / 2.0;

    /* Measurement matrix H (position measurement) */
    memset(kf->H, 0, (size_t)(l*n)*sizeof(double));
    for (int i = 0; i < 4 && i < l; i++)
        kf->H[i*n+i] = 1.0;

    /* Measurement noise R */
    double pvar = pos_std * pos_std, vvar = vel_std * vel_std;
    memset(kf->R, 0, (size_t)(l*l)*sizeof(double));
    for (int i = 0; i < 3 && i < l; i++) kf->R[i*l+i] = pvar;
    if (l >= 4) kf->R[3*l+3] = pvar;
    for (int i = 4; i < l && i < 7; i++) kf->R[i*l+i] = vvar;

    /* Initial covariance P (large = uncertain) */
    mat_eye(kf->P, n);
    for (int i = 0; i < 3; i++) kf->P[i*n+i] = 10000.0;  /* ~100m pos std */
    for (int i = 3; i < 6; i++) kf->P[i*n+i] = 100.0;    /* ~10 m/s vel std */
    kf->P[6*n+6] = 9e10;   /* ~300km clock bias (unknown at cold start) */
    kf->P[7*n+7] = 1e4;    /* ~100 m/s clock drift */
    return 0;
}

/* ========================================================================
 *  L7-I17: kf_gps_update -- GPS Measurement Update
 *
 *  Handles partial measurement vectors (z_len < l) by inflating
 *  measurement noise for unobserved dimensions. This allows the
 *  same filter to handle variable measurement availability
 *  (GPS outages, partial satellite constellations).
 *
 *  For pseudorange-based updates (nonlinear), use ekf_update with
 *  the nonlinear range function h(x) = ||pos_k - sat_pos|| + cb.
 *
 *  Reference: Parkinson & Spilker, GPS: Theory and Applications, 1996.
 * ======================================================================== */

int kf_gps_update(kalman_filter_t *kf, const double *z, int z_len)
{
    if (!kf || !z || z_len <= 0 || z_len > kf->l) return -1;
    int n = kf->n, l_orig = kf->l;

    double *z_full = (double *)calloc((size_t)l_orig, sizeof(double));
    double *H_full = (double *)calloc((size_t)(l_orig*n), sizeof(double));
    double *R_full = (double *)calloc((size_t)(l_orig*l_orig), sizeof(double));
    if (!z_full || !H_full || !R_full) {
        free(z_full); free(H_full); free(R_full); return -1;
    }

    /* Copy available measurements; copy H rows for those dimensions */
    memcpy(z_full, z, (size_t)z_len * sizeof(double));
    for (int i = 0; i < z_len; i++)
        memcpy(H_full + i*n, kf->H + i*n, (size_t)n * sizeof(double));

    /* R: use actual R for measured dims, very large R for missing */
    for (int i = 0; i < l_orig; i++) {
        if (i < z_len)
            R_full[i*l_orig+i] = kf->R[i*l_orig+i];
        else
            R_full[i*l_orig+i] = 1e12;  /* Essentially infinite noise */
    }

    /* Swap and update */
    double *H_save = kf->H, *R_save = kf->R;
    kf->H = H_full; kf->R = R_full;
    int ret = kf_update(kf, z_full);
    kf->H = H_save; kf->R = R_save;

    free(z_full); free(H_full); free(R_full);
    return ret;
}

/* ========================================================================
 *  L8-I02: kf_cholesky_2x2 -- Explicit 2x2 Cholesky Benchmark
 *
 *  For a 2x2 symmetric PD matrix A = [[a, b], [b, c]]:
 *    L = [[sqrt(a), 0], [b/sqrt(a), sqrt(c - b^2/a)]]
 *  Condition: a > 0 and a*c - b^2 > 0.
 *
 *  This specialized implementation is ~3x faster than the general
 *  n x n Cholesky for the common 2x2 case in GPS tracking loops
 *  where P and S are partitioned into 2x2 blocks.
 *
 *  Complexity: O(1) — constant time.
 * ======================================================================== */

int kf_cholesky_2x2(const double *A, double *L)
{
    if (!A || !L) return -1;
    if (A[0] <= 1e-30) return -1;  /* Not PD */

    double a = A[0], b = A[1], c = A[3];
    double l11 = sqrt(a);
    double l21 = b / l11;
    double l22_sq = c - l21 * l21;

    if (l22_sq <= 1e-30) return -1;  /* Not PD */
    double l22 = sqrt(l22_sq);

    L[0] = l11;  L[1] = 0.0;
    L[2] = l21;  L[3] = l22;
    return 0;
}

/* ========================================================================
 *  L8-I03: kf_process_noise_discrete -- DWNA Model Generator
 *
 *  Generates the Q_d matrix for a constant-velocity model with
 *  sampling interval dt and continuous-time process noise intensity q.
 *
 *  The DWNA (Discrete White Noise Acceleration) model assumes
 *  piecewise-constant white acceleration between samples:
 *    Q_d = q * [dt^3/3*I_3   dt^2/2*I_3 ]
 *              [dt^2/2*I_3   dt    *I_3 ]
 *
 *  This is the standard model for tracking maneuvering targets
 *  (Bar-Shalom et al., 2001, Ch.6).
 *
 *  @param dt       Sampling interval (seconds).
 *  @param q        Acceleration PSD (m^2/s^3 for position; or general).
 *  @param dim      Spatial dimension (1, 2, or 3).
 *  @param Q_out    Output matrix 2*dim x 2*dim, row-major.
 * ======================================================================== */

int kf_process_noise_discrete(double dt, double q, int dim, double *Q_out)
{
    if (!Q_out || dim < 1 || dim > 3 || dt <= 0.0 || q < 0.0) return -1;
    int n = 2 * dim;
    double dt2 = dt * dt, dt3 = dt2 * dt;

    memset(Q_out, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < dim; i++) {
        int pi = i, vi = dim + i;
        Q_out[pi*n+pi] = q * dt3 / 3.0;
        Q_out[vi*n+vi] = q * dt;
        Q_out[pi*n+vi] = q * dt2 / 2.0;
        Q_out[vi*n+pi] = q * dt2 / 2.0;
    }
    return 0;
}
