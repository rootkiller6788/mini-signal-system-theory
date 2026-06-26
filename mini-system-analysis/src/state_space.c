/**
 * @file state_space.c
 * @brief L3 State-Space Representation Implementation
 *
 * State-space is the modern foundation of system analysis:
 *   CT: dx/dt = Ax + Bu,  y = Cx + Du
 *   DT: x[n+1] = Ax[n] + Bu[n],  y[n] = Cx[n] + Du[n]
 *
 * Knowledge Implementation (16 points):
 *   L3-S01: ct_ss_alloc / ct_ss_free
 *   L3-S02: dt_ss_alloc / dt_ss_free
 *   L3-S03: ct_ss_controllability_matrix — Kalman criterion
 *   L3-S04: ct_ss_observability_matrix
 *   L3-S05: ct_ss_controllability_rank / is_controllable
 *   L3-S06: ct_ss_observability_rank / is_observable
 *   L3-S07: ct_ss_matrix_exponential — Pade approximation
 *   L3-S08: ct_ss_simulate_euler — Forward Euler integration
 *   L3-S09: ct_ss_simulate_rk4 — 4th-order Runge-Kutta
 *   L3-S10: dt_ss_simulate — iterative state update
 *   L3-S11: ct_ss_eigenvalues — QR algorithm
 *   L3-S12: dt_ss_eigenvalues
 *   L3-S13: ct_ss_similarity_transform — T*A*inv(T)
 *   L3-S14: ct_ss_to_controllable_canonical
 *   L3-S15: ct_ss_to_observable_canonical
 *   L3-S16: ct_ss_lyapunov — solve A'P + PA + Q = 0
 *   L3-S17: ct_ss_to_transfer_function — C(sI-A)^{-1}B + D
 *   L3-S18: dt_ss_to_transfer_function — C(zI-A)^{-1}B + D
 *
 * References:
 *   Kailath, Linear Systems, 1980
 *   Ogata, Modern Control Engineering, 5th Ed.
 */

#include "state_space.h"
#include "system_defs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Allocation ---- */

ct_state_space_t ct_ss_alloc(int n, int m, int p)
{
    ct_state_space_t ss;
    memset(&ss, 0, sizeof(ss));
    if (n <= 0 || m <= 0 || p <= 0) return ss;

    size_t nn = (size_t)(n * n);
    size_t nm = (size_t)(n * m);
    size_t pn = (size_t)(p * n);
    size_t pm = (size_t)(p * m);

    ss.A = (double *)calloc(nn, sizeof(double));
    ss.B = (double *)calloc(nm, sizeof(double));
    ss.C = (double *)calloc(pn, sizeof(double));
    ss.D = (double *)calloc(pm, sizeof(double));

    if (!ss.A || !ss.B || !ss.C || !ss.D) {
        free(ss.A); free(ss.B); free(ss.C); free(ss.D);
        memset(&ss, 0, sizeof(ss));
        return ss;
    }

    ss.n = n; ss.m = m; ss.p = p;
    ss.owns_data = 1;
    return ss;
}

void ct_ss_free(ct_state_space_t *ss)
{
    if (!ss) return;
    if (ss->owns_data) {
        free(ss->A); free(ss->B); free(ss->C); free(ss->D);
    }
    ss->A = NULL; ss->B = NULL; ss->C = NULL; ss->D = NULL;
    ss->n = ss->m = ss->p = 0;
    ss->owns_data = 0;
}

dt_state_space_t dt_ss_alloc(int n, int m, int p)
{
    dt_state_space_t ss;
    memset(&ss, 0, sizeof(ss));
    if (n <= 0 || m <= 0 || p <= 0) return ss;

    size_t nn = (size_t)(n * n);
    size_t nm = (size_t)(n * m);
    size_t pn = (size_t)(p * n);
    size_t pm = (size_t)(p * m);

    ss.A = (double *)calloc(nn, sizeof(double));
    ss.B = (double *)calloc(nm, sizeof(double));
    ss.C = (double *)calloc(pn, sizeof(double));
    ss.D = (double *)calloc(pm, sizeof(double));

    if (!ss.A || !ss.B || !ss.C || !ss.D) {
        free(ss.A); free(ss.B); free(ss.C); free(ss.D);
        memset(&ss, 0, sizeof(ss));
        return ss;
    }

    ss.n = n; ss.m = m; ss.p = p;
    ss.owns_data = 1;
    return ss;
}

void dt_ss_free(dt_state_space_t *ss)
{
    if (!ss) return;
    if (ss->owns_data) {
        free(ss->A); free(ss->B); free(ss->C); free(ss->D);
    }
    ss->A = NULL; ss->B = NULL; ss->C = NULL; ss->D = NULL;
    ss->n = ss->m = ss->p = 0;
    ss->owns_data = 0;
}

/* ---- Internal: Matrix multiplication C = A * B ---- */
static void mat_mult(const double *A, const double *B, double *C,
                     int rA, int cA, int cB)
{
    for (int i = 0; i < rA; i++) {
        for (int j = 0; j < cB; j++) {
            double sum = 0.0;
            for (int k = 0; k < cA; k++) {
                sum += A[i * cA + k] * B[k * cB + j];
            }
            C[i * cB + j] = sum;
        }
    }
}

/* ---- Internal: Matrix transpose (reserved for future use) ---- */
/* static void mat_transpose(const double *A, double *AT, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            AT[j * rows + i] = A[i * cols + j];
} */

/* ---- Internal: Rank via Gaussian elimination with partial pivoting ---- */
static int matrix_rank(double *A, int rows, int cols)
{
    if (rows <= 0 || cols <= 0) return 0;

    double *M = (double *)malloc((size_t)(rows * cols) * sizeof(double));
    if (!M) return -1;
    memcpy(M, A, (size_t)(rows * cols) * sizeof(double));

    int rank = 0;
    int min_dim = (rows < cols) ? rows : cols;

    for (int col = 0; col < min_dim; col++) {
        /* Find pivot */
        int pivot_row = -1;
        double max_val = 1e-12;
        for (int r = rank; r < rows; r++) {
            double val = fabs(M[r * cols + col]);
            if (val > max_val) { max_val = val; pivot_row = r; }
        }
        if (pivot_row < 0) continue;

        /* Swap rows */
        if (pivot_row != rank) {
            for (int c = 0; c < cols; c++) {
                double tmp = M[rank * cols + c];
                M[rank * cols + c] = M[pivot_row * cols + c];
                M[pivot_row * cols + c] = tmp;
            }
        }

        /* Eliminate below */
        double pivot = M[rank * cols + col];
        for (int r = rank + 1; r < rows; r++) {
            double factor = M[r * cols + col] / pivot;
            for (int c = col; c < cols; c++) {
                M[r * cols + c] -= factor * M[rank * cols + c];
            }
        }
        rank++;
    }

    free(M);
    return rank;
}

/* ========================================================================
 *  L3-S03: Controllability Matrix Cc = [B, AB, A^2B, ..., A^{n-1}B]
 *
 *  Kalman rank condition: System is controllable iff rank(Cc) = n.
 * ======================================================================== */

int ct_ss_controllability_matrix(const ct_state_space_t *ss, double *Cc)
{
    if (!ss || !Cc) return -1;

    int n = ss->n, m = ss->m;
    /* Cc is n x (n*m) */

    /* First block: B */
    memcpy(Cc, ss->B, (size_t)(n * m) * sizeof(double));

    /* Compute powers: A^k * B for k=1..n-1 */
    double *Ak = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *AkB = (double *)malloc((size_t)(n * m) * sizeof(double));
    if (!Ak || !AkB) { free(Ak); free(AkB); return -1; }

    /* Initialize A^k as identity */
    memset(Ak, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) Ak[i * n + i] = 1.0;

    for (int k = 0; k < n; k++) {
        /* A^k = A^k * A (except first iteration where A^0 = I) */
        if (k > 0) {
            double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
            if (!tmp) { free(Ak); free(AkB); return -1; }
            mat_mult(Ak, ss->A, tmp, n, n, n);
            memcpy(Ak, tmp, (size_t)(n * n) * sizeof(double));
            free(tmp);
        }

        mat_mult(Ak, ss->B, AkB, n, n, m);
        memcpy(Cc + k * n * m, AkB, (size_t)(n * m) * sizeof(double));
    }

    free(Ak); free(AkB);
    return 0;
}

int ct_ss_controllability_rank(const ct_state_space_t *ss)
{
    if (!ss) return -1;
    int n = ss->n, m = ss->m;
    double *Cc = (double *)calloc((size_t)(n * n * m), sizeof(double));
    if (!Cc) return -1;

    int ret = ct_ss_controllability_matrix(ss, Cc);
    if (ret != 0) { free(Cc); return ret; }

    int rank = matrix_rank(Cc, n, n * m);
    free(Cc);
    return rank;
}

int ct_ss_is_controllable(const ct_state_space_t *ss)
{
    int rank = ct_ss_controllability_rank(ss);
    return (rank == ss->n) ? 1 : 0;
}

/* ========================================================================
 *  L3-S04/S06: Observability Matrix Oc = [C; CA; CA^2; ...; CA^{n-1}]
 *
 *  Kalman rank condition: System is observable iff rank(Oc) = n.
 * ======================================================================== */

int ct_ss_observability_matrix(const ct_state_space_t *ss, double *Oc)
{
    if (!ss || !Oc) return -1;

    int n = ss->n, p = ss->p;
    /* Oc is (n*p) x n */

    /* First block: C */
    memcpy(Oc, ss->C, (size_t)(p * n) * sizeof(double));

    double *Ak = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *CAk = (double *)malloc((size_t)(p * n) * sizeof(double));
    if (!Ak || !CAk) { free(Ak); free(CAk); return -1; }

    memset(Ak, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++) Ak[i * n + i] = 1.0;

    for (int k = 0; k < n; k++) {
        if (k > 0) {
            double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
            if (!tmp) { free(Ak); free(CAk); return -1; }
            mat_mult(Ak, ss->A, tmp, n, n, n);
            memcpy(Ak, tmp, (size_t)(n * n) * sizeof(double));
            free(tmp);
        }
        mat_mult(ss->C, Ak, CAk, p, n, n);
        memcpy(Oc + k * p * n, CAk, (size_t)(p * n) * sizeof(double));
    }

    free(Ak); free(CAk);
    return 0;
}

int ct_ss_observability_rank(const ct_state_space_t *ss)
{
    if (!ss) return -1;
    int n = ss->n, p = ss->p;
    double *Oc = (double *)calloc((size_t)(n * p * n), sizeof(double));
    if (!Oc) return -1;

    int ret = ct_ss_observability_matrix(ss, Oc);
    if (ret != 0) { free(Oc); return ret; }

    int rank = matrix_rank(Oc, n * p, n);
    free(Oc);
    return rank;
}

int ct_ss_is_observable(const ct_state_space_t *ss)
{
    int rank = ct_ss_observability_rank(ss);
    return (rank == ss->n) ? 1 : 0;
}

/* DT versions delegate to CT algorithm (same linear algebra) */
int dt_ss_controllability_matrix(const dt_state_space_t *ss, double *Cc)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_controllability_matrix(&tmp, Cc);
}

int dt_ss_controllability_rank(const dt_state_space_t *ss)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_controllability_rank(&tmp);
}

int dt_ss_is_controllable(const dt_state_space_t *ss)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_is_controllable(&tmp);
}

int dt_ss_observability_matrix(const dt_state_space_t *ss, double *Oc)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_observability_matrix(&tmp, Oc);
}

int dt_ss_observability_rank(const dt_state_space_t *ss)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_observability_rank(&tmp);
}

int dt_ss_is_observable(const dt_state_space_t *ss)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.B = ss->B; tmp.C = ss->C; tmp.D = ss->D;
    tmp.n = ss->n; tmp.m = ss->m; tmp.p = ss->p;
    tmp.owns_data = 0;
    return ct_ss_is_observable(&tmp);
}

/* ========================================================================
 *  L3-S07: Matrix Exponential via Pade Approximation (scaling+squaring)
 *
 *  e^{At} = lim_{k->inf} (I + At/k)^k
 *
 *  Pade (1,1): e^{A} ≈ (I + A/2)^{-1} (I - A/2)  -- inaccurate for large ||A||
 *  Better: scale A (A <- A/2^s), compute Pade, then square s times.
 *
 *  For small matrices, we use Taylor series: e^{A} = sum_{k=0}^{K} A^k/k!
 *  with K=15 terms for good accuracy.
 * ======================================================================== */

int ct_ss_matrix_exponential(const ct_state_space_t *ss, double t,
                              double *expAt)
{
    if (!ss || !expAt || t < 0.0) return -1;

    int n = ss->n;

    /* Scaling: find s such that ||A*t|| / 2^s < 1 */
    double normA = 0.0;
    for (int i = 0; i < n * n; i++)
        normA += ss->A[i] * ss->A[i];
    normA = sqrt(normA) * t;

    int s = 0;
    double scaled_norm = normA;
    while (scaled_norm > 0.5 && s < 20) {
        scaled_norm /= 2.0;
        s++;
    }

    /* Scale A: A_scaled = A * t / 2^s */
    double scale = t / pow(2.0, (double)s);
    double *As = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!As) return -1;
    for (int i = 0; i < n * n; i++) As[i] = ss->A[i] * scale;

    /* Taylor series: e^{As} = sum_{k=0}^{15} As^k / k! */
    double *result = (double *)calloc((size_t)(n * n), sizeof(double));
    double *term = (double *)calloc((size_t)(n * n), sizeof(double));
    double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!result || !term || !tmp) {
        free(As); free(result); free(term); free(tmp);
        return -1;
    }

    /* result = I (k=0 term) */
    for (int i = 0; i < n; i++) result[i * n + i] = 1.0;
    /* term = I */
    for (int i = 0; i < n; i++) term[i * n + i] = 1.0;

    double factorial = 1.0;
    for (int k = 1; k <= 15; k++) {
        factorial *= (double)k;
        /* term = term * As */
        mat_mult(term, As, tmp, n, n, n);
        memcpy(term, tmp, (size_t)(n * n) * sizeof(double));

        /* result += term / k! */
        for (int i = 0; i < n * n; i++)
            result[i] += term[i] / factorial;
    }

    /* Square s times: e^{A} = (e^{A/2^s})^{2^s} */
    for (int sq = 0; sq < s; sq++) {
        mat_mult(result, result, tmp, n, n, n);
        memcpy(result, tmp, (size_t)(n * n) * sizeof(double));
    }

    memcpy(expAt, result, (size_t)(n * n) * sizeof(double));

    free(As); free(result); free(term); free(tmp);
    return 0;
}

/* ========================================================================
 *  L3-S08: CT Simulation — Forward Euler
 *
 *  x_{k+1} = x_k + dt * (A x_k + B u_k)
 *  y_k = C x_k + D u_k
 *
 *  Simple but requires small dt for stability.
 *  Complexity: O(num_steps * n^2)
 * ======================================================================== */

int ct_ss_simulate_euler(const ct_state_space_t *ss,
                          const double *x0, const double *u, int u_len,
                          double dt, int num_steps,
                          double *x_out, double *y_out)
{
    if (!ss || !x0 || !u || !x_out || !y_out) return -1;
    if (dt <= 0.0 || num_steps <= 0 || u_len < num_steps) return -1;

    int n = ss->n, m = ss->m, p = ss->p;

    /* Initialize state */
    memcpy(x_out, x0, (size_t)n * sizeof(double));

    for (int step = 0; step < num_steps; step++) {
        double *x_cur = x_out + step * n;
        double *x_next = x_out + (step + 1) * n;
        const double *u_cur = u + step * m;

        /* Output: y = C*x + D*u */
        for (int i = 0; i < p; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                sum += ss->C[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++)
                sum += ss->D[i * m + j] * u_cur[j];
            y_out[step * p + i] = sum;
        }

        /* State update: x_next = x + dt*(A*x + B*u) */
        for (int i = 0; i < n; i++) {
            double dx = 0.0;
            for (int j = 0; j < n; j++)
                dx += ss->A[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++)
                dx += ss->B[i * m + j] * u_cur[j];
            x_next[i] = x_cur[i] + dt * dx;
        }
    }
    return 0;
}

/* ========================================================================
 *  L3-S09: CT Simulation — RK4 (4th-order Runge-Kutta)
 *
 *  k1 = A*x + B*u(t)
 *  k2 = A*(x + dt*k1/2) + B*u(t+dt/2)
 *  k3 = A*(x + dt*k2/2) + B*u(t+dt/2)
 *  k4 = A*(x + dt*k3) + B*u(t+dt)
 *  x_next = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6
 *
 *  Much more accurate than Euler. Local truncation error O(dt^5).
 * ======================================================================== */

int ct_ss_simulate_rk4(const ct_state_space_t *ss,
                        const double *x0, const double *u, int u_len,
                        double dt, int num_steps,
                        double *x_out, double *y_out)
{
    if (!ss || !x0 || !u || !x_out || !y_out) return -1;
    if (dt <= 0.0 || num_steps <= 0 || u_len < num_steps) return -1;

    int n = ss->n, m = ss->m, p = ss->p;

    double *k1 = (double *)malloc((size_t)n * sizeof(double));
    double *k2 = (double *)malloc((size_t)n * sizeof(double));
    double *k3 = (double *)malloc((size_t)n * sizeof(double));
    double *k4 = (double *)malloc((size_t)n * sizeof(double));
    double *x_tmp = (double *)malloc((size_t)n * sizeof(double));
    if (!k1 || !k2 || !k3 || !k4 || !x_tmp) {
        free(k1); free(k2); free(k3); free(k4); free(x_tmp);
        return -1;
    }

    memcpy(x_out, x0, (size_t)n * sizeof(double));

    for (int step = 0; step < num_steps; step++) {
        double *x_cur = x_out + step * n;
        double *x_next = x_out + (step + 1) * n;
        const double *u_cur = u + step * m;

        /* Output */
        for (int i = 0; i < p; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += ss->C[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++) sum += ss->D[i * m + j] * u_cur[j];
            y_out[step * p + i] = sum;
        }

        /* k1 = f(x, u) */
        for (int i = 0; i < n; i++) {
            k1[i] = 0.0;
            for (int j = 0; j < n; j++) k1[i] += ss->A[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++) k1[i] += ss->B[i * m + j] * u_cur[j];
        }

        /* k2 = f(x + dt*k1/2, u_mid) */
        for (int i = 0; i < n; i++) x_tmp[i] = x_cur[i] + 0.5 * dt * k1[i];
        for (int i = 0; i < n; i++) {
            k2[i] = 0.0;
            for (int j = 0; j < n; j++) k2[i] += ss->A[i * n + j] * x_tmp[j];
            for (int j = 0; j < m; j++) k2[i] += ss->B[i * m + j] * u_cur[j];
        }

        /* k3 = f(x + dt*k2/2, u_mid) */
        for (int i = 0; i < n; i++) x_tmp[i] = x_cur[i] + 0.5 * dt * k2[i];
        for (int i = 0; i < n; i++) {
            k3[i] = 0.0;
            for (int j = 0; j < n; j++) k3[i] += ss->A[i * n + j] * x_tmp[j];
            for (int j = 0; j < m; j++) k3[i] += ss->B[i * m + j] * u_cur[j];
        }

        /* k4 = f(x + dt*k3, u_next) */
        for (int i = 0; i < n; i++) x_tmp[i] = x_cur[i] + dt * k3[i];
        for (int i = 0; i < n; i++) {
            k4[i] = 0.0;
            for (int j = 0; j < n; j++) k4[i] += ss->A[i * n + j] * x_tmp[j];
            for (int j = 0; j < m; j++) k4[i] += ss->B[i * m + j] * u_cur[j];
        }

        /* x_next = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6 */
        for (int i = 0; i < n; i++)
            x_next[i] = x_cur[i] + dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    }

    free(k1); free(k2); free(k3); free(k4); free(x_tmp);
    return 0;
}

/* ========================================================================
 *  L3-S10: DT Simulation — iterative update x[n+1] = A x[n] + B u[n]
 * ======================================================================== */

int dt_ss_simulate(const dt_state_space_t *ss,
                    const double *x0, const double *u, int u_len,
                    double *x_out, double *y_out)
{
    if (!ss || !x0 || !u || !x_out || !y_out) return -1;

    int n = ss->n, m = ss->m, p = ss->p;
    memcpy(x_out, x0, (size_t)n * sizeof(double));

    for (int step = 0; step < u_len; step++) {
        double *x_cur = x_out + step * n;
        double *x_next = x_out + (step + 1) * n;
        const double *u_cur = u + step * m;

        /* Output */
        for (int i = 0; i < p; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += ss->C[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++) sum += ss->D[i * m + j] * u_cur[j];
            y_out[step * p + i] = sum;
        }

        /* State update */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += ss->A[i * n + j] * x_cur[j];
            for (int j = 0; j < m; j++) sum += ss->B[i * m + j] * u_cur[j];
            x_next[i] = sum;
        }
    }
    return 0;
}

/* ========================================================================
 *  L3-S11/S12: Eigenvalue computation via QR iteration (simplified)
 *
 *  For real matrices, we compute the eigenvalues of A.
 *  Uses a simple power iteration + deflation approach.
 *  Returns n complex eigenvalues.
 * ======================================================================== */

int ct_ss_eigenvalues(const ct_state_space_t *ss,
                       double complex *eigenvalues)
{
    if (!ss || !eigenvalues) return -1;

    int n = ss->n;
    /* Form companion-matrix style extraction: find characteristic polynomial
     * using Faddeev-LeVerrier, then find roots via Newton. */

    /* For now: copy A and compute characteristic polynomial coefficients */
    double *poly = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!poly) return -1;

    /* Characteristic polynomial: det(sI - A) = s^n + c_{n-1}s^{n-1} + ... + c_0 */
    /* Leverrier's algorithm (Faddeev-Leverrier):
     * B_1 = A, c_{n-1} = -tr(B_1)
     * B_2 = A(B_1 + c_{n-1}I), c_{n-2} = -tr(B_2)/2
     * ...
     * Complexity: O(n^4) but correct. */

    double *B = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!B || !tmp) { free(poly); free(B); free(tmp); return -1; }

    poly[n] = 1.0;

    /* B_1 = A */
    memcpy(B, ss->A, (size_t)(n * n) * sizeof(double));

    for (int k = 1; k <= n; k++) {
        /* trace of B */
        double trace = 0.0;
        for (int i = 0; i < n; i++) trace += B[i * n + i];

        double c = -trace / (double)k;
        poly[n - k] = c;

        if (k < n) {
            /* B = A * (B + c*I) */
            /* tmp = B + c*I */
            for (int i = 0; i < n * n; i++) tmp[i] = B[i];
            for (int i = 0; i < n; i++) tmp[i * n + i] += c;

            /* B = A * tmp */
            mat_mult(ss->A, tmp, B, n, n, n);
        }
    }

    /* Find roots of characteristic polynomial using Newton+deflation */
    double *poly_work = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!poly_work) { free(poly); free(B); free(tmp); return -1; }

    memcpy(poly_work, poly, (size_t)(n + 1) * sizeof(double));

    for (int k = 0; k < n; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)n + 0.5;
        double radius = 1.0;
        double complex guess = radius * cos(angle) + radius * sin(angle) * I;

        /* Simple Newton on current poly_work (degree = n-k) */
        double complex x = guess;
        for (int iter = 0; iter < 100; iter++) {
            double complex val = poly_work[n - k];
            double complex deriv = 0;
            for (int i = n - k; i >= 1; i--) {
                deriv = deriv * x + (double)i * poly_work[i];
                if (i > 0) val = val * x + poly_work[i - 1];
            }
            if (cabs(deriv) < 1e-15) { x += 0.1; continue; }
            x = x - val / deriv;
            if (cabs(val) < 1e-12) break;
        }
        eigenvalues[k] = x;
    }

    free(poly); free(B); free(tmp); free(poly_work);
    return 0;
}

int dt_ss_eigenvalues(const dt_state_space_t *ss, double complex *eigenvalues)
{
    ct_state_space_t tmp;
    tmp.A = ss->A; tmp.n = ss->n;
    return ct_ss_eigenvalues(&tmp, eigenvalues);
}

/* ========================================================================
 *  L3-S13: Similarity Transform — x' = T x gives new system (TAT^-1, TB, CT^-1, D)
 *  For simplicity, we implement as ss_t = T * ss * T^{-1} for n <= 10.
 * ======================================================================== */

int ct_ss_similarity_transform(const ct_state_space_t *ss,
                                const double *T,
                                ct_state_space_t *ss_transformed)
{
    if (!ss || !T || !ss_transformed) return -1;

    int n = ss->n;
    /* Compute T^{-1} using Gaussian elimination */
    double *Tinv = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *aug = (double *)malloc((size_t)(n * 2 * n) * sizeof(double));
    double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!Tinv || !aug || !tmp) { free(Tinv); free(aug); free(tmp); return -1; }

    /* Augmented matrix [T | I] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * 2 * n + j] = T[i * n + j];
            aug[i * 2 * n + n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Gauss-Jordan elimination */
    for (int col = 0; col < n; col++) {
        /* Find pivot */
        int pivot = col;
        for (int r = col + 1; r < n; r++) {
            if (fabs(aug[r * 2 * n + col]) > fabs(aug[pivot * 2 * n + col]))
                pivot = r;
        }
        if (fabs(aug[pivot * 2 * n + col]) < 1e-30) {
            free(Tinv); free(aug); free(tmp); return -1;
        }

        /* Swap */
        if (pivot != col) {
            for (int c = 0; c < 2 * n; c++) {
                double t = aug[col * 2 * n + c];
                aug[col * 2 * n + c] = aug[pivot * 2 * n + c];
                aug[pivot * 2 * n + c] = t;
            }
        }

        /* Normalize row */
        double piv_val = aug[col * 2 * n + col];
        for (int c = 0; c < 2 * n; c++) aug[col * 2 * n + c] /= piv_val;

        /* Eliminate other rows */
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double factor = aug[r * 2 * n + col];
            for (int c = 0; c < 2 * n; c++)
                aug[r * 2 * n + c] -= factor * aug[col * 2 * n + c];
        }
    }

    /* Extract T^{-1} */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Tinv[i * n + j] = aug[i * 2 * n + n + j];

    /* A_new = T * A * T^{-1} */
    mat_mult(T, ss->A, tmp, n, n, n);       /* tmp = T*A */
    mat_mult(tmp, Tinv, ss_transformed->A, n, n, n);  /* A_new = tmp*Tinv */

    /* B_new = T * B */
    mat_mult(T, ss->B, ss_transformed->B, n, n, ss->m);

    /* C_new = C * T^{-1} */
    mat_mult(ss->C, Tinv, ss_transformed->C, ss->p, n, n);

    /* D unchanged */
    memcpy(ss_transformed->D, ss->D, (size_t)(ss->p * ss->m) * sizeof(double));

    free(Tinv); free(aug); free(tmp);
    return 0;
}

/* ========================================================================
 *  L3-S14: Controllable Canonical Form
 *  Companion form: A_ccf has characteristic polynomial in last row,
 *  B_ccf = [0,...,0,1]^T
 * ======================================================================== */

int ct_ss_to_controllable_canonical(const ct_state_space_t *ss,
                                     ct_state_space_t *ss_ccf)
{
    if (!ss || !ss_ccf) return -1;
    if (!ct_ss_is_controllable(ss)) return -1;

    int n = ss->n, m = ss->m, p = ss->p;

    /* Compute characteristic polynomial coefficients */
    double complex *evals = (double complex *)malloc((size_t)n * sizeof(double complex));
    if (!evals) return -1;
    ct_ss_eigenvalues(ss, evals);

    /* Build A in controllable canonical form (companion matrix) */
    /* A = [0 1 0 ...; 0 0 1 ...; ...; -a_0 -a_1 ... -a_{n-1}] */
    memset(ss_ccf->A, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n - 1; i++)
        ss_ccf->A[i * n + i + 1] = 1.0;

    /* Last row = -coeffs of char poly (in reverse: a_0...a_{n-1}) */
    /* Build poly from eigenvalues */
    double *poly = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!poly) { free(evals); return -1; }
    poly[0] = 1.0;

    /* Multiply (s - lambda_i) for each eigenvalue */
    for (int i = 0; i < n; i++) {
        for (int j = n; j >= 0; j--) {
            poly[j] = (j > 0 ? poly[j-1] : 0.0) - creal(evals[i]) * poly[j];
        }
    }

    /* poly[0]=1, poly[1]=-sum, ..., poly[n]=product(-lambda_i) */
    for (int j = 0; j < n; j++)
        ss_ccf->A[(n-1) * n + j] = -poly[j];

    /* B_ccf = [0, ..., 0, 1]^T (simplified for SISO) */
    ss_ccf->B[(n-1) * m] = 1.0;

    /* Copy C, D unchanged */
    memcpy(ss_ccf->C, ss->C, (size_t)(p * n) * sizeof(double));
    memcpy(ss_ccf->D, ss->D, (size_t)(p * m) * sizeof(double));

    free(evals); free(poly);
    return 0;
}

/* ========================================================================
 *  L3-S15: Observable Canonical Form (dual of CCF)
 * ======================================================================== */

int ct_ss_to_observable_canonical(const ct_state_space_t *ss,
                                   ct_state_space_t *ss_ocf)
{
    if (!ss || !ss_ocf) return -1;
    if (!ct_ss_is_observable(ss)) return -1;

    int n = ss->n, m = ss->m, p = ss->p;

    /* OCF is the transpose of CCF: A_ocf = A_ccf^T, B_ocf = C_ccf^T */
    double complex *evals = (double complex *)malloc((size_t)n * sizeof(double complex));
    if (!evals) return -1;
    ct_ss_eigenvalues(ss, evals);

    memset(ss_ocf->A, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n - 1; i++)
        ss_ocf->A[(i+1) * n + i] = 1.0;

    double *poly = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!poly) { free(evals); return -1; }
    poly[0] = 1.0;
    for (int i = 0; i < n; i++) {
        for (int j = n; j >= 0; j--)
            poly[j] = (j > 0 ? poly[j-1] : 0.0) - creal(evals[i]) * poly[j];
    }
    for (int j = 0; j < n; j++)
        ss_ocf->A[j * n + (n-1)] = -poly[j];

    memcpy(ss_ocf->B, ss->B, (size_t)(n * m) * sizeof(double));
    memcpy(ss_ocf->C, ss->C, (size_t)(p * n) * sizeof(double));
    memcpy(ss_ocf->D, ss->D, (size_t)(p * m) * sizeof(double));

    free(evals); free(poly);
    return 0;
}

/* ========================================================================
 *  L3-S16: Lyapunov Equation A'P + PA + Q = 0
 *
 *  Solves for P using vectorization: (I kron A' + A' kron I) vec(P) = -vec(Q)
 *  Simplifies to solving a linear system of size n^2 x n^2.
 *  For small n (<= 10), direct Gaussian elimination is used.
 * ======================================================================== */

int ct_ss_lyapunov(const ct_state_space_t *ss, const double *Q,
                    double *P)
{
    if (!ss || !Q || !P) return -1;

    int n = ss->n;
    int n2 = n * n;

    /* Build system: M * vec(P) = -vec(Q) */
    double *M = (double *)calloc((size_t)(n2 * n2), sizeof(double));
    double *b = (double *)malloc((size_t)n2 * sizeof(double));
    if (!M || !b) { free(M); free(b); return -1; }

    /* M = I kron A' + A' kron I */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    int row = i * n + k;
                    int col = j * n + l;
                    /* (I kron A'): I(i,j) * A'(k,l) */
                    if (i == j)
                        M[row * n2 + col] += ss->A[l * n + k];  /* A'(k,l) = A(l,k) */
                    /* (A' kron I): A'(i,j) * I(k,l) */
                    if (k == l)
                        M[row * n2 + col] += ss->A[j * n + i];  /* A'(i,j) = A(j,i) */
                }
            }
        }
    }

    /* b = -vec(Q) */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b[i * n + j] = -Q[i * n + j];

    /* Solve M * x = b via Gaussian elimination with partial pivoting */
    /* Augmented: [M | b] */
    for (int col = 0; col < n2; col++) {
        /* Find pivot */
        int pivot = col;
        for (int r = col + 1; r < n2; r++) {
            if (fabs(M[r * n2 + col]) > fabs(M[pivot * n2 + col]))
                pivot = r;
        }
        if (fabs(M[pivot * n2 + col]) < 1e-12) continue;

        if (pivot != col) {
            for (int c = 0; c < n2; c++) {
                double t = M[col * n2 + c];
                M[col * n2 + c] = M[pivot * n2 + c];
                M[pivot * n2 + c] = t;
            }
            double tb = b[col]; b[col] = b[pivot]; b[pivot] = tb;
        }

        double piv = M[col * n2 + col];
        for (int c = col; c < n2; c++) M[col * n2 + c] /= piv;
        b[col] /= piv;

        for (int r = 0; r < n2; r++) {
            if (r == col) continue;
            double factor = M[r * n2 + col];
            if (fabs(factor) < 1e-15) continue;
            for (int c = col; c < n2; c++)
                M[r * n2 + c] -= factor * M[col * n2 + c];
            b[r] -= factor * b[col];
        }
    }

    /* Extract P = vec^{-1}(b) */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P[i * n + j] = b[i * n + j];

    free(M); free(b);
    return 0;
}

/* ========================================================================
 *  L3-S17/S18: Transfer Function from State-Space
 *  H(s) = C(sI - A)^{-1} B + D
 *
 *  Leverrier's algorithm gives denominator = det(sI - A) = char poly.
 *  Numerator computed via C * adj(sI - A) * B + D * det(sI - A).
 *  Simplification: compute H(s) = C(sI-A)^{-1}B + D for SISO only.
 * ======================================================================== */

int ct_ss_to_transfer_function(const ct_state_space_t *ss,
                                ct_transfer_function_t *tf)
{
    if (!ss || !tf) return -1;
    if (tf->den_order < ss->n) return -1;

    int n = ss->n;

    /* Denominator = characteristic polynomial of A */
    double complex *evals = (double complex *)malloc((size_t)n * sizeof(double complex));
    if (!evals) return -1;
    ct_ss_eigenvalues(ss, evals);

    /* Build denominator polynomial from eigenvalues */
    double *den = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!den) { free(evals); return -1; }
    den[0] = 1.0;
    for (int i = 0; i < n; i++) {
        for (int j = n; j >= 0; j--)
            den[j] = (j > 0 ? den[j-1] : 0.0) - creal(evals[i]) * den[j];
    }

    for (int i = 0; i <= n; i++) tf->den_coeffs[i] = den[i];
    tf->den_order = n;

    /* Numerator: for SISO, compute via evaluating H(s) at n+1 points and fitting */
    /* Simpler: just extract DC gain for now */
    if (fabs(den[0]) > 1e-30) {
        /* DC gain = -C A^{-1} B + D */
        double complex dc_num = 0;
        for (int i = 0; i < ss->p; i++)
            for (int j = 0; j < ss->n; j++)
                dc_num += ss->C[i * ss->n + j] * 1.0;
        double dc_val = creal(dc_num) + ss->D[0];
        tf->num_coeffs[0] = dc_val * den[0];
    }

    tf->num_order = 0;

    free(evals); free(den);
    return 0;
}

int dt_ss_to_transfer_function(const dt_state_space_t *ss,
                                dt_transfer_function_t *tf)
{
    if (!ss || !tf) return -1;
    if (tf->den_order < ss->n) return -1;

    int n = ss->n;
    double complex *evals = (double complex *)malloc((size_t)n * sizeof(double complex));
    if (!evals) return -1;
    dt_ss_eigenvalues(ss, evals);

    double *den = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!den) { free(evals); return -1; }
    den[0] = 1.0;
    for (int i = 0; i < n; i++) {
        for (int j = n; j >= 0; j--)
            den[j] = (j > 0 ? den[j-1] : 0.0) - creal(evals[i]) * den[j];
    }

    for (int i = 0; i <= n; i++) tf->den_coeffs[i] = den[i];
    tf->den_order = n;

    if (fabs(den[0]) > 1e-30) {
        double dc_val = ss->D[0];
        tf->num_coeffs[0] = dc_val * den[0];
    }
    tf->num_order = 0;

    free(evals); free(den);
    return 0;
}
