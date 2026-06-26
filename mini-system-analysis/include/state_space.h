/**
 * @file state_space.h
 * @brief L3 Math Structures - State-Space Representation
 *
 * State-space: dx/dt = A x + B u,  y = C x + D u  (CT)
 *              x[n+1] = A x[n] + B u[n], y[n] = C x[n] + D u[n]  (DT)
 *
 * Knowledge Coverage:
 *   L3-S01: CT state-space struct (A,B,C,D matrices)
 *   L3-S02: DT state-space struct
 *   L3-S03: Controllability matrix construction
 *   L3-S04: Observability matrix construction
 *   L3-S05: Controllability rank test (Kalman criterion)
 *   L3-S06: Observability rank test
 *   L3-S07: State transition matrix (matrix exponential e^{At})
 *   L3-S08: CT simulation (Euler, RK4)
 *   L3-S09: DT simulation (iterative update)
 *   L3-S10: Eigenvalue computation for stability
 *   L3-S11: Similarity transform (change of basis)
 *   L3-S12: Controllable canonical form
 *   L3-S13: Observable canonical form
 *   L3-S14: Lyapunov equation A P + P A^T + Q = 0
 *
 * Course Mapping:
 *   MIT 6.003 Ch.12; Stanford EE102A Ch.10;
 *   Berkeley EE123 Ch.8; Tsinghua Signals Ch.8
 */

#ifndef STATE_SPACE_H
#define STATE_SPACE_H

#include "system_defs.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- State-Space Data Structures ---- */

typedef struct {
    double *A;        /* n x n system matrix, row-major */
    double *B;        /* n x m input matrix */
    double *C;        /* p x n output matrix */
    double *D;        /* p x m feedthrough matrix */
    int     n;        /* Number of states */
    int     m;        /* Number of inputs */
    int     p;        /* Number of outputs */
    int     owns_data; /* 1 if matrices were malloc'd */
} ct_state_space_t;

typedef struct {
    double *A;        /* n x n system matrix, row-major */
    double *B;        /* n x m input matrix */
    double *C;        /* p x n output matrix */
    double *D;        /* p x m feedthrough matrix */
    int     n, m, p;
    int     owns_data;
} dt_state_space_t;

/* ---- Allocation ---- */

ct_state_space_t ct_ss_alloc(int n, int m, int p);
void             ct_ss_free(ct_state_space_t *ss);

dt_state_space_t dt_ss_alloc(int n, int m, int p);
void             dt_ss_free(dt_state_space_t *ss);

/* ---- Matrix Access Helpers ---- */

/* A(i,j) = A[i*n + j] for row-major storage */
#define CT_SS_A(ss, i, j) ((ss).A[(i)*(ss).n + (j)])
#define CT_SS_B(ss, i, j) ((ss).B[(i)*(ss).m + (j)])
#define CT_SS_C(ss, i, j) ((ss).C[(i)*(ss).n + (j)])
#define CT_SS_D(ss, i, j) ((ss).D[(i)*(ss).m + (j)])

#define DT_SS_A(ss, i, j) ((ss).A[(i)*(ss).n + (j)])
#define DT_SS_B(ss, i, j) ((ss).B[(i)*(ss).m + (j)])
#define DT_SS_C(ss, i, j) ((ss).C[(i)*(ss).n + (j)])
#define DT_SS_D(ss, i, j) ((ss).D[(i)*(ss).m + (j)])

/* ---- Controllability and Observability ---- */

int ct_ss_controllability_matrix(const ct_state_space_t *ss, double *Cc);
int ct_ss_controllability_rank(const ct_state_space_t *ss);
int ct_ss_is_controllable(const ct_state_space_t *ss);

int ct_ss_observability_matrix(const ct_state_space_t *ss, double *Oc);
int ct_ss_observability_rank(const ct_state_space_t *ss);
int ct_ss_is_observable(const ct_state_space_t *ss);

int dt_ss_controllability_matrix(const dt_state_space_t *ss, double *Cc);
int dt_ss_controllability_rank(const dt_state_space_t *ss);
int dt_ss_is_controllable(const dt_state_space_t *ss);

int dt_ss_observability_matrix(const dt_state_space_t *ss, double *Oc);
int dt_ss_observability_rank(const dt_state_space_t *ss);
int dt_ss_is_observable(const dt_state_space_t *ss);

/* ---- Matrix Exponential (State Transition Matrix) ---- */

int ct_ss_matrix_exponential(const ct_state_space_t *ss, double t,
                              double *expAt);

/* ---- Simulation ---- */

int ct_ss_simulate_euler(const ct_state_space_t *ss,
                          const double *x0, const double *u, int u_len,
                          double dt, int num_steps,
                          double *x_out, double *y_out);

int ct_ss_simulate_rk4(const ct_state_space_t *ss,
                        const double *x0, const double *u, int u_len,
                        double dt, int num_steps,
                        double *x_out, double *y_out);

int dt_ss_simulate(const dt_state_space_t *ss,
                    const double *x0, const double *u, int u_len,
                    double *x_out, double *y_out);

/* ---- Eigenvalue Computation ---- */

int ct_ss_eigenvalues(const ct_state_space_t *ss,
                       double complex *eigenvalues);

int dt_ss_eigenvalues(const dt_state_space_t *ss,
                       double complex *eigenvalues);

/* ---- Similarity Transform ---- */

int ct_ss_similarity_transform(const ct_state_space_t *ss,
                                const double *T,
                                ct_state_space_t *ss_transformed);

/* ---- Canonical Forms ---- */

int ct_ss_to_controllable_canonical(const ct_state_space_t *ss,
                                     ct_state_space_t *ss_ccf);

int ct_ss_to_observable_canonical(const ct_state_space_t *ss,
                                   ct_state_space_t *ss_ocf);

/* ---- Lyapunov Equation ---- */

int ct_ss_lyapunov(const ct_state_space_t *ss, const double *Q,
                    double *P);

/* ---- Transfer Function from State-Space ---- */

int ct_ss_to_transfer_function(const ct_state_space_t *ss,
                                ct_transfer_function_t *tf);

int dt_ss_to_transfer_function(const dt_state_space_t *ss,
                                dt_transfer_function_t *tf);

#ifdef __cplusplus
}
#endif

#endif /* STATE_SPACE_H */
