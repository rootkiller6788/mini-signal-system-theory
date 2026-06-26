/**
 * @file ex3_state_space.c
 * @brief L6 Example: State-Space Simulation and Analysis
 *
 * Demonstrates state-space representation, controllability/observability
 * tests, and time-domain simulation using Euler and RK4 methods.
 *
 * Canonical Problem:
 *   Mass-spring-damper: m*x'' + c*x' + k*x = F(t)
 *   State: x1 = position, x2 = velocity
 *
 * Course: MIT 6.003 Ch.12, Stanford EE102A Ch.10
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/system_defs.h"
#include "../include/state_space.h"
#include "../include/stability.h"

int main(void)
{
    printf("=== Example 3: State-Space Analysis ===\n\n");

    /* ── Mass-Spring-Damper System ── */
    printf("--- Mass-Spring-Damper (m=1, c=0.5, k=2) ---\n");
    printf("State: x1=position, x2=velocity\n");

    double m = 1.0, c = 0.5, k = 2.0;

    ct_state_space_t ss = ct_ss_alloc(2, 1, 1);
    /* A = [0, 1; -k/m, -c/m] = [0, 1; -2, -0.5] */
    ss.A[0*2+0] = 0.0;   ss.A[0*2+1] = 1.0;
    ss.A[1*2+0] = -k/m;  ss.A[1*2+1] = -c/m;
    /* B = [0; 1/m] = [0; 1] */
    ss.B[0] = 0.0;  ss.B[1] = 1.0/m;
    /* C = [1, 0] (measure position) */
    ss.C[0] = 1.0;  ss.C[1] = 0.0;
    /* D = [0] */
    ss.D[0] = 0.0;

    /* Controllability test */
    printf("  Controllable: %s\n",
           ct_ss_is_controllable(&ss) ? "YES" : "NO");

    /* Observability test */
    printf("  Observable: %s\n",
           ct_ss_is_observable(&ss) ? "YES" : "NO");

    /* Eigenvalue analysis */
    double complex evals[2];
    ct_ss_eigenvalues(&ss, evals);
    printf("  Eigenvalues: %.3f%+.3fj, %.3f%+.3fj\n",
           creal(evals[0]), cimag(evals[0]),
           creal(evals[1]), cimag(evals[1]));

    /* Natural frequency and damping */
    double wn = sqrt(k/m);  /* sqrt(2) = 1.414 */
    double zeta = c / (2.0 * sqrt(k*m));  /* 0.5/(2*sqrt(2)) = 0.177 */
    printf("  Natural frequency: %.3f rad/s\n", wn);
    printf("  Damping ratio: %.3f\n", zeta);

    /* ── Time-Domain Simulation ── */
    printf("\n--- RK4 Simulation (step response, F=1) ---\n");

    int n_steps = 200;
    double dt = 0.1;
    double x0[2] = {0.0, 0.0};
    double u[200];
    for (int i = 0; i < 200; i++) u[i] = 1.0;  /* Step input */

    double *x = (double *)calloc((size_t)(2 * (n_steps + 1)), sizeof(double));
    double *y = (double *)calloc((size_t)(n_steps + 1), sizeof(double));

    ct_ss_simulate_rk4(&ss, x0, u, n_steps, dt, n_steps, x, y);

    /* Print selected time points */
    printf("  t=1.0: pos=%.4f, vel=%.4f\n", x[10*2+0], x[10*2+1]);
    printf("  t=5.0: pos=%.4f, vel=%.4f\n", x[50*2+0], x[50*2+1]);
    printf("  t=10.0: pos=%.4f, vel=%.4f\n", x[100*2+0], x[100*2+1]);
    printf("  Steady-state position ≈ %.4f (expected = 1/k = 0.5)\n",
           x[n_steps*2+0]);

    printf("  Settling time estimate: %.3f sec\n",
           ct_settling_time_estimate(zeta, wn));
    printf("  Overshoot estimate: %.2f%%\n",
           ct_overshoot_estimate(zeta));

    free(x); free(y);
    ct_ss_free(&ss);

    /* ── Discrete-Time State-Space ── */
    printf("\n--- DT State-Space (accumulator) ---\n");
    dt_state_space_t dt_ss = dt_ss_alloc(1, 1, 1);
    dt_ss.A[0] = 1.0;   /* x[n+1] = x[n] + u[n] */
    dt_ss.B[0] = 1.0;
    dt_ss.C[0] = 1.0;
    dt_ss.D[0] = 0.0;

    double x0_dt[1] = {0.0};
    double u_dt[10] = {1,1,1,1,1, 1,1,1,1,1};
    double x_dt[11] = {0};
    double y_dt[10] = {0};

    dt_ss_simulate(&dt_ss, x0_dt, u_dt, 10, x_dt, y_dt);

    printf("  x[10] after 10 unit steps = %.1f (expected 10.0)\n", x_dt[10]);
    printf("  Observable: %s\n",
           dt_ss_is_observable(&dt_ss) ? "YES" : "NO");

    dt_ss_free(&dt_ss);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
