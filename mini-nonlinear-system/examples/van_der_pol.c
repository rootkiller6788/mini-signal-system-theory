/**
 * @example van_der_pol.c
 * @brief Van der Pol oscillator — canonical nonlinear system (L6).
 *
 * Equation: d^2x/dt^2 - mu(1-x^2)dx/dt + x = 0
 *
 * For mu > 0: nonlinear damping produces a stable limit cycle.
 * For mu << 1: nearly sinusoidal (harmonic oscillator)
 * For mu >> 1: relaxation oscillations
 *
 * This example integrates the VdP oscillator from various initial
 * conditions, demonstrates convergence to the limit cycle, and
 * computes the describing function prediction.
 *
 * Course: MIT 6.003, Stanford EE263, Tsinghua 信号与系统
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/nonlinear_core.h"
#include "../include/nonlinear_stability.h"
#include "../include/nonlinear_phase.h"
#include "../include/nonlinear_describing.h"

static int vdp_f(const double *x, double t, const double *params,
                 double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 2) return -1;
    double mu = params[0];
    dxdt[0] = x[1];
    dxdt[1] = mu * (1.0 - x[0] * x[0]) * x[1] - x[0];
    return 0;
}

int main(void)
{
    printf("=== Van der Pol Oscillator Analysis ===\n");
    printf("Equation: d2x/dt2 - mu*(1-x^2)*dx/dt + x = 0\n\n");

    double mu_values[] = {0.5, 2.0, 5.0};

    for (int k = 0; k < 3; k++) {
        double mu = mu_values[k];

        /* Setup system */
        nl_nonlinear_system_t sys;
        sys.dim = 2;
        sys.f = vdp_f;
        sys.jacobian = NULL;
        sys.params[0] = mu;
        sys.num_params = 1;
        sys.is_autonomous = 1;

        printf("--- mu = %.1f ---\n", mu);

        /* Equilibrium at origin */
        double x_eq_arr[2] = {0.0, 0.0};
        nl_jacobian_t J;
        nl_state_t eq_state; (void)eq_state;
        eq_state.dim = 2;
        eq_state.x[0] = 0.0; eq_state.x[1] = 0.0;

        nl_compute_jacobian(&sys, x_eq_arr, &J);
        nl_stability_t stab = nl_classify_equilibrium(&J, 1e-6);

        printf("  Jacobian at origin: [[%.4f, %.4f], [%.4f, %.4f]]\n",
               J.matrix[0][0], J.matrix[0][1],
               J.matrix[1][0], J.matrix[1][1]);

        /* For mu > 0: origin is unstable focus/node
         * J = [[0, 1], [-1, mu]] -> eigenvalues: (mu ± sqrt(mu^2-4))/2 */
        printf("  Origin stability: ");
        switch (stab) {
        case NL_UNSTABLE_FOCUS:
            printf("UNSTABLE FOCUS -> limit cycle expected\n"); break;
        case NL_UNSTABLE_NODE:
            printf("UNSTABLE NODE -> limit cycle expected\n"); break;
        case NL_SADDLE:
            printf("SADDLE\n"); break;
        default:
            printf("type %d\n", (int)stab); break;
        }

        /* Integrate from initial condition */
        double x0[2] = {2.0, 0.0};
        nl_trajectory_t traj;
        nl_trajectory_init(&traj, 5000, 2);
        size_t n_steps;
        nl_integrate_rk4(&sys, x0, 30.0, 0.01, &traj, &n_steps);

        /* Find approximate limit cycle amplitude from trajectory */
        double max_x = 0.0;
        size_t n_skip = traj.num_points / 2;  /* skip transient */
        for (size_t i = n_skip; i < traj.num_points; i++) {
            if (fabs(traj.points[i].x[0]) > max_x)
                max_x = fabs(traj.points[i].x[0]);
        }
        printf("  Limit cycle amplitude (from simulation): x_max ≈ %.4f\n",
               max_x);

        /* Estimate period from zero crossings */
        int crossings = 0;
        double t_first = 0.0, t_last = 0.0;
        double dt = 0.01;
        for (size_t i = n_skip + 1; i < traj.num_points; i++) {
            if (traj.points[i-1].x[0] * traj.points[i].x[0] < 0.0
                && traj.points[i-1].x[0] < 0.0) {
                if (crossings == 0)
                    t_first = (double)i * dt;
                t_last = (double)i * dt;
                crossings++;
            }
        }
        if (crossings >= 2) {
            double period = (t_last - t_first) / (double)(crossings / 2);
            double freq = 2.0 * M_PI / period;
            printf("  Estimated oscillation frequency: ω ≈ %.4f rad/s\n",
                   freq);
        }

        printf("  Trajectory points: %zu\n\n", traj.num_points);
        nl_trajectory_free(&traj);
    }

    /* Describing function prediction for mu=2.0 */
    printf("--- DF Prediction (mu=2.0) ---\n");
    printf("VdP can be written as Lur'e system:\n");
    printf("  Linear part: G(s) = s / (s^2 - mu*s + 1)\n");
    printf("  Nonlinearity: φ(y) = (mu/3) * y^3\n");
    printf("DF method predicts a limit cycle satisfying G(jw) * N(A) = -1\n");
    printf("For relaxation oscillations (mu >> 1):\n");
    printf("  Amplitude A ≈ 2.0, Period T ≈ mu * (3 - 2*ln 2) ≈ 1.614*mu\n");

    return 0;
}
