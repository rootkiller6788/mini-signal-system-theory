/**
 * @example lorenz_attractor.c
 * @brief Lorenz attractor — canonical chaotic system (L6, L8).
 *
 * Lorenz equations (1963):
 *   dx/dt = σ(y - x)
 *   dy/dt = x(ρ - z) - y
 *   dz/dt = xy - βz
 *
 * Classic parameters for chaos: σ=10, ρ=28, β=8/3.
 *
 * The Lorenz system was discovered during meteorological research
 * and was one of the first recognized examples of deterministic chaos.
 *
 * Course: MIT 18.385, Tsinghua 非线性动力学
 * Reference: Lorenz, "Deterministic Nonperiodic Flow",
 *            J. Atmospheric Sciences, 1963.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/nonlinear_core.h"
#include "../include/nonlinear_stability.h"
#include "../include/nonlinear_bifurcation.h"
#include "../include/nonlinear_phase.h"
#include "../include/nonlinear_chaos.h"

static int lorenz_f(const double *x, double t, const double *params,
                    double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 3) return -1;
    double sigma = params[0];
    double rho   = params[1];
    double beta  = params[2];
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
    return 0;
}

int main(void)
{
    printf("=== Lorenz Attractor Analysis ===\n");
    printf("Equations:\n");
    printf("  dx/dt = sigma*(y - x)\n");
    printf("  dy/dt = x*(rho - z) - y\n");
    printf("  dz/dt = x*y - beta*z\n\n");

    /* Classic chaotic parameters */
    double sigma = 10.0, rho = 28.0, beta = 8.0/3.0;

    printf("Parameters: sigma=%.1f, rho=%.1f, beta=%.4f\n\n", sigma, rho, beta);

    nl_nonlinear_system_t sys;
    sys.dim = 3;
    sys.f = lorenz_f;
    sys.jacobian = NULL;
    sys.params[0] = sigma;
    sys.params[1] = rho;
    sys.params[2] = beta;
    sys.num_params = 3;
    sys.is_autonomous = 1;

    /* Find equilibria:
     * (1) Origin: (0,0,0)
     * (2,3) C±: (±sqrt(beta*(rho-1)), ±sqrt(beta*(rho-1)), rho-1)
     * For rho > 1: origin becomes unstable, C± are the "butterfly" centers */
    printf("--- Equilibria ---\n");

    double x0_arr[3];

    /* Origin */
    x0_arr[0] = 0.1; x0_arr[1] = 0.1; x0_arr[2] = 0.1;
    nl_state_t eq;
    int ret = nl_find_equilibrium(&sys, x0_arr, NULL, 0, &eq, 1e-6);
    if (ret == 0) {
        printf("Origin: (%.6f, %.6f, %.6f)\n",
               eq.x[0], eq.x[1], eq.x[2]);
        nl_jacobian_t J;
        nl_compute_jacobian(&sys, eq.x, &J);
        nl_stability_t stab = nl_classify_equilibrium(&J, 1e-6);
        printf("  Stability: ");
        switch (stab) {
        case NL_SADDLE: printf("SADDLE (1D unstable manifold + 2D stable)\n"); break;
        case NL_UNSTABLE: printf("UNSTABLE\n"); break;
        default: printf("type %d\n", (int)stab); break;
        }
    }

    /* C+ equilibrium */
    printf("\nC+: (%.4f, %.4f, %.1f) [analytic]\n",
           sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), rho-1);
    printf("C-: (%.4f, %.4f, %.1f) [analytic]\n",
           -sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1);

    /* Integrate trajectory from initial condition near C+ */
    printf("\n--- Trajectory ---\n");
    double x0[3] = {1.0, 1.0, 1.0};
    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 20000, 3);
    size_t n_steps;
    nl_integrate_rk4(&sys, x0, 50.0, 0.005, &traj, &n_steps);
    printf("Integrated %zu steps\n", n_steps);

    /* Extract x-coordinate time series for chaos analysis */
    double *ts = (double *)malloc(traj.num_points * sizeof(double));
    if (ts) {
        for (size_t i = 0; i < traj.num_points; i++)
            ts[i] = traj.points[i].x[0];

        nl_chaos_metrics_t metrics;
        nl_chaos_analyze(&sys, ts, traj.num_points, &metrics);

        printf("\n--- Chaos Metrics ---\n");
        printf("  Largest Lyapunov exponent: %.6f\n",
               metrics.lyapunov_exponent);
        printf("  Kaplan-Yorke dimension: %d\n",
               metrics.lyapunov_dimension);
        printf("  K-S entropy (Pesin): %.6f\n",
               metrics.ks_entropy);

        /* 0-1 test for chaos */
        double K;
        if (nl_test_chaos_01(ts, traj.num_points, 50, &K) == 0) {
            printf("  0-1 test for chaos: K = %.4f ", K);
            if (K > 0.5) printf("(CHAOTIC)\n");
            else printf("(REGULAR)\n");
        }

        /* The Lorenz attractor has:
         *   λ₁ ≈ 0.906 (positive → chaos)
         *   λ₂ ≈ 0.0
         *   λ₃ ≈ -14.57
         *   D_KY ≈ 2.062
         */
        printf("\n  Reference values (Sprott, 2003):\n");
        printf("    λ₁ ≈ 0.906, λ₂ ≈ 0.0, λ₃ ≈ -14.57\n");
        printf("    D_KY ≈ 2.062 (fractal)\n");

        free(ts);
    }

    nl_trajectory_free(&traj);
    printf("\n=== Analysis Complete ===\n");
    return 0;
}
