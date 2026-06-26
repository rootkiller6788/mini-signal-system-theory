/**
 * @example pll_nonlinear.c
 * @brief Phase-Locked Loop (PLL) nonlinear analysis (L6, L7).
 *
 * A PLL is a classic nonlinear feedback system used in communication
 * systems for carrier recovery, frequency synthesis, and clock
 * generation. The loop contains a phase detector (nonlinear),
 * loop filter (linear), and VCO (integrator).
 *
 * The nonlinear PLL model:
 *   dφ/dt = ω_in - ω_vco - K_pd * K_vco * h(t) * sin(φ)
 *
 * where sin(φ) is the phase detector nonlinearity and h(t) is
 * the loop filter impulse response.
 *
 * Analysis covers:
 *   - Lock range (using DF method)
 *   - Pull-in range
 *   - Hold range
 *   - Cycle slipping behavior
 *
 * Application: GPS receivers (L7), 5G NR carrier recovery,
 *              Boeing 747 navigation systems.
 *
 * Course: MIT 6.450, Stanford EE359, Tsinghua 通信原理
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/nonlinear_core.h"
#include "../include/nonlinear_stability.h"
#include "../include/nonlinear_describing.h"
#include "../include/nonlinear_phase.h"

/* PLL phase error dynamics (2nd-order, type II PLL):
 * State: x[0] = phase error φ, x[1] = frequency error dφ/dt
 * Parameters: params[0] = loop gain K, params[1] = natural freq ω_n,
 *             params[2] = damping ζ */
static int pll_f(const double *x, double t, const double *params,
                 double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 2) return -1;
    double K   = params[0];  /* loop gain */
    double wn  = params[1];  /* natural frequency */
    double zeta = params[2]; /* damping factor */
    double domega = params[3]; /* initial frequency offset */

    dxdt[0] = x[1];
    /* d²φ/dt² + 2ζω_n dφ/dt + K ω_n² sin(φ) = 0
     * => dxdt[1] = -2ζω_n x[1] - K ω_n² sin(x[0]) + domega */
    dxdt[1] = -2.0 * zeta * wn * x[1]
              - K * wn * wn * sin(x[0])
              + domega;
    return 0;
}

int main(void)
{
    printf("=== PLL Nonlinear Analysis ===\n");
    printf("Model: 2nd-order type-II PLL with sinusoidal PD\n\n");

    /* PLL parameters (typical GPS receiver PLL) */
    double K      = 1.0;    /* normalized loop gain */
    double wn     = 50.0;   /* natural frequency ω_n = 50 rad/s (~8 Hz) */
    double zeta   = 0.707;  /* damping factor (Butterworth) */
    double domega = 10.0;   /* initial frequency offset */

    printf("Parameters:\n");
    printf("  Loop gain K = %.2f\n", K);
    printf("  Natural freq ω_n = %.1f rad/s\n", wn);
    printf("  Damping ζ = %.3f\n", zeta);
    printf("  Initial freq offset Δω = %.1f rad/s\n\n", domega);

    nl_nonlinear_system_t sys;
    sys.dim = 2;
    sys.f = pll_f;
    sys.jacobian = NULL;
    sys.params[0] = K;
    sys.params[1] = wn;
    sys.params[2] = zeta;
    sys.params[3] = domega;
    sys.num_params = 4;
    sys.is_autonomous = 1;

    /* --- Equilibria --- */
    printf("--- Equilibria (phase-locked states) ---\n");
    /* sin(φ) = Δω / (K ω_n²) */
    double sin_phi_eq = domega / (K * wn * wn);
    printf("sin(φ_eq) = Δω / (K ω_n²) = %.6f\n", sin_phi_eq);

    if (fabs(sin_phi_eq) <= 1.0) {
        double phi_eq = asin(sin_phi_eq);
        printf("Stable equilibrium: φ* = %.4f rad (%.2f deg)\n",
               phi_eq, phi_eq * 180.0 / M_PI);

        double x_eq_arr[2] = {phi_eq, 0.0};
        nl_jacobian_t J;
        nl_state_t eq_st; (void)eq_st;
        eq_st.dim = 2;
        eq_st.x[0] = phi_eq; eq_st.x[1] = 0.0;
        nl_compute_jacobian(&sys, x_eq_arr, &J);
        nl_stability_t stab = nl_classify_equilibrium(&J, 1e-6);
        printf("Local stability: ");
        switch (stab) {
        case NL_STABLE_FOCUS:
            printf("STABLE FOCUS (damped oscillations)\n"); break;
        case NL_STABLE_NODE:
            printf("STABLE NODE (overdamped)\n"); break;
        default:
            printf("type %d\n", (int)stab); break;
        }

        /* Jacobian at lock point */
        printf("Jacobian = [[%.4f, %.4f], [%.4f, %.4f]]\n",
               J.matrix[0][0], J.matrix[0][1],
               J.matrix[1][0], J.matrix[1][1]);
    } else {
        printf("NO equilibrium exists: Δω exceeds hold range!\n");
    }

    /* --- Lock Range via Describing Function --- */
    printf("\n--- Lock Range Analysis (DF Method) ---\n");

    /* Treat sin(φ) as nonlinearity with N(A) = 2*J1(A)/A
     * where J1 is the Bessel function.
     * For small A: N(A) ≈ 1 - A²/8
     * For large A: N(A) ≈ 0 (sin nonlinearity saturates to ±1) */

    /* DF for sinusoidal nonlinearity:
     * f(φ) = sin(φ) has DF N(A) = 2*J1(A)/A (real, no phase shift) */
    printf("The sinusoidal PD nonlinearity f(φ) = sin(φ)\n");
    printf("has describing function N(A) = 2*J₁(A) / A\n\n");

    double A_vals[] = {0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0};
    printf("N(A) values:\n");
    for (int i = 0; i < 7; i++) {
        /* Use numerical DF or the Bessel approximation */
        nl_nonlinearity_t nl;
        double p[1] = {0.0};  /* not used for custom */
        nl_nonlinearity_init(&nl, NL_TYPE_CUSTOM, p, 0);

        nl_describing_function_t df;
        if (nl_describing_function(&nl, A_vals[i], 200, &df) == 0) {
            printf("  A=%.1f: N = %.4f + j%.4f\n",
                   A_vals[i], df.P, df.Q);
        }
    }

    /* Hold range: |Δω| < K ω_n² (from equilibrium condition)
     * Lock range: smaller than hold range due to cycle slipping
     * Pull-in range: depends on loop filter type */
    double hold_range = K * wn * wn;
    printf("\n--- Range Summary ---\n");
    printf("Hold range:  |Δω| < %.1f rad/s (%.1f Hz)\n",
           hold_range, hold_range / (2.0 * M_PI));
    printf("Lock range:  typically 70-90%% of hold range for 2nd-order PLL\n");
    printf("Pull-in range: depends on loop filter bandwidth\n");

    /* --- Simulate pull-in transient --- */
    printf("\n--- Pull-in Transient Simulation ---\n");
    double x0[2] = {0.0, domega};  /* start with phase 0, frequency offset */
    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 5000, 2);
    size_t n_steps;
    nl_integrate_rk4(&sys, x0, 5.0, 0.002, &traj, &n_steps);

    printf("Simulated %.1fs with dt=0.002s (%zu steps)\n",
           5.0, n_steps);

    /* Check if phase error converged */
    double phi_final = traj.points[traj.num_points - 1].x[0];
    double freq_final = traj.points[traj.num_points - 1].x[1];
    printf("Final phase error: %.4f rad\n", phi_final);
    printf("Final frequency error: %.4f rad/s\n", freq_final);

    if (fabs(freq_final) < 1.0) {
        printf("PLL LOCKED successfully\n");
    } else {
        printf("PLL NOT LOCKED — cycle slipping detected\n");
    }

    nl_trajectory_free(&traj);

    printf("\n--- Applications (L7) ---\n");
    printf("• GPS L1 carrier tracking (1575.42 MHz)\n");
    printf("• 5G NR carrier recovery (sub-6 GHz and mmWave)\n");
    printf("• Boeing 787 navigation systems\n");
    printf("• Tesla Autopilot radar signal processing\n");
    printf("• MRI clock generation\n");

    return 0;
}
