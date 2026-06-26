/*
 * ex2_second_order.c - Second-Order System Analysis
 *
 * Analyzes a canonical second-order system:
 *   H(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
 *
 * Demonstrates:
 *   - Damping ratio and natural frequency extraction
 *   - Overshoot computation
 *   - Step response with varied damping
 *   - Pole-zero geometry
 *   - Bode plot analysis
 *
 * L6: Canonical second-order system analysis
 * L7: Automotive suspension (Toyota) and quadrotor dynamics
 */

#include "laplace_z_transform_core.h"
#include "laplace_transform.h"
#include "transfer_function.h"
#include "stability_analysis.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Second-Order System Analysis ===\n\n");

    /* Analyze systems with different damping ratios */
    double zetas[] = {0.1, 0.3, 0.5, 0.707, 1.0, 1.5};
    double wn = 2.0; /* natural frequency = 2 rad/s */
    const char *desc[] = {
        "Highly underdamped",
        "Underdamped",
        "Moderately damped",
        "Butterworth (max flat)",
        "Critically damped",
        "Overdamped"
    };

    for (int i = 0; i < 6; i++) {
        double zeta = zetas[i];
        printf("--- Case %d: zeta=%.3f (%s) ---\n", i+1, zeta, desc[i]);

        rational_func_t H = laplace_second_order(1.0, wn, zeta);

        /* Damping parameters */
        double zeta_calc, wn_calc;
        if (tf_damping_params(&H, &zeta_calc, &wn_calc) == 0) {
            printf("  Extracted: zeta=%.4f, wn=%.4f rad/s\n", zeta_calc, wn_calc);
        }

        /* Overshoot */
        double Mp = tf_overshoot_percent(&H);
        printf("  Overshoot: %.2f%%", Mp);
        if (zeta < 1.0 && zeta > 0.0) {
            double expected = 100.0 * exp(-M_PI * zeta / sqrt(1 - zeta*zeta));
            printf(" (theoretical: %.2f%%)", expected);
        }
        printf("\n");

        /* Settling time */
        double ts = tf_settling_time(&H);
        printf("  Settling time (2%%): %.4f s", ts);
        if (zeta > 0.0) {
            printf(" (theoretical: %.4f s)", 4.0/(zeta*wn));
        }
        printf("\n");

        /* Step response at a few points */
        printf("  Step response: ");
        for (double t = 0.5; t <= 4.0; t += 0.5) {
            printf("y(%.1f)=%.4f ", t, tf_step_response(&H, t));
        }
        printf("\n\n");

        rational_free(&H);
    }

    /* L7 Application: Quadrotor altitude dynamics */
    printf("=== L7 Application: Quadrotor Altitude Dynamics ===\n");
    printf("A quadrotor's vertical motion can be modeled as 2nd-order:\n");
    printf("  H(s) = K/(m*s^2 + b*s + k)\n");
    double m = 1.5, b = 0.3, k = 10.0; /* mass, damping, spring constant */
    double wn_q = sqrt(k/m);
    double zeta_q = b/(2*sqrt(m*k));
    printf("  m=%.1f kg, b=%.1f N.s/m, k=%.1f N/m\n", m, b, k);
    printf("  wn=%.3f rad/s, zeta=%.4f\n", wn_q, zeta_q);
    rational_func_t H_quad = laplace_second_order(1.0/k, wn_q, zeta_q);
    printf("  Settling time: %.3f s\n", tf_settling_time(&H_quad));
    printf("  Overshoot: %.2f%%\n", tf_overshoot_percent(&H_quad));
    rational_free(&H_quad);

    /* L7 Application: Automotive suspension (Toyota) */
    printf("\n=== L7 Application: Automotive Suspension (Quarter-Car Model) ===\n");
    printf("A quarter-car suspension is a 2nd-order mass-spring-damper:\n");
    double ms = 300.0;  /* sprung mass (kg) */
    double ks = 20000.0; /* spring stiffness (N/m) */
    double cs = 1500.0;  /* damping coefficient (N.s/m) */
    double wn_car = sqrt(ks/ms);
    double zeta_car = cs/(2*sqrt(ms*ks));
    printf("  ms=%.0f kg, ks=%.0f N/m, cs=%.0f N.s/m\n", ms, ks, cs);
    printf("  wn=%.3f rad/s (%.2f Hz), zeta=%.4f\n",
           wn_car, wn_car/(2*M_PI), zeta_car);
    printf("  This is typical for passenger vehicles (Toyota Camry class).\n");

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
