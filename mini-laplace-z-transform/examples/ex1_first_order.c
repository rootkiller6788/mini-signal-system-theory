/*
 * ex1_first_order.c - First-Order System Analysis
 *
 * Analyzes a first-order RC lowpass filter:
 *   H(s) = 1/(tau*s + 1) = 1/(RC*s + 1)
 *
 * Demonstrates:
 *   - Laplace transform of first-order system
 *   - Step response computation
 *   - Frequency response (Bode plot data)
 *   - Settling time estimation
 *   - Pole-zero analysis
 *
 * L6: Canonical first-order system analysis
 */

#include "laplace_z_transform_core.h"
#include "laplace_transform.h"
#include "transfer_function.h"
#include "stability_analysis.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== First-Order System Analysis ===\n");
    printf("System: H(s) = 1/(tau*s + 1), tau = 1.0\n\n");

    double tau = 1.0;
    rational_func_t H = laplace_first_order(1.0, tau);

    /* 1. Transfer function display */
    printf("1. Transfer Function:\n");
    rational_print(&H, "s");

    /* 2. DC gain */
    double dc_gain = tf_dc_gain(&H, 0);
    printf("\n2. DC Gain: %.6g (expected: 1.0)\n", dc_gain);

    /* 3. Pole location */
    pole_zero_t pz = pz_from_rational(&H);
    printf("\n3. Pole-Zero Analysis:\n");
    pz_print(&pz);

    /* 4. Step response at key time points */
    printf("\n4. Step Response y(t) = 1 - exp(-t/tau):\n");
    double times[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
    for (int i = 0; i < 7; i++) {
        double y = tf_step_response(&H, times[i]);
        double expected = 1.0 - exp(-times[i] / tau);
        printf("  y(%.1f) = %.6f (expected: %.6f, err: %.2e)\n",
               times[i], y, expected, fabs(y - expected));
    }

    /* 5. Settling time (2% criterion): ts = 4*tau */
    double ts = tf_settling_time(&H);
    printf("\n5. Settling Time (2%%): %.4f (theoretical: %.4f)\n", ts, 4.0*tau);

    /* 6. Frequency response */
    printf("\n6. Frequency Response (Bode):\n");
    printf("  omega (rad/s) |  |H(jw)|   |  Phase (deg)\n");
    printf("  ---------------+-------------+--------------\n");
    double freqs[] = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    for (int i = 0; i < 6; i++) {
        double mag_db = tf_bode_mag_db(&H, freqs[i]);
        double phase_deg = tf_bode_phase_deg(&H, freqs[i]);
        printf("  %12.4f | %9.4f dB | %9.2f deg\n", freqs[i], mag_db, phase_deg);
    }

    /* 7. Cutoff frequency (-3dB) */
    double wc = tf_cutoff_frequency(&H, 0.01, 10.0);
    printf("\n7. -3dB Cutoff Frequency: %.4f rad/s (expected: %.4f)\n", wc, 1.0/tau);

    /* 8. Stability check */
    printf("\n8. Stability: %s\n", laplace_is_stable(&H) ? "STABLE" : "UNSTABLE");

    /* 9. DC motor speed model (L7 application: Detroit, Toyota) */
    printf("\n9. Application: DC Motor Speed Model\n");
    printf("   A DC motor speed can be modeled as 1st-order:\n");
    printf("   H(s) = K/(tau*s + 1) where tau = J/B (inertia/damping)\n");
    double J = 0.01, B = 0.1; /* inertia and damping */
    double tau_motor = J / B;
    printf("   J=%.2f kg.m^2, B=%.2f N.m.s, tau=%.2f s\n", J, B, tau_motor);

    rational_free(&H);
    pz_free(&pz);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
