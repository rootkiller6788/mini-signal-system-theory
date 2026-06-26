/**
 * @file ex2_stability.c
 * @brief L6 Example: Stability Analysis with Routh-Hurwitz and Margins
 *
 * Demonstrates the Routh-Hurwitz criterion, BIBO stability tests,
 * and gain/phase margin computation for feedback control systems.
 *
 * Canonical Problem:
 *   Given open-loop L(s), determine closed-loop stability and margins.
 *
 * Course: MIT 6.003 Ch.11, Stanford EE102A Ch.9, Berkeley EE123 Ch.7
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/system_defs.h"
#include "../include/transfer_function.h"
#include "../include/stability.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Example 2: Stability Analysis ===\n\n");

    /* ── Routh-Hurwitz: Verify stability of s^3 + 2s^2 + 3s + 4 ── */
    printf("--- Routh-Hurwitz Test ---\n");
    printf("Polynomial: s^3 + 2s^2 + 3s + 4\n");

    double poly[] = {4.0, 3.0, 2.0, 1.0};  /* a0...a3 */
    int special;
    double rh[30];
    int rhp = routh_hurwitz(poly, 3, rh, &special);

    printf("  RHP poles: %d\n", rhp);
    printf("  Special case: %d (0=none, 1=zero 1st elem, 2=row of zeros)\n", special);
    printf("  Verdict: %s\n", (rhp == 0) ? "STABLE" : "UNSTABLE");

    /* ── BIBO Stability: Compare stable vs unstable ── */
    printf("\n--- BIBO Stability Comparison ---\n");

    /* Stable: H1(s) = 1/(s^2 + 2s + 2), poles = -1 +/- j */
    ct_transfer_function_t tf_stable = ct_tf_alloc(0, 2);
    tf_stable.num_coeffs[0] = 2.0;
    tf_stable.den_coeffs[0] = 2.0;
    tf_stable.den_coeffs[1] = 2.0;
    tf_stable.den_coeffs[2] = 1.0;

    int stable1 = ct_is_bibo_stable(&tf_stable);
    printf("  System 1 (1/(s^2+2s+2)): %s\n",
           stable1 ? "BIBO STABLE" : "UNSTABLE");

    /* Unstable: H2(s) = 1/(s^2 - 2s + 2), poles = 1 +/- j (RHP) */
    ct_transfer_function_t tf_unstable = ct_tf_alloc(0, 2);
    tf_unstable.num_coeffs[0] = 2.0;
    tf_unstable.den_coeffs[0] = 2.0;
    tf_unstable.den_coeffs[1] = -2.0;
    tf_unstable.den_coeffs[2] = 1.0;

    int stable2 = ct_is_bibo_stable(&tf_unstable);
    printf("  System 2 (1/(s^2-2s+2)): %s\n",
           stable2 ? "BIBO STABLE" : "UNSTABLE");

    /* ── Gain and Phase Margins ── */
    printf("\n--- Stability Margins ---\n");
    printf("Open-loop: L(s) = 10/(s(s+1)(s+5))\n");

    /* Loop transfer function: L(s) = 10/(s^3 + 6s^2 + 5s) */
    ct_transfer_function_t L = ct_tf_alloc(0, 3);
    L.num_coeffs[0] = 10.0;
    L.den_coeffs[0] = 0.0;
    L.den_coeffs[1] = 5.0;
    L.den_coeffs[2] = 6.0;
    L.den_coeffs[3] = 1.0;

    frequency_response_t fr = freq_resp_alloc(200, 0.01, 20.0, 1);
    ct_frequency_response(&L, &fr);

    double gm_db, pm_deg, gm_freq, pm_freq;
    compute_stability_margins(&fr, &gm_db, &pm_deg, &gm_freq, &pm_freq);

    printf("  Gain margin: %.2f dB at %.3f Hz\n", gm_db, gm_freq);
    printf("  Phase margin: %.2f deg at %.3f Hz\n", pm_deg, pm_freq);

    if (gm_db > 0 && pm_deg > 0) {
        printf("  Verdict: CLOSED-LOOP STABLE (GM>0, PM>0)\n");
    } else {
        printf("  Verdict: MAY BE UNSTABLE (check Nyquist)\n");
    }

    /* ── Feedback Stability Check ── */
    printf("\n--- Feedback Stability ---\n");
    ct_transfer_function_t G = ct_tf_alloc(0, 1);
    G.num_coeffs[0] = 1.0;
    G.den_coeffs[0] = 0.0; G.den_coeffs[1] = 1.0;  /* G(s) = 1/s */

    ct_transfer_function_t H = ct_tf_alloc(0, 0);
    H.num_coeffs[0] = 1.0;
    H.den_coeffs[0] = 1.0;  /* H(s) = 1 (unity feedback) */

    int fb_stable = feedback_stability_check(&G, &H);
    printf("  Unity feedback of 1/s: %s\n",
           fb_stable ? "STABLE" : "UNSTABLE");
    printf("  (1/s has pole at origin, closed-loop = 1/(s+1) -> STABLE)\n");

    /* Cleanup */
    freq_resp_free(&fr);
    ct_tf_free(&tf_stable);
    ct_tf_free(&tf_unstable);
    ct_tf_free(&L);
    ct_tf_free(&G);
    ct_tf_free(&H);

    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
