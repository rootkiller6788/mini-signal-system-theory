/**
 * @file ex1_step_response.c
 * @brief L6 Example: Step Response of First and Second-Order Systems
 *
 * Demonstrates computation of step response from transfer functions
 * and analysis of rise time, overshoot, and settling time.
 *
 * Canonical Problem:
 *   Given H(s), compute y(t) = L^{-1}{H(s)/s} and extract time-domain metrics.
 *
 * Course: MIT 6.003 Ch.3, Stanford EE102A Ch.5, Berkeley EE123 Ch.4
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/system_defs.h"
#include "../include/convolution.h"
#include "../include/transfer_function.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Example 1: Step Response Analysis ===\n\n");

    /* ── First-Order System: H(s) = 1/(tau*s + 1) ── */
    printf("--- First-Order System: H(s) = 1/(s + 1) ---\n");
    printf("Time constant tau = 1 sec\n");

    ct_transfer_function_t tf1 = ct_tf_alloc(0, 1);
    tf1.num_coeffs[0] = 1.0;
    tf1.den_coeffs[0] = 1.0;
    tf1.den_coeffs[1] = 1.0;

    /* Simulate step response via convolution with step input */
    int N = 1000;
    double T_end = 10.0;
    double dt = T_end / (double)(N - 1);
    double Fs = 1.0 / dt;

    /* Generate impulse response h(t) via inverse Laplace (sampled) */
    /* h(t) = (1/tau) * exp(-t/tau) = exp(-t) */
    ct_signal_t h_ct = ct_signal_alloc(N, 0.0, T_end, Fs);
    dt_signal_t h_dt = dt_signal_alloc(N);
    for (int i = 0; i < N; i++) {
        double t = (double)i * dt;
        h_ct.samples[i] = exp(-t);  /* h(t) = e^{-t} */
        h_dt.samples[i] = h_ct.samples[i];
    }

    /* Step input u[n] = 1 for n >= 0 */
    dt_signal_t u = dt_signal_alloc(N);
    dt_signal_set_unit_step(&u, 0);

    /* Output: y = u * h */
    dt_signal_t y = dt_signal_alloc(2 * N);
    dt_convolution_direct(&u, &h_dt, &y);

    printf("  t=1.0: y = %.4f (expected 0.6321)\n", y.samples[100]);
    printf("  t=3.0: y = %.4f (expected 0.9502)\n", y.samples[300]);
    printf("  t=5.0: y = %.4f (expected 0.9933)\n", y.samples[500]);

    /* Step response characteristics from impulse: s(t) = integral_0^t h(tau) dtau */
    dt_signal_t s = dt_signal_alloc(N);
    dt_step_from_impulse(&h_dt, &s);
    printf("  Step response s[100] (t=1.0): %.4f\n", s.samples[100]);
    printf("  Steady-state (t=10): %.4f (DC gain = 1.0)\n",
           s.samples[N-1]);

    dt_signal_free(&u); dt_signal_free(&y); dt_signal_free(&h_dt);
    dt_signal_free(&s);
    ct_signal_free(&h_ct);
    ct_tf_free(&tf1);

    /* ── Second-Order System: H(s) = wn^2/(s^2 + 2*zeta*wn*s + wn^2) ── */
    printf("\n--- Second-Order System: zeta=0.3, wn=2 rad/s ---\n");
    double zeta = 0.3, wn = 2.0;

    ct_transfer_function_t tf2 = ct_tf_alloc(0, 2);
    tf2.num_coeffs[0] = wn * wn;
    tf2.den_coeffs[0] = wn * wn;
    tf2.den_coeffs[1] = 2.0 * zeta * wn;
    tf2.den_coeffs[2] = 1.0;

    /* Evaluate frequency response */
    frequency_response_t fr = freq_resp_alloc(100, 0.01, 10.0, 1);
    ct_frequency_response(&tf2, &fr);

    printf("  DC gain: %.4f (%.2f dB)\n", fr.points[0].magnitude,
           fr.points[0].magnitude_db);
    printf("  Peak magnitude: %.4f at %.3f Hz\n",
           fr.points[50].magnitude,
           fr.points[50].frequency);

    /* Resonance peak: M_r = 1/(2*zeta*sqrt(1-zeta^2)) for zeta < 1/sqrt(2) */
    double Mr_expected = 1.0 / (2.0 * zeta * sqrt(1.0 - zeta * zeta));
    printf("  Expected resonance peak: %.4f\n", Mr_expected);

    /* Settling time estimate (2%%): ts = 4/(zeta*wn) */
    double ts = 4.0 / (zeta * wn);
    printf("  Settling time (2%% criterion): %.3f sec\n", ts);

    /* Peak overshoot: Mp = exp(-pi*zeta/sqrt(1-zeta^2)) * 100%% */
    double Mp = 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
    printf("  Peak overshoot: %.2f%%\n", Mp);

    freq_resp_free(&fr);
    ct_tf_free(&tf2);

    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
