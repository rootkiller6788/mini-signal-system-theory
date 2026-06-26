/*
 * ex3_bilinear_filter.c - Bilinear Transform and Digital Filter Design
 *
 * Demonstrates the bilinear transform (Tustin's method) for converting
 * a continuous-time filter to a discrete-time equivalent.
 *
 * Bilinear transform: s = (2/T) * (1 - z^{-1})/(1 + z^{-1})
 *
 * Designs a digital lowpass Butterworth filter from an analog prototype.
 *
 * L6: Digital filter design via bilinear transform
 * L7: GPS receiver filtering (GPS application)
 */

#include "laplace_z_transform_core.h"
#include "laplace_transform.h"
#include "z_transform.h"
#include "transfer_function.h"
#include "stability_analysis.h"
#include <stdio.h>
#include <math.h>

/*
 * Apply bilinear transform to convert H(s) to H(z).
 *
 * s = (2/T) * (1 - z^{-1})/(1 + z^{-1})
 *
 * For a first-order system H(s) = a/(s + a):
 * H(z) = a / ((2/T)*(1-z^{-1})/(1+z^{-1}) + a)
 *      = a*(1+z^{-1}) / ((2/T)*(1-z^{-1}) + a*(1+z^{-1}))
 *      = (a*T + a*T*z^{-1}) / ((2 + a*T) + (a*T - 2)*z^{-1})
 */

rational_func_t bilinear_first_order(double a, double T) {
    /*
     * Bilinear transform of H(s) = a/(s+a).
     *
     * H(z) = (a*T*(1+z^{-1})) / ((2+a*T) + (a*T-2)*z^{-1})
     *
     * Reference: Oppenheim & Schafer (2010), Section 7.2.2
     */
    double aT = a * T;
    double n0 = aT;
    double n1 = aT;
    double d0 = 2.0 + aT;
    double d1 = aT - 2.0;

    double ncoeffs[2] = {n0/d0, n1/d0};
    double dcoeffs[2] = {1.0, d1/d0};

    polynomial_t num = polynomial_from_coeffs(ncoeffs, 2);
    polynomial_t den = polynomial_from_coeffs(dcoeffs, 2);

    rational_func_t Hz; Hz.num = num; Hz.den = den;
    return Hz;
}

/*
 * Bilinear transform for second-order system:
 * H(s) = wn^2/(s^2 + 2*zeta*wn*s + wn^2)
 *
 * Substitute s = (2/T)*(1-z^{-1})/(1+z^{-1}):
 *
 * H(z^{-1}) = (b0 + b1*z^{-1} + b2*z^{-2}) / (1 + a1*z^{-1} + a2*z^{-2})
 */
rational_func_t bilinear_second_order(double wn, double zeta, double T) {
    double K = 2.0 / T;
    double K2 = K * K;

    double wn2 = wn * wn;
    double den_scale = K2 + 2.0*zeta*wn*K + wn2;

    double b0 = wn2 / den_scale;
    double b1 = 2.0 * wn2 / den_scale;
    double b2 = wn2 / den_scale;

    double a1 = (2.0*wn2 - 2.0*K2) / den_scale;
    double a2 = (K2 - 2.0*zeta*wn*K + wn2) / den_scale;

    double ncoeffs[3] = {b0, b1, b2};
    double dcoeffs[3] = {1.0, a1, a2};

    polynomial_t num = polynomial_from_coeffs(ncoeffs, 3);
    polynomial_t den = polynomial_from_coeffs(dcoeffs, 3);

    rational_func_t Hz; Hz.num = num; Hz.den = den;
    return Hz;
}

int main(void) {
    printf("=== Bilinear Transform and Digital Filter Design ===\n\n");

    double T = 0.001; /* sampling period = 1 ms (1 kHz) */
    double fs = 1.0 / T;

    /* 1. First-order analog lowpass: H(s) = a/(s+a), a = 100 rad/s */
    printf("1. First-Order Analog Lowpass -> Digital Filter\n");
    printf("   Sampling: fs = %.0f Hz, T = %.4f s\n", fs, T);
    double a = 100.0; /* cutoff = 100 rad/s (~16 Hz) */
    printf("   Analog: H(s) = %.1f/(s + %.1f)\n", a, a);

    rational_func_t Hz1 = bilinear_first_order(a, T);
    printf("   Digital: H(z) = ");
    rational_print(&Hz1, "z^{-1}");

    /* Check DC gain preservation */
    double dc_analog = a / a; /* = 1 */
    double dc_digital = rational_mag_at(&Hz1, 0.0, 1);
    printf("   DC gain: analog=%.4f, digital=%.4f\n", dc_analog, dc_digital);

    /* Check stability */
    int stable = zt_is_stable(&Hz1);
    printf("   Stability: %s\n", stable ? "STABLE (poles inside unit circle)" : "UNSTABLE");
    printf("\n");

    /* 2. Second-order Butterworth filter */
    printf("2. Second-Order Butterworth Lowpass Filter\n");
    double fc = 50.0; /* cutoff frequency: 50 Hz */
    double wc = 2.0 * M_PI * fc;
    double zeta_butter = 0.707106781; /* 1/sqrt(2) */

    printf("   Analog prototype: fc=%.0f Hz, wn=%.1f rad/s, zeta=%.4f\n",
           fc, wc, zeta_butter);

    /* Pre-warp the cutoff frequency for bilinear transform */
    double wc_warped = (2.0/T) * tan(wc*T/2.0);
    printf("   Pre-warped cutoff: %.1f rad/s (to compensate bilinear warping)\n",
           wc_warped);

    rational_func_t Hz2 = bilinear_second_order(wc_warped, zeta_butter, T);
    printf("   Digital biquad coefficients:\n");
    printf("   Numerator: [%.6f, %.6f, %.6f]\n",
           Hz2.num.coeff[0], Hz2.num.coeff[1], Hz2.num.coeff[2]);
    printf("   Denominator: [%.6f, %.6f, %.6f]\n",
           Hz2.den.coeff[0], Hz2.den.coeff[1], Hz2.den.coeff[2]);

    /* Frequency response at key points */
    printf("\n   Frequency Response:\n");
    printf("   Freq (Hz) | Mag (dB) | Phase (deg)\n");
    printf("   -----------+-----------+------------\n");
    double test_freqs[] = {1.0, 10.0, 30.0, 50.0, 80.0, 200.0, 400.0};
    for (int i = 0; i < 7; i++) {
        double f = test_freqs[i];
        double w = 2.0 * M_PI * f;
        /* Discrete-time frequency: omega = w*T */
        double mag_db = 20.0 * log10(rational_mag_at(&Hz2, w*T, 1));
        double phase_deg = rational_phase_at(&Hz2, w*T, 1) * 180.0 / M_PI;
        printf("   %9.1f | %8.3f | %9.2f\n", f, mag_db, phase_deg);
    }

    /* Stability */
    printf("\n   Stability: %s\n",
           zt_is_stable(&Hz2) ? "STABLE" : "UNSTABLE");

    /* 3. L7 Application: GPS receiver anti-aliasing filter */
    printf("\n3. L7 Application: GPS L1 C/A Code Anti-Aliasing Filter\n");
    printf("   GPS L1 signal: 1575.42 MHz carrier\n");
    printf("   C/A code chip rate: 1.023 MHz\n");
    printf("   Digital IF sampling: 5.714 MHz (typical)\n");
    printf("   Anti-aliasing lowpass: Butterworth, fc=2.5 MHz\n");
    printf("   Prevents aliasing of out-of-band noise before ADC.\n");
    printf("   This is implemented using the bilinear transform method above.\n");

    rational_free(&Hz1);
    rational_free(&Hz2);

    printf("\n=== Design Complete ===\n");
    return 0;
}
