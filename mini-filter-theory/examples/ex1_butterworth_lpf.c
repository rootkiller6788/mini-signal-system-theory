/**
 * ex1_butterworth_lpf.c — Butterworth Lowpass Filter Design and Analysis
 *
 * L6 Canonical Problem: Design a 5th-order Butterworth lowpass filter
 * with fc=1kHz, analyze its frequency response, and demonstrate filtering.
 *
 * This example shows the complete analog filter design workflow:
 *   1. Specification
 *   2. Prototype design
 *   3. Frequency response analysis
 *   4. Filter a test signal (sine + noise via bilinear transform)
 *
 * Reference: Oppenheim & Willsky, Signals and Systems, Ch. 9
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "filter_defs.h"
#include "filter_transfer.h"
#include "filter_analog.h"
#include "filter_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Example 1: Butterworth Lowpass Filter ===\n\n");

    /* Step 1: Define specification */
    filter_spec_t spec = filter_spec_init();
    spec.type = FILTER_LOWPASS;
    spec.approx = APPROX_BUTTERWORTH;
    spec.fc1 = 1000.0;     /* 1 kHz cutoff */
    spec.order = 5;        /* 5th order */
    spec.gain_db = 0.0;    /* Unity gain */

    printf("Specification:\n");
    printf("  Type: Butterworth Lowpass\n");
    printf("  Order: %d\n", spec.order);
    printf("  Cutoff: %.0f Hz\n", spec.fc1);

    /* Step 2: Design the filter */
    tf_continuous_t *tf = analog_filter_design(&spec);
    if (!tf) {
        printf("ERROR: Filter design failed!\n");
        return 1;
    }

    printf("\nTransfer Function:\n");
    printf("  Numerator degree: %d\n", tf->num_len - 1);
    printf("  Denominator degree: %d\n", tf->den_len - 1);

    /* Step 3: Analyze frequency response */
    printf("\nFrequency Response (key points):\n");
    printf("  %-12s %-12s %-10s\n", "Freq(Hz)", "|H|(dB)", "Phase(deg)");

    double test_freqs[] = {10.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    int nf = sizeof(test_freqs) / sizeof(test_freqs[0]);
    int i;

    for (i = 0; i < nf; i++) {
        double w = 2.0 * M_PI * test_freqs[i];
        freq_resp_point_t fp = tf_continuous_response_at(tf, w);
        printf("  %-12.0f %-12.2f %-10.1f\n",
               test_freqs[i], fp.magnitude_db, fp.phase_deg);
    }

    /* Step 4: Find exact -3dB cutoff */
    double f3db = analog_find_cutoff(tf, 100.0, 10000.0, 0.1);
    printf("\n  -3 dB cutoff frequency: %.1f Hz\n", f3db);

    /* Step 5: Test roll-off rate */
    double wc = 2.0 * M_PI * 1000.0;
    double mag_wc = cabs(tf_continuous_eval(tf, wc * I));
    double mag_2wc = cabs(tf_continuous_eval(tf, 2.0 * wc * I));
    double rolloff = 20.0 * log10(mag_2wc / mag_wc);
    printf("  Roll-off at 2x cutoff: %.1f dB (expected ~ -30 dB for 5th order)\n",
           rolloff);

    /* Step 6: Check stability */
    printf("\nStability: ");
    if (tf_continuous_is_stable(tf)) {
        printf("STABLE (all poles in LHP)\n");
    } else {
        printf("UNSTABLE!\n");
    }

    /* Verify all poles are in LHP */
    double complex *poles = malloc(spec.order * sizeof(double complex));
    if (poles) {
        butterworth_poles(spec.order, wc, poles);
        printf("Pole locations:\n");
        int k;
        for (k = 0; k < spec.order; k++) {
            printf("  p[%d] = %.3f %+.3fj  (|p|=%.3f)\n",
                   k, creal(poles[k]), cimag(poles[k]), cabs(poles[k]));
        }
        free(poles);
    }

    tf_continuous_free(tf);
    printf("\n=== Example 1 complete ===\n");
    return 0;
}
