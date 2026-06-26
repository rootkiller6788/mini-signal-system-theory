/**
 * demo_response.c — Frequency Response Visualization Demo
 *
 * Demonstrates filter magnitude and phase response for a 5th-order
 * Butterworth lowpass, outputting CSV for plotting.
 *
 * Usage: ./build/demo_response > response.csv
 * Plot with: gnuplot, matplotlib, Excel, etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "filter_defs.h"
#include "filter_transfer.h"
#include "filter_analog.h"
#include "filter_analysis.h"

int main(void) {
    printf("# Frequency Response Demo: 5th-order Butterworth Lowpass\n");
    printf("# fc = 1 kHz\n");
    printf("# Freq_Hz, Mag_dB, Phase_deg, GroupDelay_ms\n");

    filter_spec_t spec = filter_spec_init();
    spec.type = FILTER_LOWPASS;
    spec.approx = APPROX_BUTTERWORTH;
    spec.fc1 = 1000.0;
    spec.order = 5;

    tf_continuous_t *tf = analog_filter_design(&spec);
    if (!tf) { fprintf(stderr, "Design failed\n"); return 1; }

    freq_resp_t *resp = tf_continuous_freq_response(tf, 10.0, 10000.0, 200);
    if (!resp) { tf_continuous_free(tf); return 1; }

    int i;
    for (i = 0; i < resp->num_points; i++) {
        double f = resp->points[i].freq_hz;
        double mag = resp->points[i].magnitude_db;
        double phase = resp->points[i].phase_deg;
        printf("%.2f, %.3f, %.2f, %.3f\n", f, mag, phase,
               resp->points[i].group_delay_s * 1000.0);
    }

    freq_resp_free(resp);
    tf_continuous_free(tf);
    return 0;
}
