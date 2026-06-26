/**
 * ex4_bessel_pulse.c — Bessel Filter: Pulse Preservation
 *
 * L6 Canonical Problem: Compare Butterworth vs Bessel filter for
 * pulse/waveform fidelity.
 *
 * The Bessel filter has maximally flat group delay — it preserves
 * the shape of pulses and transient waveforms with minimal overshoot
 * and ringing. This is critical for:
 *   - Oscilloscope front-ends (Tektronix, Keysight)
 *   - Digital communications (pulse shaping)
 *   - ECG/EKG medical instrumentation (waveform fidelity)
 *   - Audio DAC reconstruction filters
 *
 * This example compares:
 *   - 4th-order Butterworth: overshoot ~11%, significant ringing
 *   - 4th-order Bessel: overshoot < 1%, minimal ringing
 *
 * L7 Application: Medical instrumentation (ECG waveform preservation)
 * Reference: Thomson, "Delay Networks Having Maximally Flat Frequency
 *            Characteristics", Proc. IEE, 1949.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "filter_defs.h"
#include "filter_transfer.h"
#include "filter_analog.h"
#include "filter_digital.h"
#include "filter_design.h"
#include "filter_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Example 4: Bessel Filter — Pulse Preservation ===\n\n");

    double fs = 1000.0;     /* 1 kHz sample rate */
    double fc = 100.0;      /* 100 Hz cutoff */
    int order = 4;

    printf("Comparing %dth-order filters at fc=%.0f Hz:\n\n", order, fc);

    /* Design Butterworth IIR */
    filter_spec_t bw_spec = filter_spec_init();
    bw_spec.type = FILTER_LOWPASS;
    bw_spec.approx = APPROX_BUTTERWORTH;
    bw_spec.fc1 = fc;
    bw_spec.order = order;
    bw_spec.sample_rate = fs;

    iir_filter_t *bw_iir = iir_design_bilinear(&bw_spec);

    /* Design Bessel IIR */
    filter_spec_t bs_spec = filter_spec_init();
    bs_spec.type = FILTER_LOWPASS;
    bs_spec.approx = APPROX_BESSEL;
    bs_spec.fc1 = fc;
    bs_spec.order = order;
    bs_spec.sample_rate = fs;

    iir_filter_t *bs_iir = iir_design_bilinear(&bs_spec);

    if (!bw_iir || !bs_iir) {
        printf("ERROR: Filter design failed!\n");
        if (bw_iir) iir_filter_free(bw_iir);
        if (bs_iir) iir_filter_free(bs_iir);
        return 1;
    }

    /* Frequency response at key points */
    printf("Frequency Response:\n");
    printf("  %-12s %-16s %-16s %-16s %-16s\n",
           "Freq(Hz)", "BW|H|(dB)", "BS|H|(dB)", "BW GD(ms)", "BS GD(ms)");
    double f_tests[] = {10, 50, 100, 150, 200, 300, 500};
    int nf = 7;
    int i;
    for (i = 0; i < nf; i++) {
        double w = 2.0 * M_PI * f_tests[i] / fs;
        freq_resp_point_t fp_bw = iir_freq_response_at(bw_iir, w);
        freq_resp_point_t fp_bs = iir_freq_response_at(bs_iir, w);
        double gd_bw = group_delay_at_iir(bw_iir, w, 0.001);
        double gd_bs = group_delay_at_iir(bs_iir, w, 0.001);
        printf("  %-12.0f %-16.2f %-16.2f %-16.3f %-16.3f\n",
               f_tests[i], fp_bw.magnitude_db, fp_bs.magnitude_db,
               gd_bw * 1000.0 / fs, gd_bs * 1000.0 / fs);
    }

    /* Impulse response comparison */
    int N_out = 128;
    double *h_bw = (double *)malloc(N_out * sizeof(double));
    double *h_bs = (double *)malloc(N_out * sizeof(double));
    if (h_bw && h_bs) {
        iir_impulse_response(bw_iir, N_out, h_bw);
        iir_impulse_response(bs_iir, N_out, h_bs);

        /* Find peak and ringing */
        double max_bw = 0.0, max_bs = 0.0;
        int i;
        for (i = 0; i < N_out; i++) {
            if (fabs(h_bw[i]) > max_bw) max_bw = fabs(h_bw[i]);
            if (fabs(h_bs[i]) > max_bs) max_bs = fabs(h_bs[i]);
        }

        /* Measure ringing: ratio of late samples to peak */
        double ring_bw = 0.0, ring_bs = 0.0;
        for (i = N_out / 4; i < N_out; i++) {
            ring_bw += h_bw[i] * h_bw[i];
            ring_bs += h_bs[i] * h_bs[i];
        }
        ring_bw = sqrt(ring_bw) / max_bw * 100.0;
        ring_bs = sqrt(ring_bs) / max_bs * 100.0;

        printf("\nImpulse Response Analysis:\n");
        printf("  Butterworth: peak=%.4f, late ringing=%.1f%%\n", max_bw, ring_bw);
        printf("  Bessel:      peak=%.4f, late ringing=%.1f%%\n", max_bs, ring_bs);
        printf("  Bessel has %.1fx less ringing — better waveform fidelity\n",
               ring_bw / (ring_bs + 1e-10));

        free(h_bw);
        free(h_bs);
    }

    printf("\nKey insight: Bessel filters are preferred when waveform\n");
    printf("  shape matters more than sharp frequency cutoff.\n");

    iir_filter_free(bw_iir);
    iir_filter_free(bs_iir);
    printf("\n=== Example 4 complete ===\n");
    return 0;
}
