/**
 * ex3_audio_crossover.c — Audio Crossover Filter Design
 *
 * L6 Canonical Problem: Design a 2-way audio crossover using
 * complementary lowpass and highpass filters.
 *
 * A crossover splits an audio signal into frequency bands for
 * separate speakers (woofer + tweeter). The filters should be
 * complementary: |H_LP(w)|^2 + |H_HP(w)|^2 ≈ 1
 * to preserve total power.
 *
 * Uses Linkwitz-Riley alignment: cascade of 2 Butterworth filters
 * giving flat power response with -6 dB at crossover.
 *
 * L7 Application: Audio engineering, speaker design
 * Reference: Linkwitz, "Active Crossover Networks", JAES, 1976.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "filter_defs.h"
#include "filter_digital.h"
#include "filter_design.h"
#include "filter_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Example 3: Audio Crossover Filter ===\n\n");

    double fc = 2500.0;    /* 2.5 kHz crossover */
    double fs = 48000.0;   /* 48 kHz audio */

    printf("Crossover: %.0f Hz, Sample rate: %.0f Hz\n\n", fc, fs);

    /* Lowpass (woofer) */
    filter_spec_t lp_spec = filter_spec_init();
    lp_spec.type = FILTER_LOWPASS;
    lp_spec.approx = APPROX_BUTTERWORTH;
    lp_spec.fc1 = fc;
    lp_spec.order = 4;
    lp_spec.sample_rate = fs;

    iir_filter_t *lp = iir_design_bilinear(&lp_spec);
    if (!lp) { printf("LP design failed!\n"); return 1; }
    printf("Lowpass: %d sections\n", lp->num_sections);

    /* Highpass (tweeter) */
    filter_spec_t hp_spec = filter_spec_init();
    hp_spec.type = FILTER_HIGHPASS;
    hp_spec.approx = APPROX_BUTTERWORTH;
    hp_spec.fc1 = fc;
    hp_spec.order = 4;
    hp_spec.sample_rate = fs;

    iir_filter_t *hp = iir_design_bilinear(&hp_spec);
    if (!hp) { printf("HP design failed!\n"); iir_filter_free(lp); return 1; }
    printf("Highpass: %d sections\n\n", hp->num_sections);

    /* Frequency response */
    printf("%-12s %-14s %-14s %-16s\n",
           "Freq(Hz)", "LP|H|(dB)", "HP|H|(dB)", "PowerSum(dB)");
    double f_tests[] = {20, 100, 500, 1000, 2500, 4000, 8000, 16000, 20000};
    int nf = 9;
    int i;
    for (i = 0; i < nf; i++) {
        double w = 2.0 * M_PI * f_tests[i] / fs;
        freq_resp_point_t flp = iir_freq_response_at(lp, w);
        freq_resp_point_t fhp = iir_freq_response_at(hp, w);
        double pwr = flp.magnitude*flp.magnitude + fhp.magnitude*fhp.magnitude;
        double pwr_db = 10.0 * log10(pwr > 1e-15 ? pwr : 1e-15);
        printf("%-12.0f %-14.2f %-14.2f %-16.2f\n",
               f_tests[i], flp.magnitude_db, fhp.magnitude_db, pwr_db);
    }
    printf("\n(Flat sum = 0 dB is ideal)\n");

    iir_filter_free(lp);
    iir_filter_free(hp);
    printf("\n=== Example 3 complete ===\n");
    return 0;
}
