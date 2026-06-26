/**
 * ex2_fir_design.c — FIR Filter Design by Windowing
 *
 * L6 Canonical Problem: Design an FIR lowpass filter using the window method,
 * compare different windows (Hann, Hamming, Blackman, Kaiser),
 * and analyze the trade-offs.
 *
 * Key trade-offs demonstrated:
 *   - Hann:  good side-lobe roll-off (18 dB/octave)
 *   - Hamming: optimized first side lobe (-42.7 dB)
 *   - Blackman: deep side lobes (-58 dB)
 *   - Kaiser:  adjustable via beta parameter
 *
 * Reference: Harris, "On the Use of Windows for Harmonic Analysis", 1978
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
    printf("=== Example 2: FIR Filter Design by Windowing ===\n\n");

    double fs = 10000.0;         /* 10 kHz sample rate */
    double fc = 1000.0;          /* 1 kHz cutoff */
    int num_taps = 51;           /* 51-tap FIR */
    double wc_norm = fc / (fs / 2.0); /* Normalized cutoff (0 to 1) */

    printf("Parameters:\n");
    printf("  Sample rate: %.0f Hz\n", fs);
    printf("  Cutoff freq: %.0f Hz\n", fc);
    printf("  Normalized cutoff: %.3f (fraction of Nyquist)\n", wc_norm);
    printf("  Number of taps: %d\n\n", num_taps);

    /* Design with different windows */
    window_type_t windows[] = {
        WINDOW_HANN, WINDOW_HAMMING, WINDOW_BLACKMAN, WINDOW_KAISER
    };
    const char *win_names[] = {"Hann", "Hamming", "Blackman", "Kaiser(beta=5)"};
    double win_alphas[] = {0.0, 0.0, 0.0, 5.0};
    int nw = 4;
    int w;

    for (w = 0; w < nw; w++) {
        /* Design filter */
        fir_filter_t *fir = fir_design_window(wc_norm / 2.0, num_taps,
                                               windows[w], win_alphas[w]);
        if (!fir) {
            printf("  %s: DESIGN FAILED\n", win_names[w]);
            continue;
        }

        /* Analyze frequency response */
        freq_resp_t *resp = fir_freq_response(fir, 0.0, fs / 2.0, 256, fs);
        if (!resp) {
            printf("  %s: RESPONSE FAILED\n", win_names[w]);
            fir_filter_free(fir);
            continue;
        }

        /* Find key metrics */
        double passband_flat = passband_flatness(resp, 0.0, fc);
        double stopband_att = stopband_attenuation(resp, fc * 1.5, fs / 2.0);
        double f3db = find_3db_cutoff(resp, resp->points[0].magnitude_db);
        double gd = fir->group_delay;

        printf("  %-20s:\n", win_names[w]);
        printf("    DC gain:      %.3f\n", fir->gain_dc);
        printf("    -3dB cutoff:  %.1f Hz\n", f3db);
        printf("    Passband flatness: %.3f dB\n", passband_flat);
        printf("    Stopband atten:    %.1f dB\n", stopband_att);
        printf("    Group delay:  %.1f samples (%.3f ms)\n",
               gd, gd * 1000.0 / fs);
        printf("    Linear phase: Type %d\n\n", fir->lp_type);

        freq_resp_free(resp);
        fir_filter_free(fir);
    }

    /* Order estimation */
    double trans_bw = 2.0 * M_PI * 500.0 / fs; /* 500 Hz transition */
    double atten_db = 60.0;
    int est_order = fir_estimate_length(trans_bw, atten_db);
    printf("Estimated FIR length for %.0f dB attenuation: %d taps\n",
           atten_db, est_order);

    double beta_opt = kaiser_beta_from_attenuation(atten_db);
    printf("Optimal Kaiser beta for %.0f dB: %.3f\n", atten_db, beta_opt);

    printf("\n=== Example 2 complete ===\n");
    return 0;
}
