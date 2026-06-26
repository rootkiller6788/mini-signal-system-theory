/**
 * @file    example_vibration_analysis.c
 * @brief   Example: Vibration Analysis for Bearing Fault Detection (L7 Application)
 *
 * Demonstrates the use of Fourier analysis for rotating machinery
 * condition monitoring, a core L7 application studied at Michigan
 * (automotive) and TU Munich (industrial engineering).
 *
 * The pipeline:
 *   1. Simulate vibration signal from a healthy bearing
 *   2. Add fault signature (inner race defect at BPFI frequency)
 *   3. Compute vibration spectrum via Welch's method
 *   4. Detect characteristic bearing fault frequencies
 *   5. Compute Total Harmonic Distortion (for electrical motor health)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fourier_core.h"
#include "fourier_window.h"
#include "fourier_spectrum.h"
#include "fourier_applications.h"

int main(void) {
    printf("=== Example: Vibration Analysis for Bearing Diagnostics ===\n\n");

    /* ── Bearing parameters (SKF 6205 deep groove ball bearing) ── */
    double shaft_freq   = 30.0;     /* 1800 RPM */
    int    n_balls      = 9;
    double ball_diam    = 7.94e-3;  /* 7.94 mm */
    double pitch_diam   = 39.04e-3; /* 39.04 mm */
    double contact_angle = 0.0;     /* radians */

    /* Characteristic frequencies */
    double ratio = ball_diam / pitch_diam;
    double bpfo = ((double)n_balls / 2.0) * shaft_freq * (1.0 - ratio);
    double bpfi = ((double)n_balls / 2.0) * shaft_freq * (1.0 + ratio);
    double bsf  = (pitch_diam / (2.0 * ball_diam)) * shaft_freq * (1.0 - ratio * ratio);
    double ftf  = (shaft_freq / 2.0) * (1.0 - ratio);

    printf("Bearing: SKF 6205\n");
    printf("  Shaft speed:       %.0f RPM (%.2f Hz)\n", shaft_freq * 60.0, shaft_freq);
    printf("  BPFO (outer race): %.2f Hz\n", bpfo);
    printf("  BPFI (inner race): %.2f Hz\n", bpfi);
    printf("  BSF (ball spin):   %.2f Hz\n", bsf);
    printf("  FTF (cage):        %.2f Hz\n\n", ftf);

    /* ── Simulate vibration signal ── */
    double fs = 10000.0;    /* 10 kHz sampling */
    double duration = 1.0;  /* 1 second */
    size_t N = (size_t)(fs * duration);

    double *vibration = (double*)malloc(N * sizeof(double));
    if (!vibration) { fprintf(stderr, "malloc failed\n"); return 1; }

    for (size_t i = 0; i < N; i++) {
        double t = (double)i / fs;
        /* Fundamental rotation */
        vibration[i] = 1.0 * sin(2.0 * M_PI * shaft_freq * t);
        /* Shaft harmonics (unbalance, misalignment) */
        vibration[i] += 0.5 * sin(2.0 * M_PI * 2.0 * shaft_freq * t);
        vibration[i] += 0.3 * sin(2.0 * M_PI * 3.0 * shaft_freq * t);
        /* Inner race fault (key feature) */
        vibration[i] += 0.4 * sin(2.0 * M_PI * bpfi * t);
        vibration[i] += 0.2 * sin(2.0 * M_PI * 2.0 * bpfi * t);
        /* Gaussian noise (sensor + background) */
        vibration[i] += 0.15 * ((double)rand() / (double)RAND_MAX - 0.5);
    }

    printf("Vibration signal: %zu samples at %.0f Hz\n", N, fs);

    /* ── Spectral analysis ── */
    size_t seg_len = 4096;
    psd_result_t psd = welch_psd(vibration, N, fs, seg_len, 0.5, WINDOW_HANN);

    printf("PSD: %zu segments of %zu points, resolution %.2f Hz\n\n",
           psd.num_segments, psd.N, psd.freq_resolution_hz);

    /* ── Extract magnitude spectrum for bearing fault detection ── */
    double *magnitude = (double*)malloc(psd.N * sizeof(double));
    double *freq_axis = (double*)malloc(psd.N * sizeof(double));
    for (size_t k = 0; k < psd.N; k++) {
        magnitude[k] = sqrt(psd.psd[k]);
        freq_axis[k] = (double)k * fs / (double)psd.N;
    }

    /* ── Detect bearing faults ── */
    double bpfo_peak, bpfi_peak, severity;
    bearing_fault_detect(magnitude, freq_axis, psd.N,
                          shaft_freq, n_balls, ball_diam, pitch_diam,
                          contact_angle,
                          &bpfo_peak, &bpfi_peak, &severity);

    printf("Fault detection results:\n");
    printf("  BPFO peak:   %.4f\n", bpfo_peak);
    printf("  BPFI peak:   %.4f\n", bpfi_peak);
    printf("  Severity:    %.2f (%.1f = critical)\n\n", severity,
           (severity > 6.0) ? 6.0 : severity > 0 ? severity : 0.0);

    if (bpfi_peak > bpfo_peak * 2.0) {
        printf("  → DIAGNOSIS: Inner race fault detected (BPFI dominant)\n\n");
    } else if (bpfo_peak > bpfi_peak * 2.0) {
        printf("  → DIAGNOSIS: Outer race fault detected (BPFO dominant)\n\n");
    } else {
        printf("  → DIAGNOSIS: No clear bearing fault signature\n\n");
    }

    /* ── THD analysis for motor electrical health ── */
    double thd = total_harmonic_distortion(vibration, N, fs, shaft_freq, 15);
    printf("Electrical health (THD at shaft frequency):\n");
    printf("  THD = %.2f%%\n", thd * 100.0);
    printf("  IEEE 519 limit: < 5%% for low-voltage systems\n");
    if (thd * 100.0 < 5.0) {
        printf("  → PASS: THD within IEEE 519 limits\n");
    } else {
        printf("  → WARNING: THD exceeds IEEE 519 limits\n");
    }

    /* ── Harmonic analysis ── */
    double *harmonics = (double*)malloc(20 * sizeof(double));
    double thd_percent;
    power_harmonics(vibration, N, fs, shaft_freq, 20, harmonics, &thd_percent);

    printf("\nHarmonic amplitudes (normalized to fundamental):\n");
    for (size_t h = 0; h < 10; h++) {
        if (harmonics[h] > 0.01) {
            printf("  H%-2zu: %.4f\n", h + 1, harmonics[h] / harmonics[0]);
        }
    }

    /* ── Bandwidth analysis of the vibration spectrum ── */
    bandwidth_metrics_t bw = spectrum_bandwidth(psd.psd, psd.N, fs, 0.99);
    printf("\nVibration bandwidth metrics:\n");
    printf("  3 dB bandwidth:     %.1f Hz\n", bw.bandwidth_3dB_hz);
    printf("  Equivalent noise BW: %.1f Hz\n", bw.bandwidth_enb_hz);
    printf("  99%% power bandwidth: %.1f Hz\n", bw.bandwidth_occupied_hz);
    printf("  Center frequency:   %.1f Hz\n", bw.center_freq_hz);

    /* ── Cleanup ── */
    free(vibration);
    free(magnitude);
    free(freq_axis);
    free(harmonics);
    psd_result_free(&psd);

    printf("\n=== Example Complete ===\n");
    return 0;
}
