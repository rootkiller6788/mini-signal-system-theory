/**
 * @file    example_spectrum_analyzer.c
 * @brief   Example: Audio Spectrum Analyzer (L6 Canonical Problem)
 *
 * Demonstrates a complete spectral analysis pipeline:
 *   1. Generate a multi-tone test signal (simulating audio)
 *   2. Apply window functions to mitigate spectral leakage
 *   3. Compute Welch's PSD estimate
 *   4. Extract octave-band energy (human auditory model)
 *   5. Verify Parseval's theorem
 *
 * This is a canonical L6 problem: frequency-domain characterization
 * of a measured signal using the DFT.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fourier_core.h"
#include "fourier_window.h"
#include "fourier_spectrum.h"
#include "fourier_convolution.h"
#include "fourier_applications.h"

int main(void) {
    printf("=== Example: Audio Spectrum Analyzer ===\n\n");

    /* ── Signal generation: 3-tone audio signal ── */
    double fs = 44100.0;         /* CD-quality sampling rate */
    double duration = 0.5;       /* 500 ms */
    size_t N = (size_t)(fs * duration);  /* 22050 samples */

    double *signal = (double*)malloc(N * sizeof(double));
    if (!signal) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* Tones: A4 (440 Hz), C#5 (554.37 Hz), E5 (659.25 Hz) — A major chord */
    double f1 = 440.0, f2 = 554.37, f3 = 659.25;
    double a1 = 1.0, a2 = 0.7, a3 = 0.5;  /* amplitudes */

    for (size_t i = 0; i < N; i++) {
        double t = (double)i / fs;
        signal[i] = a1 * sin(2.0 * M_PI * f1 * t)
                  + a2 * sin(2.0 * M_PI * f2 * t)
                  + a3 * sin(2.0 * M_PI * f3 * t);
        /* Add low-level white noise (simulating measurement noise) */
        signal[i] += 0.01 * ((double)rand() / (double)RAND_MAX - 0.5);
    }

    printf("Signal: %zu samples at %.0f Hz (%.2f s)\n", N, fs, duration);
    printf("Tones: %.1f Hz, %.2f Hz, %.2f Hz\n", f1, f2, f3);
    printf("Expected peaks at bins: %.0f, %.0f, %.0f\n\n",
           f1 * N / fs, f2 * N / fs, f3 * N / fs);

    /* ── Verify Parseval: check energy conservation ── */
    double parseval_err = verify_parseval_theorem(signal, N);
    printf("Parseval theorem error: %.2e (should be ~0)\n\n", parseval_err);

    /* ── Spectral Estimation: Welch's method ── */
    size_t seg_len = 4096;
    double overlap = 0.5;
    psd_result_t psd = welch_psd(signal, N, fs, seg_len, overlap, WINDOW_HANN);

    /* Print top 10 spectral peaks */
    printf("Welch PSD (%zu segments, %zu-point FFT, Hann window):\n",
           psd.num_segments, psd.N);
    printf("  Freq resolution: %.2f Hz\n", psd.freq_resolution_hz);
    printf("  Variance reduction factor: %.1f\n\n", psd.variance_reduction);

    /* Find and display peaks */
    printf("Top 10 frequency peaks:\n");
    printf("  %-12s %-15s %-15s\n", "Freq (Hz)", "PSD (V²/Hz)", "dB rel max");

    /* Simple peak picking: local maxima above threshold */
    double max_psd = 0.0;
    for (size_t k = 0; k < psd.N; k++) {
        if (psd.psd[k] > max_psd) max_psd = psd.psd[k];
    }

    double threshold = max_psd * 0.01;
    size_t peaks_found = 0;
    for (size_t k = 2; k < psd.N - 1 && peaks_found < 10; k++) {
        if (psd.psd[k] > threshold
            && psd.psd[k] > psd.psd[k-1]
            && psd.psd[k] > psd.psd[k+1]) {
            double freq = (double)k * fs / (double)psd.N;
            double db = 10.0 * log10(psd.psd[k] / max_psd);
            printf("  %-12.2f %-15.6e %-15.2f\n", freq, psd.psd[k], db);
            peaks_found++;
        }
    }
    printf("\n");

    /* ── Octave-band analysis ── */
    double *band_energy = NULL, *band_centers = NULL;
    int num_bands = 0;
    int ret = audio_octave_spectrum(signal, N, fs, 1,
                                     &band_energy, &band_centers, &num_bands);
    if (ret == 0 && num_bands > 0) {
        printf("Octave-band spectrum (%d bands):\n", num_bands);
        printf("  %-12s %-15s\n", "Center (Hz)", "Energy");
        for (int b = 0; b < num_bands; b++) {
            printf("  %-12.1f %-15.6e\n", band_centers[b], band_energy[b]);
        }
        free(band_energy);
        free(band_centers);
    }

    /* ── Compare window functions on the same signal ── */
    printf("\nWindow comparison (freq bin 20, rel to max):\n");
    printf("  %-16s %-12s %-12s\n", "Window", "Mainlobe(bins)", "Sidelobe(dB)");

    window_type_t wins[] = {WINDOW_RECTANGULAR, WINDOW_HANN, WINDOW_HAMMING,
                            WINDOW_BLACKMAN, WINDOW_BARTLETT};
    const char *wnames[] = {"Rectangular", "Hann", "Hamming", "Blackman", "Bartlett"};

    for (int w = 0; w < 5; w++) {
        window_config_t wc;
        switch (wins[w]) {
            case WINDOW_RECTANGULAR: wc = window_rectangular(seg_len); break;
            case WINDOW_HANN:        wc = window_hann(seg_len); break;
            case WINDOW_HAMMING:     wc = window_hamming(seg_len); break;
            case WINDOW_BLACKMAN:    wc = window_blackman(seg_len); break;
            case WINDOW_BARTLETT:    wc = window_bartlett(seg_len); break;
            default: continue;
        }
        printf("  %-16s %-12.2f %-12.1f\n",
               wnames[w], wc.mainlobe_3dB_width, wc.sidelobe_level_db);
        window_config_free(&wc);
    }

    /* ── Cleanup ── */
    psd_result_free(&psd);
    free(signal);

    printf("\n=== Example Complete ===\n");
    return 0;
}
