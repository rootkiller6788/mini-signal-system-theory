/**
 * @file    example_radar_pulse.c
 * @brief   Example: Radar LFM Pulse Compression (L6 Canonical Problem)
 *
 * Demonstrates linear frequency modulated (LFM / chirp) pulse generation
 * and matched-filter pulse compression — a canonical problem in radar
 * signal processing.
 *
 * Pulse compression achieves high range resolution without requiring
 * high peak power by transmitting a long coded pulse and compressing
 * the received signal via matched filtering.
 *
 * Compression ratio = time-bandwidth product (B·T).
 * Range resolution = c / (2·B)  (c = speed of light).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fourier_core.h"
#include "fourier_convolution.h"
#include "fourier_applications.h"

int main(void) {
    printf("=== Example: Radar LFM Pulse Compression ===\n\n");

    /* ── Radar parameters ── */
    double pulse_width  = 10e-6;      /* 10 µs pulse */
    double bandwidth    = 50e6;       /* 50 MHz bandwidth */
    double fs           = 100e6;      /* 100 MHz sampling (Nyquist for 50 MHz BW) */
    double c            = 3e8;        /* Speed of light (m/s) */

    double range_res = c / (2.0 * bandwidth);
    printf("Pulse width:       %.1f µs\n", pulse_width * 1e6);
    printf("Bandwidth:         %.0f MHz\n", bandwidth / 1e6);
    printf("Time-BW product:   %.0f\n", pulse_width * bandwidth);
    printf("Range resolution:  %.2f m\n\n", range_res);

    /* ── Generate LFM chirp ── */
    size_t N_tx = 0;
    double *tx_real = NULL, *tx_imag = NULL;
    chirp_generate_lfm(pulse_width, bandwidth, fs, &N_tx, &tx_real, &tx_imag);

    fourier_complex_t *tx = (fourier_complex_t*)malloc(N_tx * sizeof(fourier_complex_t));
    for (size_t i = 0; i < N_tx; i++) {
        tx[i] = tx_real[i] + I * tx_imag[i];
    }
    printf("Transmit pulse: %zu samples\n", N_tx);

    /* ── Simulate radar returns: 3 targets at different ranges ── */
    /* Target ranges: 1500 m, 3000 m, 4500 m → delays: 10 µs, 20 µs, 30 µs */
    double ranges[] = {1500.0, 3000.0, 4500.0};
    double delays[] = {10e-6, 20e-6, 30e-6};  /* round-trip delays */
    double rcs[]    = {1.0, 0.7, 0.3};         /* radar cross-section (relative) */

    /* Build received signal: superposition of delayed, attenuated chirps */
    size_t rx_len = (size_t)((delays[2] + pulse_width) * fs) + 100;
    double *received = (double*)calloc(rx_len, sizeof(double));

    for (int t = 0; t < 3; t++) {
        size_t delay_samples = (size_t)(delays[t] * fs);
        for (size_t i = 0; i < N_tx && (delay_samples + i) < rx_len; i++) {
            received[delay_samples + i] += rcs[t] * tx_real[i];
        }
    }
    printf("Received signal: %zu samples (3 targets)\n", rx_len);

    /* ── Pulse compression via matched filter ── */
    double *compressed = (double*)malloc((rx_len + N_tx - 1) * sizeof(double));
    pulse_compression_matched_filter(received, rx_len, tx, N_tx, compressed);

    /* ── Find peaks in the compressed output ── */
    size_t comp_len = rx_len + N_tx - 1;
    printf("\nPulse compression result (%zu samples):\n", comp_len);

    /* Peak detection: find local maxima above threshold */
    double max_val = 0.0;
    for (size_t i = 0; i < comp_len; i++) {
        if (compressed[i] > max_val) max_val = compressed[i];
    }
    double threshold = max_val * 0.3;

    printf("  %-10s %-15s %-18s %-15s\n", "Delay(µs)", "Range(m)", "Amplitude", "Expected RCS");
    for (size_t i = 1; i < comp_len - 1; i++) {
        if (compressed[i] > threshold
            && compressed[i] > compressed[i-1]
            && compressed[i] > compressed[i+1]) {
            double delay = (double)i / fs;
            double range = delay * c / 2.0;
            printf("  %-10.2f %-15.1f %-18.4f ", delay * 1e6, range, compressed[i]);
            /* Match to nearest expected range */
            for (int t = 0; t < 3; t++) {
                if (fabs(range - ranges[t]) < 100.0) {
                    printf("%-15.2f", rcs[t]);
                    break;
                }
            }
            printf("\n");
        }
    }

    /* ── Compression ratio verification ── */
    /* The compressed pulse width ≈ 1/B = 20 ns, compared to 10 µs original */
    double compression_ratio = pulse_width * bandwidth;
    printf("\nCompression ratio (T·B): %.0f\n", compression_ratio);
    printf("Uncompressed pulse: %.2f µs → Compressed: %.2f ns\n",
           pulse_width * 1e6, (1.0 / bandwidth) * 1e9);

    /* ── Cleanup ── */
    free(tx_real); free(tx_imag); free(tx);
    free(received); free(compressed);

    printf("\n=== Example Complete ===\n");
    return 0;
}
