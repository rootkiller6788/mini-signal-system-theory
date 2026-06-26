/**
 * bench_filter.c — Performance Benchmarks
 *
 * Measures throughput of FIR and IIR filter implementations
 * in samples per second. Useful for selecting the right
 * implementation for real-time applications.
 *
 * L6: Performance characterization of filter implementations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "filter_defs.h"
#include "filter_digital.h"
#include "filter_design.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== Filter Performance Benchmarks ===\n\n");

    int trials = 5;
    int block_len = 100000;
    double *x = (double *)malloc(block_len * sizeof(double));
    double *y = (double *)malloc(block_len * sizeof(double));
    if (!x || !y) { free(x); free(y); return 1; }
    int i;
    for (i = 0; i < block_len; i++) x[i] = sin(2.0 * M_PI * 100.0 * i / 48000.0);

    printf("%-20s %-10s %-12s %-12s\n", "Filter", "Taps", "Time(ms)", "MSamples/s");
    printf("----------------------------------------------------------\n");

    /* FIR benchmarks at various lengths */
    int taps[] = {16, 32, 64, 128, 256};
    int nt = 5;
    int t;
    for (t = 0; t < nt; t++) {
        fir_filter_t *fir = fir_design_window(0.125, taps[t], WINDOW_HAMMING, 0.0);
        if (!fir) continue;
        filter_state_t *st = fir_state_alloc(taps[t]);
        if (!st) { fir_filter_free(fir); continue; }

        double t_start = get_time_sec();
        int tr;
        for (tr = 0; tr < trials; tr++) {
            fir_process_block(fir, st, x, y, block_len);
            filter_state_reset(st);
        }
        double t_elapsed = (get_time_sec() - t_start) / trials;
        double msps = (block_len / 1e6) / t_elapsed;

        printf("FIR(%-12s) %-10d %-12.3f %-12.2f\n",
               "lowpass", taps[t], t_elapsed * 1000.0, msps);

        filter_state_free(st);
        fir_filter_free(fir);
    }

    /* IIR benchmark */
    filter_spec_t spec = filter_spec_init();
    spec.sample_rate = 48000.0;
    spec.fc1 = 1000.0;
    spec.order = 4;

    iir_filter_t *iir = iir_design_bilinear(&spec);
    if (iir) {
        filter_state_t *st = (filter_state_t *)calloc(1, sizeof(filter_state_t));
        if (st) {
            double t_start = get_time_sec();
            int tr;
            for (tr = 0; tr < trials; tr++) {
                iir_process_block(iir, st, x, y, block_len);
                filter_state_reset(st);
            }
            double t_elapsed = (get_time_sec() - t_start) / trials;
            double msps = (block_len / 1e6) / t_elapsed;

            printf("IIR(bilinear,4th)   %-10d %-12.3f %-12.2f\n",
                   iir->num_sections * 2, t_elapsed * 1000.0, msps);
            filter_state_free(st);
        }
        iir_filter_free(iir);
    }

    free(x);
    free(y);
    printf("\n=== Benchmark complete ===\n");
    return 0;
}
