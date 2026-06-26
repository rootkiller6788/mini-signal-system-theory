/**
 * @file bench_convolution.c
 * @brief Performance Benchmarks: Conv. Algorithms and Matrix Operations
 *
 * Benchmarks key signal processing primitives:
 *   - DT direct convolution O(N^2) vs FFT-based O(N log N)
 *   - CT convolution (trapezoidal integration)
 *   - Overlap-add / Overlap-save block convolution
 *   - Matrix exponential (scaling+squaring)
 *   - Pole-zero finding (Newton+deflation)
 *   - DFT/IDFT direct computation
 *
 * Each benchmark measures wall-clock time and reports
 * operations per second where applicable.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/system_defs.h"
#include "../include/convolution.h"
#include "../include/transfer_function.h"
#include "../include/state_space.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double get_time_sec(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

#define BENCH(name, N, block) do { \
    double t0 = get_time_sec(); \
    for (int _r = 0; _r < (N); _r++) { block; } \
    double t1 = get_time_sec(); \
    double elapsed = (t1 - t0) / (double)(N); \
    printf("  %-40s %10.6f sec/op", name, elapsed); \
    if (elapsed > 0.0) printf("  (%10.0f ops/sec)", 1.0/elapsed); \
    printf("\n"); \
} while(0)

int main(void) {
    printf("=== Convolution and DSP Performance Benchmarks ===\n\n");
    srand(42);

    /* Prepare test signals */
    int Nx=2000, Nh=200;
    dt_signal_t x = dt_signal_alloc(Nx);
    dt_signal_t h = dt_signal_alloc(Nh);
    dt_signal_t y = dt_signal_alloc(Nx+Nh-1);
    for (int i=0;i<Nx;i++) x.samples[i]=(double)rand()/(double)RAND_MAX - 0.5;
    for (int i=0;i<Nh;i++) h.samples[i]=exp(-(double)i*0.05);

    printf("--- DT Linear Convolution (len(x)=%d, len(h)=%d) ---\n", Nx, Nh);
    BENCH("Direct O(N^2)", 5, dt_convolution_direct(&x, &h, &y));
    BENCH("FFT-based O(N log N)", 5, dt_convolution_fft(&x, &h, &y));
    BENCH("Overlap-Add (block=512)", 5, dt_overlap_add(&x, &h, 512, &y));
    BENCH("Overlap-Save (block=512)", 5, dt_overlap_save(&x, &h, 512, &y));

    dt_signal_free(&x); dt_signal_free(&h); dt_signal_free(&y);

    /* Circular convolution */
    int Nc=1024;
    dt_signal_t xc = dt_signal_alloc(Nc);
    dt_signal_t hc = dt_signal_alloc(Nc);
    dt_signal_t yc = dt_signal_alloc(Nc);
    for (int i=0;i<Nc;i++) { xc.samples[i]=sin(2.0*M_PI*(double)i/64.0); hc.samples[i]=exp(-(double)i/100.0); }
    printf("\n--- DT Circular Convolution (N=%d) ---\n", Nc);
    BENCH("Circular Direct O(N^2)", 10, dt_circular_convolution_direct(&xc, &hc, &yc));
    BENCH("Circular FFT-based", 10, dt_circular_convolution_fft(&xc, &hc, &yc));
    dt_signal_free(&xc); dt_signal_free(&hc); dt_signal_free(&yc);

    /* Correlation */
    dt_signal_t xcorr = dt_signal_alloc(500);
    dt_signal_t ycorr = dt_signal_alloc(500);
    dt_signal_t rxy = dt_signal_alloc(999);
    for (int i=0;i<500;i++) { xcorr.samples[i]=(double)rand()/(double)RAND_MAX; ycorr.samples[i]=xcorr.samples[i]+0.1*((double)rand()/(double)RAND_MAX-0.5); }
    printf("\n--- Correlation (N=500) ---\n");
    BENCH("Cross-correlation", 20, dt_cross_correlation(&xcorr, &ycorr, &rxy));
    BENCH("Auto-correlation", 20, dt_auto_correlation(&xcorr, &rxy));
    dt_signal_free(&xcorr); dt_signal_free(&ycorr); dt_signal_free(&rxy);

    /* Matrix exponential */
    ct_state_space_t ss = ct_ss_alloc(4, 1, 1);
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) ss.A[i*4+j]=(i==j)?-1.0:((i+1==j)?1.0:0.0);
    double *expAt = (double*)malloc(16*sizeof(double));
    printf("\n--- Matrix Exponential (4x4, scaling+squaring) ---\n");
    BENCH("expm(At) for t=1.0", 100, ct_ss_matrix_exponential(&ss, 1.0, expAt));
    BENCH("expm(At) for t=5.0", 100, ct_ss_matrix_exponential(&ss, 5.0, expAt));
    free(expAt); ct_ss_free(&ss);

    /* Pole-zero finding */
    ct_transfer_function_t tf = ct_tf_alloc(0, 8);
    tf.den_coeffs[8]=1.0; tf.den_coeffs[7]=5.0; tf.den_coeffs[6]=15.0; tf.den_coeffs[5]=35.0;
    tf.den_coeffs[4]=70.0; tf.den_coeffs[3]=126.0; tf.den_coeffs[2]=210.0; tf.den_coeffs[1]=330.0; tf.den_coeffs[0]=429.0;
    pole_zero_collection_t pz = pz_collection_alloc(8, 0);
    printf("\n--- Pole Finding (8th-order polynomial) ---\n");
    BENCH("Newton+deflation (8 roots)", 20, ct_tf_find_poles(&tf, &pz));
    pz_collection_free(&pz); ct_tf_free(&tf);

    /* Eigenvalue computation */
    ct_state_space_t ss2 = ct_ss_alloc(5, 1, 1);
    for (int i=0;i<5;i++) for (int j=0;j<5;j++) ss2.A[i*5+j]=(i==j)?(double)(-i-1)*0.5:(i>j?1.0:-0.5);
    double complex evals[5];
    printf("\n--- Eigenvalue Computation (5x5, Faddeev-LeVerrier) ---\n");
    BENCH("Char poly + Newton (5 evals)", 50, ct_ss_eigenvalues(&ss2, evals));
    ct_ss_free(&ss2);

    /* Frequency response sweep */
    ct_transfer_function_t tf_fr = ct_tf_alloc(0, 4);
    tf_fr.num_coeffs[0]=1.0; tf_fr.den_coeffs[0]=1.0; tf_fr.den_coeffs[1]=2.613; tf_fr.den_coeffs[2]=3.414; tf_fr.den_coeffs[3]=2.613; tf_fr.den_coeffs[4]=1.0;
    frequency_response_t fr = freq_resp_alloc(10000, 0.01, 1000.0, 1);
    printf("\n--- Frequency Response (10000 points, 4th-order) ---\n");
    BENCH("Bode sweep 10000 pts", 10, ct_frequency_response(&tf_fr, &fr));
    freq_resp_free(&fr); ct_tf_free(&tf_fr);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
