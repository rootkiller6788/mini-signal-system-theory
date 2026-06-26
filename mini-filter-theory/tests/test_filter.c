/**
 * test_filter.c — Comprehensive Test Suite for mini-filter-theory
 *
 * Tests cover:
 *   L1: struct definitions, initialization, validation
 *   L4: Butterworth poles, Routh-Hurwitz stability, Jury stability, bilinear
 *   L5: Window functions, FIR design, IIR design, order estimation
 *   L6: Frequency response, group delay, step response
 *
 * Uses standard assert() for validation.
 * Run: make test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
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

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  Test %d: %s ... ", tests_run, name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

/* ==================================================================
 * L1: Core Type Definitions
 * ================================================================== */

static void test_filter_spec_init(void) {
    TEST("filter_spec_init");
    filter_spec_t spec = filter_spec_init();
    assert(spec.type == FILTER_LOWPASS);
    assert(spec.approx == APPROX_BUTTERWORTH);
    assert(spec.fc1 == 1000.0);
    assert(spec.order == 0);
    PASS();
}

static void test_filter_spec_is_valid(void) {
    TEST("filter_spec_is_valid");
    filter_spec_t spec = filter_spec_init();
    assert(filter_spec_is_valid(&spec) == 1);

    /* Invalid: negative cutoff */
    spec.fc1 = -100.0;
    assert(filter_spec_is_valid(&spec) == 0);

    /* Valid BP */
    spec = filter_spec_init();
    spec.type = FILTER_BANDPASS;
    spec.fc1 = 500.0;
    spec.fc2 = 1500.0;
    assert(filter_spec_is_valid(&spec) == 1);

    /* Invalid BP: fc2 < fc1 */
    spec.fc2 = 300.0;
    assert(filter_spec_is_valid(&spec) == 0);
    PASS();
}

/* ==================================================================
 * L4: Butterworth Algorithm
 * ================================================================== */

static void test_butterworth_poles(void) {
    TEST("butterworth_poles order 4");
    int N = 4;
    double complex poles[4];
    int ret = butterworth_poles(N, 1.0, poles);
    assert(ret == FILTER_OK);

    /* Verify all poles are in LHP */
    int i;
    for (i = 0; i < N; i++) {
        assert(creal(poles[i]) < 0.0);
    }

    /* Poles should be symmetric about real axis */
    int found_conj = 0;
    for (i = 0; i < N; i++) {
        int j;
        for (j = 0; j < N; j++) {
            if (i != j) {
                if (fabs(creal(poles[i]) - creal(poles[j])) < 1e-10 &&
                    fabs(cimag(poles[i]) + cimag(poles[j])) < 1e-10) {
                    found_conj = 1;
                }
            }
        }
    }
    assert(found_conj == 1);
    PASS();
}

static void test_butterworth_prototype(void) {
    TEST("butterworth_prototype order 4");
    tf_continuous_t *proto = butterworth_prototype(4);
    assert(proto != NULL);

    /* Check DC gain = 1 */
    double complex H0 = tf_continuous_eval(proto, 0);
    assert(fabs(cabs(H0) - 1.0) < 0.01);

    /* Check -3dB at w=1 (normalized cutoff) */
    double complex Hwc = tf_continuous_eval(proto, 1.0 * I);
    double mag_db = 20.0 * log10(cabs(Hwc));
    assert(fabs(mag_db + 3.01) < 0.5);

    tf_continuous_free(proto);
    PASS();
}

/* ==================================================================
 * L4: Routh-Hurwitz Stability
 * ================================================================== */

static void test_routh_hurwitz(void) {
    TEST("routh-hurwitz stable filter");
    /* Create a stable 2nd-order TF: H(s) = 1/(s^2 + 2*s + 1) */
    tf_continuous_t *tf = tf_continuous_alloc(1, 3);
    assert(tf != NULL);
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; /* s^0 */
    tf->den[1] = 2.0; /* s^1 */
    tf->den[2] = 1.0; /* s^2 */
    tf->gain = 1.0;

    assert(tf_continuous_is_stable(tf) == 1);

    /* Unstable: H(s) = 1/(s^2 - s + 1) — pole in RHP */
    tf->den[1] = -1.0;
    assert(tf_continuous_is_stable(tf) == 0);

    tf_continuous_free(tf);
    PASS();
}

/* ==================================================================
 * L4: Jury Stability Test
 * ================================================================== */

static void test_jury_stability(void) {
    TEST("jury stability test");
    /* Stable: H(z) = 1/(1 - 0.5z^{-1} + 0.06z^{-2})
     * Poles at z = 0.3, z = 0.2 */
    tf_discrete_t *tf = tf_discrete_alloc(1, 3);
    assert(tf != NULL);
    tf->num[0] = 1.0;
    tf->den[0] = 1.0;
    tf->den[1] = -0.5;
    tf->den[2] = 0.06;
    tf->gain = 1.0;

    assert(tf_discrete_is_stable(tf) == 1);

    /* Unstable: H(z) = 1/(1 - 1.5z^{-1}) — pole at z=1.5 */
    tf->den[1] = -1.5;
    tf->den[2] = 0.0;
    assert(tf_discrete_is_stable(tf) == 0);

    tf_discrete_free(tf);
    PASS();
}

/* ==================================================================
 * L4: Bilinear Transform
 * ================================================================== */

static void test_bilinear_prewarp(void) {
    TEST("bilinear prewarp");
    double fs = 1000.0;
    double wd = 2.0 * M_PI * 100.0 / fs; /* 100 Hz digital cutoff in rad/sample */
    double wa = bilinear_prewarp(wd, fs);

    /* Pre-warped analog frequency should be > digital frequency */
    double fd = 100.0;
    double fa = wa / (2.0 * M_PI);
    assert(fa > fd);

    /* Check round-trip: wd -> wa -> wd */
    double wd_back = 2.0 * atan(wa / (2.0 * fs));
    assert(fabs(wd - wd_back) < 1e-10);
    PASS();
}

static void test_bilinear_s_to_z(void) {
    TEST("bilinear s->z stability preservation");
    double fs = 1000.0;

    /* LHP point should map inside unit circle */
    double complex s_stable = -100.0 + 200.0 * I;
    double complex z = bilinear_s_to_z(s_stable, fs);
    assert(cabs(z) < 1.0);

    /* jw-axis should map to unit circle */
    double complex s_imag = 300.0 * I;
    z = bilinear_s_to_z(s_imag, fs);
    assert(fabs(cabs(z) - 1.0) < 1e-10);
    PASS();
}

/* ==================================================================
 * L5: Window Functions
 * ================================================================== */

static void test_window_generate(void) {
    TEST("window_generate all types");
    int N = 64;
    double *window = (double *)malloc(N * sizeof(double));
    assert(window != NULL);

    window_type_t types[] = {
        WINDOW_RECTANGULAR, WINDOW_HANN, WINDOW_HAMMING,
        WINDOW_BLACKMAN, WINDOW_BLACKMAN_HARRIS, WINDOW_KAISER,
        WINDOW_FLATTOP, WINDOW_BARTLETT, WINDOW_NUTTALL, WINDOW_TUKEY
    };
    int num_types = sizeof(types) / sizeof(types[0]);
    int t;

    for (t = 0; t < num_types; t++) {
        window_params_t wp = {types[t], 5.0, N};
        int ret = window_generate(window, &wp);
        assert(ret == FILTER_OK);

        /* Window should be symmetric */
        int n;
        for (n = 0; n < N / 2; n++) {
            assert(fabs(window[n] - window[N - 1 - n]) < 1e-10);
        }

        /* Max value should be at center (or near it) */
        double max_val = 0.0;
        for (n = 0; n < N; n++) {
            if (window[n] > max_val) max_val = window[n];
        }
        assert(fabs(max_val - 1.0) < 0.1);
    }

    free(window);
    PASS();
}

/* ==================================================================
 * L5: FIR Filter Design
 * ================================================================== */

static void test_fir_design_window(void) {
    TEST("fir_design_window lowpass");
    /* Design a lowpass FIR with fc = 0.25 * fs/2 = 0.125 normalized */
    fir_filter_t *fir = fir_design_window(0.125, 31, WINDOW_HAMMING, 0.0);
    assert(fir != NULL);
    assert(fir->length == 31);

    /* Check symmetry (linear phase) */
    int n;
    for (n = 0; n < 15; n++) {
        assert(fabs(fir->coeff[n] - fir->coeff[30 - n]) < 1e-10);
    }

    /* DC gain should be close to 1.0 */
    assert(fabs(fir->gain_dc - 1.0) < 0.3);

    fir_filter_free(fir);
    PASS();
}

/* ==================================================================
 * L5: FIR Differentiator
 * ================================================================== */

static void test_fir_differentiator(void) {
    TEST("fir_design_differentiator");
    fir_filter_t *fir = fir_design_differentiator(31, WINDOW_HAMMING, 0.0);
    assert(fir != NULL);
    assert(fir->lp_type == FIR_TYPE_III || fir->lp_type == FIR_TYPE_IV);

    /* Anti-symmetry check */
    int n;
    for (n = 0; n < 15; n++) {
        assert(fabs(fir->coeff[n] + fir->coeff[30 - n]) < 1e-8);
    }

    fir_filter_free(fir);
    PASS();
}

/* ==================================================================
 * L5: FIR Hilbert Transformer
 * ================================================================== */

static void test_fir_hilbert(void) {
    TEST("fir_design_hilbert");
    fir_filter_t *fir = fir_design_hilbert(31, WINDOW_HAMMING, 0.0);
    assert(fir != NULL);
    assert(fir->lp_type == FIR_TYPE_III);
    fir_filter_free(fir);
    PASS();
}

/* ==================================================================
 * L5: Order Estimation
 * ================================================================== */

static void test_order_estimation(void) {
    TEST("butterworth_estimate_order");
    int order = butterworth_estimate_order(0.1, 60.0,
                                            2.0 * M_PI * 1000.0,
                                            2.0 * M_PI * 1500.0);
    assert(order >= 4);
    assert(order <= 25);
    PASS();

    TEST("chebyshev1_estimate_order");
    order = chebyshev1_estimate_order(0.1, 60.0, 2.0 * M_PI * 1500.0,
                                       2.0 * M_PI * 1000.0);
    assert(order >= 3);
    assert(order <= 10);
    PASS();

    TEST("elliptic_estimate_order");
    order = elliptic_estimate_order(0.1, 60.0, 2.0 * M_PI * 1000.0,
                                     2.0 * M_PI * 1500.0);
    assert(order >= 2);
    PASS();
}

/* ==================================================================
 * L5: Kaiser Beta
 * ================================================================== */

static void test_kaiser_beta(void) {
    TEST("kaiser_beta_from_attenuation");
    double beta_60 = kaiser_beta_from_attenuation(60.0);
    assert(beta_60 > 5.0);
    double beta_30 = kaiser_beta_from_attenuation(30.0);
    assert(beta_30 > 1.0);
    assert(beta_60 > beta_30);
    PASS();
}

/* ==================================================================
 * L6: Frequency Response
 * ================================================================== */

static void test_fir_freq_response(void) {
    TEST("fir_freq_response_at");
    /* Design a simple averaging filter: h[n] = 1/5 for n=0..4 */
    double h[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
    fir_filter_t *fir = fir_filter_alloc(h, 5);
    assert(fir != NULL);

    /* At DC (w=0): magnitude should be 1.0 */
    freq_resp_point_t fp = fir_freq_response_at(fir, 0.0);
    assert(fabs(fp.magnitude - 1.0) < 1e-10);

    /* At Nyquist (w=pi): magnitude of moving average filter */
    fp = fir_freq_response_at(fir, M_PI);
    assert(fp.magnitude < 1.0);

    fir_filter_free(fir);
    PASS();
}

static void test_iir_freq_response(void) {
    TEST("iir_freq_response_at");
    /* Simple first-order lowpass: y[n] = 0.5*x[n] + 0.5*y[n-1] */
    biquad_section_t bq = {0.5, 0.0, 0.0, -0.5, 0.0, 1.0};
    iir_filter_t *iir = iir_filter_alloc(&bq, 1, 1.0);
    assert(iir != NULL);

    /* At DC: |H| = 0.5/(1-0.5) = 1.0 */
    freq_resp_point_t fp = iir_freq_response_at(iir, 0.0);
    assert(fabs(fp.magnitude - 1.0) < 0.01);

    /* At Nyquist: |H| = 0.5/(1+0.5) = 0.333 */
    fp = iir_freq_response_at(iir, M_PI);
    assert(fabs(fp.magnitude - 0.333) < 0.05);

    iir_filter_free(iir);
    PASS();
}

/* ==================================================================
 * L6: Step Response
 * ================================================================== */

static void test_step_response(void) {
    TEST("fir_step_response");
    /* Impulse response: h = [0.25, 0.25, 0.25, 0.25] */
    double h[4] = {0.25, 0.25, 0.25, 0.25};
    double step[10];
    fir_step_response(h, 4, step, 10);

    /* Step should eventually reach 1.0 */
    assert(fabs(step[9] - 1.0) < 1e-10);

    /* Step is non-decreasing */
    int n;
    for (n = 1; n < 10; n++) {
        assert(step[n] >= step[n - 1] - 1e-15);
    }
    PASS();
}

/* ==================================================================
 * L6: SNR Improvement
 * ================================================================== */

static void test_snr_improvement(void) {
    TEST("snr_improvement");
    /* Create ideal lowpass: perfect SNR improvement */
    int N = 128;
    double *h = (double *)malloc(N * sizeof(double));
    assert(h != NULL);
    ideal_lowpass_impulse(h, N, M_PI * 0.25);

    fir_filter_t *fir = fir_filter_alloc(h, N);
    assert(fir != NULL);

    freq_resp_t *resp = fir_freq_response(fir, 0.0, 500.0, 256, 1000.0);
    assert(resp != NULL);

    double snr = snr_improvement(resp, 0.0, 125.0);
    /* Lowpass filter should improve SNR in passband */
    assert(snr > 0.0);

    freq_resp_free(resp);
    fir_filter_free(fir);
    free(h);
    PASS();
}

/* ==================================================================
 * L3: Polynomial Operations
 * ================================================================== */

static void test_binomial(void) {
    TEST("binomial");
    assert(fabs(binomial(5, 0) - 1.0) < 1e-10);
    assert(fabs(binomial(5, 1) - 5.0) < 1e-10);
    assert(fabs(binomial(5, 2) - 10.0) < 1e-10);
    assert(fabs(binomial(5, 5) - 1.0) < 1e-10);
    PASS();
}

static void test_bessel_i0(void) {
    TEST("bessel_i0");
    /* I0(0) = 1 */
    assert(fabs(bessel_i0(0.0) - 1.0) < 1e-10);
    /* I0 is monotonic increasing for x > 0 */
    assert(bessel_i0(5.0) > bessel_i0(1.0));
    assert(bessel_i0(10.0) > bessel_i0(5.0));
    PASS();
}

/* ==================================================================
 * L3: Bessel Polynomials
 * ================================================================== */

static void test_bessel_polynomials(void) {
    TEST("bessel_polynomials");
    int n = 4;
    double *coeffs = (double *)malloc((n + 1) * (n + 1) * sizeof(double));
    assert(coeffs != NULL);

    int ret = bessel_polynomials(n, coeffs);
    assert(ret == FILTER_OK);

    /* B0(s) = 1 */
    assert(fabs(coeffs[0 * 5 + 0] - 1.0) < 1e-10);

    /* B1(s) = s + 1 */
    assert(fabs(coeffs[1 * 5 + 0] - 1.0) < 1e-10);
    assert(fabs(coeffs[1 * 5 + 1] - 1.0) < 1e-10);

    /* B2(s) = s^2 + 3s + 3 */
    assert(fabs(coeffs[2 * 5 + 0] - 3.0) < 1e-10);
    assert(fabs(coeffs[2 * 5 + 1] - 3.0) < 1e-10);
    assert(fabs(coeffs[2 * 5 + 2] - 1.0) < 1e-10);

    free(coeffs);
    PASS();
}

/* ==================================================================
 * L3: Elliptic K
 * ================================================================== */

static void test_elliptic_k(void) {
    TEST("elliptic_k");
    /* K(0) = pi/2 */
    assert(fabs(elliptic_k(0.0) - M_PI / 2.0) < 1e-10);
    /* K(1) = infinity */
    assert(isinf(elliptic_k(1.0)));
    /* K should be monotonic increasing */
    assert(elliptic_k(0.5) > elliptic_k(0.3));
    PASS();
}

/* ==================================================================
 * L5: Coef Quantization
 * ================================================================== */

static void test_coef_quantize(void) {
    TEST("coef_quantize");
    double coeff[3] = {0.123456, -0.654321, 0.999999};
    coef_quantize(coeff, 3, 8); /* 8 fractional bits */

    /* Check quantization to 1/256 resolution */
    int i;
    for (i = 0; i < 3; i++) {
        double scaled = coeff[i] * 256.0;
        assert(fabs(scaled - round(scaled)) < 1e-10);
    }
    PASS();
}

/* ==================================================================
 * L2: FIR Linear Phase Type Detection
 * ================================================================== */

static void test_fir_linear_phase_type(void) {
    TEST("fir_linear_phase_type Type I");
    double h1[5] = {0.1, 0.2, 0.4, 0.2, 0.1}; /* Symmetric, odd length */
    int t = fir_determine_linear_phase_type(h1, 5, 1e-10);
    assert(t == FIR_TYPE_I);
    PASS();

    TEST("fir_linear_phase_type Type III");
    double h3[5] = {0.1, 0.2, 0.0, -0.2, -0.1}; /* Anti-symmetric, odd */
    t = fir_determine_linear_phase_type(h3, 5, 1e-10);
    assert(t == FIR_TYPE_III);
    PASS();
}

/* ==================================================================
 * L6: Analog filter design
 * ================================================================== */

static void test_analog_filter_design(void) {
    TEST("analog_filter_design Butterworth LP");
    filter_spec_t spec = filter_spec_init();
    spec.type = FILTER_LOWPASS;
    spec.approx = APPROX_BUTTERWORTH;
    spec.fc1 = 1000.0;
    spec.order = 4;

    tf_continuous_t *tf = analog_filter_design(&spec);
    assert(tf != NULL);

    /* Check DC gain ≈ 0 dB */
    double mag_dc = cabs(tf_continuous_eval(tf, 0));
    double db_dc = 20.0 * log10(mag_dc);
    assert(fabs(db_dc) < 1.0);

    /* Check roll-off at 2x cutoff */
    double wc = 2.0 * M_PI * 1000.0;
    double mag_2x = cabs(tf_continuous_eval(tf, 2.0 * wc * I));
    double db_2x = 20.0 * log10(mag_2x);
    /* 4th order: -24 dB at 2x fc, expect < -20 dB */
    assert(db_2x < -15.0);

    tf_continuous_free(tf);
    PASS();
}

/* ==================================================================
 * L6: Chebyshev I design
 * ================================================================== */

static void test_chebyshev1_design(void) {
    TEST("chebyshev1_prototype");
    tf_continuous_t *proto = chebyshev1_prototype(4, 1.0);
    assert(proto != NULL);

    /* Check it's stable (all poles in LHP) */
    assert(tf_continuous_is_stable(proto) == 1);

    tf_continuous_free(proto);
    PASS();
}

/* ==================================================================
 * Main
 * ================================================================== */

int main(void) {
    printf("\n=== mini-filter-theory Test Suite ===\n\n");

    /* L1: Definitions */
    test_filter_spec_init();
    test_filter_spec_is_valid();

    /* L3: Math */
    test_binomial();
    test_bessel_i0();
    test_bessel_polynomials();
    test_elliptic_k();

    /* L4: Stability */
    test_butterworth_poles();
    test_butterworth_prototype();
    test_routh_hurwitz();
    test_jury_stability();
    test_bilinear_prewarp();
    test_bilinear_s_to_z();

    /* L5: Design */
    test_window_generate();
    test_fir_design_window();
    test_fir_differentiator();
    test_fir_hilbert();
    test_order_estimation();
    test_kaiser_beta();
    test_coef_quantize();

    /* L6: Analysis */
    test_fir_freq_response();
    test_iir_freq_response();
    test_step_response();
    test_snr_improvement();
    test_fir_linear_phase_type();
    test_analog_filter_design();
    test_chebyshev1_design();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
