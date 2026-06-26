/**
 * @file test_system.c
 * @brief Test suite for system_defs, convolution, and frequency modules
 *
 * All tests use standard assert() from <assert.h>.
 * Run: make test
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/system_defs.h"
#include "../include/convolution.h"

#define TOL 1e-10
#define EPS 1e-6

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

/* ---- L1: Signal Lifecycle ---- */
static void test_ct_signal_alloc_free(void)
{
    TEST("ct_signal_alloc/free");
    ct_signal_t sig = ct_signal_alloc(100, 0.0, 1.0, 100.0);
    CHECK(sig.samples != NULL, "allocation failed");
    CHECK(sig.num_samples == 100, "wrong num_samples");
    CHECK(fabs(sig.t_start - 0.0) < TOL, "wrong t_start");
    CHECK(fabs(sig.t_end - 1.0) < TOL, "wrong t_end");
    CHECK(sig.owns_buffer == 1, "owns_buffer not set");
    ct_signal_free(&sig);
    CHECK(sig.samples == NULL, "free did not nullify pointer");
    PASS();
}

static void test_dt_signal_alloc_free(void)
{
    TEST("dt_signal_alloc/free");
    dt_signal_t sig = dt_signal_alloc(50);
    CHECK(sig.samples != NULL, "allocation failed");
    CHECK(sig.num_samples == 50, "wrong num_samples");
    dt_signal_free(&sig);
    CHECK(sig.samples == NULL, "free did not nullify pointer");
    PASS();
}

/* ---- L1: Signal Energy ---- */
static void test_signal_energy(void)
{
    TEST("dt_signal_energy");
    dt_signal_t sig = dt_signal_alloc(4);
    sig.samples[0] = 1.0;
    sig.samples[1] = 2.0;
    sig.samples[2] = 3.0;
    sig.samples[3] = 4.0;
    double energy = dt_signal_energy(&sig);
    double expected = 1.0 + 4.0 + 9.0 + 16.0;  /* 30.0 */
    CHECK(fabs(energy - expected) < TOL, "wrong energy");
    dt_signal_free(&sig);
    PASS();
}

/* ---- L1: Unit Impulse/Step ---- */
static void test_unit_impulse(void)
{
    TEST("unit impulse");
    dt_signal_t sig = dt_signal_alloc(5);
    dt_signal_set_unit_impulse(&sig, 2);
    CHECK(sig.samples[2] == 1.0, "impulse not at correct index");
    CHECK(sig.samples[0] == 0.0, "non-zero before impulse");
    dt_signal_free(&sig);
    PASS();
}

static void test_unit_step(void)
{
    TEST("unit step");
    dt_signal_t sig = dt_signal_alloc(5);
    dt_signal_set_unit_step(&sig, 2);
    CHECK(sig.samples[0] == 0.0, "non-zero before step");
    CHECK(sig.samples[1] == 0.0, "non-zero before step");
    CHECK(sig.samples[2] == 1.0, "step not at correct index");
    CHECK(sig.samples[4] == 1.0, "not sustained");
    dt_signal_free(&sig);
    PASS();
}

/* ---- L2: DT Convolution ---- */
static void test_dt_convolution_direct(void)
{
    TEST("dt_convolution_direct (simple)");
    dt_signal_t x = dt_signal_alloc(3);
    dt_signal_t h = dt_signal_alloc(2);
    dt_signal_t y = dt_signal_alloc(4);
    x.samples[0] = 1.0; x.samples[1] = 2.0; x.samples[2] = 3.0;
    h.samples[0] = 1.0; h.samples[1] = 1.0;

    int ret = dt_convolution_direct(&x, &h, &y);
    CHECK(ret == 0, "convolution failed");

    /* y = [1, 3, 5, 3] */
    CHECK(fabs(y.samples[0] - 1.0) < TOL, "y[0] wrong");
    CHECK(fabs(y.samples[1] - 3.0) < TOL, "y[1] wrong");
    CHECK(fabs(y.samples[2] - 5.0) < TOL, "y[2] wrong");
    CHECK(fabs(y.samples[3] - 3.0) < TOL, "y[3] wrong");

    dt_signal_free(&x); dt_signal_free(&h); dt_signal_free(&y);
    PASS();
}

/* ---- L2: Convolution Commutativity ---- */
static void test_convolution_commutativity(void)
{
    TEST("convolution commutativity: x*h == h*x");
    dt_signal_t x = dt_signal_alloc(20);
    dt_signal_t h = dt_signal_alloc(10);
    for (size_t i = 0; i < 20; i++) x.samples[i] = sin((double)i * 0.5);
    for (size_t i = 0; i < 10; i++) h.samples[i] = exp(-(double)i * 0.3);

    double diff = dt_convolution_commutativity_check(&x, &h);
    CHECK(diff < 1e-6, "commutativity violation");

    dt_signal_free(&x); dt_signal_free(&h);
    PASS();
}

/* ---- L2: Step from Impulse ---- */
static void test_step_from_impulse(void)
{
    TEST("step response from impulse response");
    dt_signal_t h = dt_signal_alloc(5);
    dt_signal_t s = dt_signal_alloc(5);
    h.samples[0]=0.5; h.samples[1]=0.3; h.samples[2]=0.1; h.samples[3]=0.05; h.samples[4]=0.02;

    dt_step_from_impulse(&h, &s);
    /* s[0]=0.5, s[1]=0.8, s[2]=0.9, s[3]=0.95, s[4]=0.97 */
    CHECK(fabs(s.samples[0] - 0.5) < TOL, "s[0]");
    CHECK(fabs(s.samples[1] - 0.8) < TOL, "s[1]");
    CHECK(fabs(s.samples[2] - 0.9) < TOL, "s[2]");

    dt_signal_free(&h); dt_signal_free(&s);
    PASS();
}

/* ---- L2: Auto-Correlation ---- */
static void test_auto_correlation(void)
{
    TEST("auto-correlation: Rxx[0] = energy");
    dt_signal_t x = dt_signal_alloc(3);
    x.samples[0]=1.0; x.samples[1]=2.0; x.samples[2]=3.0;
    dt_signal_t rxx = dt_signal_alloc(5);

    dt_auto_correlation(&x, &rxx);
    /* Rxx[0] = 1^2+2^2+3^2 = 14 */
    CHECK(fabs(rxx.samples[2] - 14.0) < TOL, "Rxx[0] != energy");

    dt_signal_free(&x); dt_signal_free(&rxx);
    PASS();
}

/* ---- L2: Cross-Correlation ---- */
static void test_cross_correlation(void)
{
    TEST("cross-correlation of identical signals");
    dt_signal_t x = dt_signal_alloc(4);
    x.samples[0]=1.0; x.samples[1]=0.0; x.samples[2]=0.0; x.samples[3]=0.0;

    dt_signal_t r = dt_signal_alloc(7);
    dt_cross_correlation(&x, &x, &r);

    /* Rxx[0] = 1, others=0 */
    CHECK(fabs(r.samples[3] - 1.0) < TOL, "cross-correlation wrong");

    dt_signal_free(&x); dt_signal_free(&r);
    PASS();
}

int main(void)
{
    printf("=== mini-system-analysis Test Suite ===\n\n");

    test_ct_signal_alloc_free();
    test_dt_signal_alloc_free();
    test_signal_energy();
    test_unit_impulse();
    test_unit_step();
    test_dt_convolution_direct();
    test_convolution_commutativity();
    test_step_from_impulse();
    test_auto_correlation();
    test_cross_correlation();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return (tests_failed > 0) ? 1 : 0;
}
