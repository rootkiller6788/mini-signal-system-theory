/**
 * @file test_stability.c
 * @brief Tests for transfer function, stability, and state-space modules
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/system_defs.h"
#include "../include/transfer_function.h"
#include "../include/stability.h"
#include "../include/state_space.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TOL 1e-8

static int passed = 0, failed = 0;
#define T(name) printf("  TEST: %s ... ", name)
#define P() do { printf("PASS\n"); passed++; } while(0)
#define F(m) do { printf("FAIL: %s\n", m); failed++; } while(0)
#define C(c,m) do { if(!(c)){F(m);return;} } while(0)

/* ---- TF: First-order system evaluation ---- */
static void test_tf_evaluate(void)
{
    T("TF evaluate H(s)=1/(s+1) at s=0");
    ct_transfer_function_t tf = ct_tf_alloc(0, 1);
    tf.num_coeffs[0] = 1.0;
    tf.den_coeffs[0] = 1.0; tf.den_coeffs[1] = 1.0;

    double complex H0 = ct_tf_evaluate(&tf, 0.0 + 0.0 * I);
    C(fabs(creal(H0) - 1.0) < TOL && fabs(cimag(H0)) < TOL, "DC gain != 1");

    double complex H1 = ct_tf_evaluate(&tf, 1.0 * I);
    /* H(j1) = 1/(1+j) = (1-j)/2 = 0.5 - 0.5j */
    C(fabs(creal(H1) - 0.5) < TOL, "Re{H(j1)} wrong");
    C(fabs(cimag(H1) + 0.5) < TOL, "Im{H(j1)} wrong");

    ct_tf_free(&tf);
    P();
}

/* ---- TF: Frequency response ---- */
static void test_frequency_response(void)
{
    T("frequency response of 1/(s+1)");
    ct_transfer_function_t tf = ct_tf_alloc(0, 1);
    tf.num_coeffs[0] = 1.0;
    tf.den_coeffs[0] = 1.0; tf.den_coeffs[1] = 1.0;

    frequency_response_t fr = freq_resp_alloc(3, 0.01, 10.0, 1);
    ct_frequency_response(&tf, &fr);

    /* At DC: 0 dB */
    C(fabs(fr.points[0].magnitude_db) < 1.0, "DC gain not 0 dB");
    /* At high freq: should be attenuated */
    C(fr.points[2].magnitude_db < fr.points[0].magnitude_db, "not attenuating");

    freq_resp_free(&fr);
    ct_tf_free(&tf);
    P();
}

/* ---- TF: Pole-Zero analysis ---- */
static void test_pole_zero(void)
{
    T("pole-zero: 1/(s^2 + 2*zeta*wn*s + wn^2)");
    double zeta = 0.5, wn = 2.0;
    ct_transfer_function_t tf = ct_tf_alloc(0, 2);
    tf.num_coeffs[0] = wn * wn;
    tf.den_coeffs[0] = wn * wn;
    tf.den_coeffs[1] = 2.0 * zeta * wn;
    tf.den_coeffs[2] = 1.0;

    pole_zero_collection_t pz = pz_collection_alloc(2, 0);
    ct_tf_find_poles(&tf, &pz);

    C(pz.num_poles == 2, "wrong number of poles");
    C(pz.is_stable == 1, "stable system flagged unstable");
    /* Poles at -zeta*wn +/- j*wn*sqrt(1-zeta^2) = -1 +/- j*sqrt(3) */
    double expected_wn = wn;  /* 2.0 */
    C(fabs(pz.poles[0].wn - expected_wn) < 0.2, "wrong natural frequency");

    pz_collection_free(&pz);
    ct_tf_free(&tf);
    P();
}

/* ---- Stability: Routh-Hurwitz ---- */
static void test_routh_hurwitz(void)
{
    T("Routh-Hurwitz: s^3 + s^2 + 2s + 24 (unstable)");
    double coeffs[] = {24.0, 2.0, 1.0, 1.0};  /* a0...a3 */
    int special;
    double rh[20];
    int rhp = routh_hurwitz(coeffs, 3, rh, &special);
    /* s^3+s^2+2s+24 has RHP roots */
    C(rhp > 0, "should detect RHP poles");
    P();
}

static void test_routh_hurwitz_stable(void)
{
    T("Routh-Hurwitz: s^2 + 3s + 2 = (s+1)(s+2) (stable)");
    double coeffs[] = {2.0, 3.0, 1.0};
    double rh[15];
    int special;
    int rhp = routh_hurwitz(coeffs, 2, rh, &special);
    C(rhp == 0, "stable system flagged unstable");
    P();
}

/* ---- Stability: BIBO ---- */
static void test_bibo_stable(void)
{
    T("BIBO: 1/(s+1) is stable");
    ct_transfer_function_t tf = ct_tf_alloc(0, 1);
    tf.num_coeffs[0] = 1.0;
    tf.den_coeffs[0] = 1.0; tf.den_coeffs[1] = 1.0;
    C(ct_is_bibo_stable(&tf) == 1, "1/(s+1) should be stable");
    ct_tf_free(&tf);
    P();
}

static void test_bibo_unstable(void)
{
    T("BIBO: 1/(s-1) is unstable");
    ct_transfer_function_t tf = ct_tf_alloc(0, 1);
    tf.num_coeffs[0] = 1.0;
    tf.den_coeffs[0] = -1.0; tf.den_coeffs[1] = 1.0;
    /* This is 1/(s-1) = pole at +1 (RHP) */
    C(ct_is_bibo_stable(&tf) == 0, "1/(s-1) should be unstable");
    ct_tf_free(&tf);
    P();
}

/* ---- Stability: Settling time and overshoot ---- */
static void test_settling_overshoot(void)
{
    T("settling time and overshoot estimates");
    double zeta = 0.5, wn = 2.0;
    double ts = ct_settling_time_estimate(zeta, wn);
    double Mp = ct_overshoot_estimate(zeta);

    /* ts = 4/(0.5*2) = 4.0 */
    C(fabs(ts - 4.0) < TOL, "wrong settling time");
    /* Mp = exp(-pi*0.5/sqrt(0.75)) = exp(-1.814) ≈ 16.3% */
    double expected_Mp = 100.0 * exp(-M_PI * 0.5 / sqrt(0.75));
    C(fabs(Mp - expected_Mp) < 1e-6, "wrong overshoot");

    P();
}

/* ---- State-Space: Controllability ---- */
static void test_controllability(void)
{
    T("state-space controllability");
    ct_state_space_t ss = ct_ss_alloc(2, 1, 1);
    /* Double integrator: A=[0 1; 0 0], B=[0; 1] */
    ss.A[0*2+1] = 1.0;
    ss.B[1] = 1.0;

    int ctrl = ct_ss_is_controllable(&ss);
    C(ctrl == 1, "double integrator should be controllable");

    ct_ss_free(&ss);
    P();
}

/* ---- Jury Test ---- */
static void test_jury_stable(void)
{
    T("Jury: z^2 - 0.5z + 0.06 (stable)");
    double coeffs[] = {0.06, -0.5, 1.0};  /* a0, a1, a2 */
    int stable = jury_stability_test(coeffs, 2);
    C(stable == 1, "stable DT system flagged unstable");
    P();
}

int main(void)
{
    printf("=== mini-system-analysis Stability Test Suite ===\n\n");

    test_tf_evaluate();
    test_frequency_response();
    test_pole_zero();
    test_routh_hurwitz();
    test_routh_hurwitz_stable();
    test_bibo_stable();
    test_bibo_unstable();
    test_settling_overshoot();
    test_controllability();
    test_jury_stable();

    printf("\n=== Results: %d passed, %d failed ===\n", passed, failed);
    return (failed > 0) ? 1 : 0;
}
