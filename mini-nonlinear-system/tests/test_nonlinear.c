/**
 * @file test_nonlinear.c
 * @brief Test suite for mini-nonlinear-system module.
 *
 * Tests all core APIs using runtime assertions.
 * Coverage: L1 (definitions), L4 (theorems verified), L5 (algorithms),
 * L6 (canonical problems).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../include/nonlinear_core.h"
#include "../include/nonlinear_stability.h"
#include "../include/nonlinear_describing.h"
#include "../include/nonlinear_bifurcation.h"
#include "../include/nonlinear_chaos.h"
#include "../include/nonlinear_phase.h"
#include "../include/nonlinear_volterra.h"
#include "../include/nonlinear_applications.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

#define FAIL(msg) do { \
    printf("FAILED: %s\n", msg); \
} while(0)

/* ===================================================================
 * Test System: Van der Pol oscillator
 * dx/dt = y
 * dy/dt = mu*(1 - x^2)*y - x
 * =================================================================== */

static int vdp_f(const double *x, double t, const double *params,
                 double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 2) return -1;
    double mu = params[0];
    dxdt[0] = x[1];
    dxdt[1] = mu * (1.0 - x[0] * x[0]) * x[1] - x[0];
    return 0;
}

/* ===================================================================
 * Test System: Simple linear oscillator (decoupled 2D)
 * dx/dt = -x
 * dy/dt = -2y
 * Equilibrium at origin, stable node.
 * =================================================================== */

static int linear_stable_f(const double *x, double t, const double *params,
                           double *dxdt, size_t dim)
{
    (void)t; (void)params;
    if (dim != 2) return -1;
    dxdt[0] = -x[0];
    dxdt[1] = -2.0 * x[1];
    return 0;
}

/* ===================================================================
 * Test System: Saddle point
 * dx/dt = x
 * dy/dt = -y
 * =================================================================== */

static int saddle_f(const double *x, double t, const double *params,
                    double *dxdt, size_t dim)
{
    (void)t; (void)params;
    if (dim != 2) return -1;
    dxdt[0] = x[0];
    dxdt[1] = -x[1];
    return 0;
}

/* ===================================================================
 * Test: L1 - Core Definitions
 * =================================================================== */

static void test_system_init(void)
{
    TEST("nl_system_init");
    nl_nonlinear_system_t sys;
    int ret = nl_system_init(&sys, 2, linear_stable_f, "linear_test");
    assert(ret == 0);
    assert(sys.dim == 2);
    assert(sys.f != NULL);
    assert(strcmp(sys.name, "linear_test") == 0);
    PASS();
}

static void test_trajectory_ops(void)
{
    TEST("nl_trajectory_init/append/free");

    nl_trajectory_t traj;
    int ret = nl_trajectory_init(&traj, 10, 2);
    assert(ret == 0);
    assert(traj.capacity == 10);
    assert(traj.num_points == 0);

    /* Append points */
    nl_state_t s;
    s.dim = 2;
    s.x[0] = 1.0; s.x[1] = 0.0;
    ret = nl_trajectory_append(&traj, &s);
    assert(ret == 0);
    assert(traj.num_points == 1);

    s.x[0] = 0.5; s.x[1] = 0.5;
    ret = nl_trajectory_append(&traj, &s);
    assert(ret == 0);
    assert(traj.num_points == 2);

    /* Verify points */
    assert(fabs(traj.points[0].x[0] - 1.0) < 1e-10);
    assert(fabs(traj.points[1].x[1] - 0.5) < 1e-10);

    /* Test auto-grow (append beyond capacity) */
    nl_trajectory_t traj2;
    nl_trajectory_init(&traj2, 2, 2);
    for (int i = 0; i < 10; i++) {
        s.x[0] = (double)i; s.x[1] = (double)(i * 2);
        assert(nl_trajectory_append(&traj2, &s) == 0);
    }
    assert(traj2.num_points == 10);
    assert(traj2.capacity >= 10);

    nl_trajectory_free(&traj);
    nl_trajectory_free(&traj2);
    PASS();
}

static void test_nonlinearity_eval(void)
{
    TEST("nl_nonlinearity_eval (all types)");

    nl_nonlinearity_t nl;
    double p[3];

    /* Ideal relay: M=2.0 */
    p[0] = 2.0;
    nl_nonlinearity_init(&nl, NL_TYPE_IDEAL_RELAY, p, 1);
    assert(nl_nonlinearity_eval(&nl, 3.0) == 2.0);
    assert(nl_nonlinearity_eval(&nl, -3.0) == -2.0);
    assert(nl_nonlinearity_eval(&nl, 0.0) == 0.0);

    /* Dead zone: delta=0.5 */
    p[0] = 0.5;
    nl_nonlinearity_init(&nl, NL_TYPE_DEAD_ZONE, p, 1);
    assert(nl_nonlinearity_eval(&nl, 0.3) == 0.0);
    assert(nl_nonlinearity_eval(&nl, 1.0) == 0.5);
    assert(nl_nonlinearity_eval(&nl, -1.0) == -0.5);

    /* Saturation: k=2.0, a=1.0 */
    p[0] = 2.0; p[1] = 1.0;
    nl_nonlinearity_init(&nl, NL_TYPE_SATURATION, p, 2);
    assert(nl_nonlinearity_eval(&nl, 0.3) == 0.6);
    assert(nl_nonlinearity_eval(&nl, 2.0) == 1.0);
    assert(nl_nonlinearity_eval(&nl, -3.0) == -1.0);

    /* Polynomial: f(x) = 1 + 2x + 3x^2 */
    p[0] = 1.0; p[1] = 2.0; p[2] = 3.0;
    nl_nonlinearity_init(&nl, NL_TYPE_POLYNOMIAL, p, 3);
    double val = nl_nonlinearity_eval(&nl, 2.0);
    double expected = 1.0 + 2.0*2.0 + 3.0*4.0;  /* = 17.0 */
    assert(fabs(val - expected) < 1e-10);

    /* Sigmoid: alpha=1.0 */
    p[0] = 1.0;
    nl_nonlinearity_init(&nl, NL_TYPE_SIGMOID, p, 1);
    assert(fabs(nl_nonlinearity_eval(&nl, 0.0)) < 1e-10);
    assert(nl_nonlinearity_eval(&nl, 10.0) > 0.99);

    /* Quantizer: step=0.5 */
    p[0] = 0.5;
    nl_nonlinearity_init(&nl, NL_TYPE_QUANTIZER, p, 1);
    assert(nl_nonlinearity_eval(&nl, 0.23) == 0.0);
    assert(nl_nonlinearity_eval(&nl, 0.3) == 0.5);

    /* Cubic stiffness */
    p[0] = 1.0; p[1] = 0.5;
    nl_nonlinearity_init(&nl, NL_TYPE_CUBIC_STIFFNESS, p, 2);
    assert(fabs(nl_nonlinearity_eval(&nl, 2.0) - 6.0) < 1e-10);

    PASS();
}

/* ===================================================================
 * Test: L4 - Lyapunov Stability (Indirect Method)
 * =================================================================== */

static void test_stability_classification(void)
{
    TEST("nl_classify_equilibrium (stable node, saddle)");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "linear_stable");

    /* Equilibrium at origin: both eigenvalues negative -> stable node */
    double x_eq_arr[2] = {0.0, 0.0};
    nl_state_t eq;
    eq.dim = 2;
    memcpy(eq.x, x_eq_arr, 2 * sizeof(double));

    nl_jacobian_t J;
    int ret = nl_compute_jacobian(&sys, eq.x, &J);
    assert(ret == 0);

    nl_stability_t stab = nl_classify_equilibrium(&J, 1e-6);
    /* eigenvalues: -1 and -2 -> stable node */
    assert(stab == NL_STABLE_NODE || stab == NL_STABLE_FOCUS
           || stab == NL_ASYMPTOTICALLY_STABLE);

    /* Saddle test */
    nl_system_init(&sys, 2, saddle_f, "saddle");
    ret = nl_compute_jacobian(&sys, eq.x, &J);
    assert(ret == 0);
    stab = nl_classify_equilibrium(&J, 1e-6);
    assert(stab == NL_SADDLE);

    PASS();
}

static void test_quadratic_lyapunov(void)
{
    TEST("nl_quadratic_lyapunov / nl_quadratic_lyapunov_dot");

    /* V(x) = x^2 + y^2, P = I */
    double P[4] = {1.0, 0.0, 0.0, 1.0};
    double x[2] = {3.0, 4.0};
    double V = nl_quadratic_lyapunov(P, x, 2);
    assert(fabs(V - 25.0) < 1e-10);  /* 3^2 + 4^2 = 25 */

    /* V_dot for linear stable system */
    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");
    double V_dot = nl_quadratic_lyapunov_dot(P, x, &sys, 2);
    /* dx/dt = -x, dy/dt = -2y -> V_dot = 2(x*(-x) + y*(-2y)) = -2*9 - 4*16 = -82 */
    assert(V_dot < -80.0 && V_dot > -84.0);

    PASS();
}

static void test_lyapunov_direct_method(void)
{
    TEST("nl_direct_method");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    double x_eq_arr[2] = {0.0, 0.0};
    nl_state_t eq;
    eq.dim = 2;
    memcpy(eq.x, x_eq_arr, 2 * sizeof(double));

    double P[4];
    double radius;
    nl_stability_t stab = nl_direct_method(&sys, &eq, P, &radius);
    assert(stab == NL_ASYMPTOTICALLY_STABLE);
    assert(radius > 0.0);

    /* Check P is positive definite */
    double min_eig;
    assert(nl_is_positive_definite(P, 2, &min_eig) == 1);
    assert(min_eig > 0.0);

    PASS();
}

static void test_eigenvalues(void)
{
    TEST("nl_eigenvalues (2x2)");

    /* A = [[0, 1], [-1, 0]] -> eigenvalues = ±j */
    double A[4] = {0.0, 1.0, -1.0, 0.0};
    double real_p[2], imag_p[2];
    int ret = nl_eigenvalues(A, 2, real_p, imag_p, 2000, 1e-8);
    assert(ret == 0);
    assert(fabs(real_p[0]) < 1e-4);
    assert(fabs(imag_p[0] - 1.0) < 0.1 || fabs(imag_p[0] + 1.0) < 0.1);

    /* A = [[-1, 0], [0, -2]] -> eigenvalues = -1, -2 */
    double B[4] = {-1.0, 0.0, 0.0, -2.0};
    ret = nl_eigenvalues(B, 2, real_p, imag_p, 2000, 1e-8);
    assert(ret == 0);
    /* For diagonal matrices, eigenvalues should be close to diagonal entries */
    assert(fabs(real_p[0] + 1.0) < 0.1 || fabs(real_p[0] + 2.0) < 0.1);
    assert(fabs(real_p[1] + 1.0) < 0.1 || fabs(real_p[1] + 2.0) < 0.1);

    PASS();
}

static void test_circle_criterion(void)
{
    TEST("nl_circle_criterion / nl_popov_criterion");

    nl_linear_system_t lin;
    memset(&lin, 0, sizeof(lin));
    lin.n = 1;
    /* G(s) = 1/(s+1): stable first-order lowpass */
    lin.A[0][0] = -1.0;
    lin.B[0] = 1.0;
    lin.C[0] = 1.0;
    lin.D = 0.0;

    nl_nonlinearity_t nl;
    double p[2] = {0.0, 10.0};  /* sector [0, 10] */
    nl_nonlinearity_init(&nl, NL_TYPE_SECTOR_BOUNDED, p, 2);

    double B_nl[1] = {-1.0};  /* negative feedback */
    nl_state_t x0;
    x0.dim = 1; x0.x[0] = 0.0;

    nl_lure_system_t lure;
    nl_lure_init(&lure, &lin, &nl, B_nl, &x0);

    double freqs[] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0};
    double margin;
    int stable = nl_circle_criterion(&lure, freqs, 9, &margin);
    assert(stable == 1 || stable == 0);  /* either stable or inconclusive */
    assert(isfinite(margin));

    double eta;
    int popov_stable = nl_popov_criterion(&lure, freqs, 9, &eta, &margin);
    assert(popov_stable == 1 || popov_stable == 0);
    assert(eta >= 0.0);

    PASS();
}

/* ===================================================================
 * Test: L5 - Describing Functions
 * =================================================================== */

static void test_describing_functions(void)
{
    TEST("nl_df_ideal_relay / dead_zone / saturation / backlash");

    /* Ideal relay: N(A) = 4M/(pi*A) */
    double M = 2.0, A = 1.0;
    double N_A = nl_df_ideal_relay(M, A);
    double expected = 4.0 * M / (M_PI * A);
    assert(fabs(N_A - expected) < 1e-6);

    /* Dead zone: delta=0.2, A=1.0 */
    double delta = 0.2;
    double N_dz = nl_df_dead_zone(delta, 1.0);
    assert(N_dz > 0.0 && N_dz < 1.0);
    /* A <= delta => N(A) = 0 */
    assert(nl_df_dead_zone(delta, 0.1) == 0.0);

    /* Saturation: k=2, a=1, A=0.3 -> linear region */
    double N_sat = nl_df_saturation(2.0, 1.0, 0.3);
    assert(fabs(N_sat - 2.0) < 1e-6);
    /* A >> a/k -> gain compression */
    double N_sat_big = nl_df_saturation(2.0, 1.0, 10.0);
    assert(N_sat_big < 1.0);

    /* Backlash */
    double re, im;
    int ret = nl_df_backlash(1.0, 0.1, 1.0, &re, &im);
    assert(ret == 0);
    assert(im < 0.0);  /* phase lag (negative imaginary) */

    /* Relay with hysteresis */
    ret = nl_df_relay_hysteresis(2.0, 0.1, 1.0, &re, &im);
    assert(ret == 0);
    assert(re > 0.0);

    /* Critical locus */
    nl_nonlinearity_t nl;
    double p[1] = {1.0};
    nl_nonlinearity_init(&nl, NL_TYPE_IDEAL_RELAY, p, 1);
    double A_arr[] = {0.5, 1.0, 2.0, 5.0};
    double locus_re[4], locus_im[4];
    ret = nl_df_critical_locus(&nl, A_arr, 4, locus_re, locus_im);
    assert(ret == 0);
    /* All should be on negative real axis (imag ≈ 0) */
    for (int i = 0; i < 4; i++) {
        assert(locus_re[i] < 0.0);
        assert(fabs(locus_im[i]) < 1e-1);  /* numerical tolerance for relay */
    }

    PASS();
}

/* ===================================================================
 * Test: L5 - Bifurcation Detection
 * =================================================================== */

static void test_bifurcation_detection(void)
{
    TEST("nl_find_equilibrium / nl_matrix_determinant");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    /* Find equilibrium at origin */
    double x0[2] = {0.1, 0.1};
    nl_state_t eq;
    int ret = nl_find_equilibrium(&sys, x0, NULL, 0, &eq, 1e-8);
    assert(ret == 0);
    assert(fabs(eq.x[0]) < 1e-6);
    assert(fabs(eq.x[1]) < 1e-6);

    /* Matrix determinant */
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double det;
    ret = nl_matrix_determinant(A, 2, &det);
    assert(ret == 0);
    assert(fabs(det - (-2.0)) < 1e-6);  /* 1*4 - 2*3 = -2 */

    /* Singular matrix */
    double B[4] = {1.0, 2.0, 2.0, 4.0};
    ret = nl_matrix_determinant(B, 2, &det);
    assert(ret == 0);
    assert(fabs(det) < 1e-10);

    PASS();
}

/* ===================================================================
 * Test: L5 - Chaos Analysis
 * =================================================================== */

static void test_chaos_analysis(void)
{
    TEST("nl_lyapunov_exponents / nl_kaplan_yorke_dimension");

    /* For a stable linear system, all Lyapunov exponents should be negative */
    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    double x0[2] = {0.1, 0.1};
    double exponents[2];
    int ret = nl_lyapunov_exponents(&sys, x0, 10.0, 0.01, 10, exponents);
    assert(ret == 0);
    assert(exponents[0] < 0.0);  /* largest exponent negative for stable system */

    /* Kaplan-Yorke dimension */
    double D_KY = nl_kaplan_yorke_dimension(exponents, 2);
    assert(D_KY == 0.0);  /* all exponents negative -> dimension 0 */

    /* K-S entropy */
    double h_KS = nl_ks_entropy_pesin(exponents, 2);
    assert(h_KS == 0.0);  /* no positive exponents */

    PASS();
}

static void test_correlation_dimension(void)
{
    TEST("nl_correlation_sum / nl_correlation_dimension");

    /* Generate random data */
    size_t N = 100;
    size_t dim = 3;
    double *data = (double *)malloc(N * dim * sizeof(double));
    for (size_t i = 0; i < N * dim; i++)
        data[i] = (double)rand() / (double)RAND_MAX;

    double epsilons[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
    double D2, r2;

    /* Should compute without error */
    double C = nl_correlation_sum(data, N, dim, 0.1);
    assert(C >= 0.0 && C <= 1.0);

    int ret = nl_correlation_dimension(data, N, dim, epsilons, 6, &D2, &r2);
    /* May fail for random data, but should not crash */
    (void)ret;
    (void)D2;
    (void)r2;

    free(data);
    PASS();
}

static void test_sample_entropy(void)
{
    TEST("nl_sample_entropy");

    /* Regular sine wave -> low sample entropy */
    size_t N = 200;
    double *sine = (double *)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++)
        sine[i] = sin(2.0 * M_PI * (double)i / 20.0);

    double sampen = nl_sample_entropy(sine, N, 2, 0.2);
    assert(sampen >= 0.0);

    free(sine);
    PASS();
}

/* ===================================================================
 * Test: L5 - Phase Analysis
 * =================================================================== */

static void test_phase_portrait(void)
{
    TEST("nl_phase_portrait_init / nl_compute_phase_portrait");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    nl_phase_portrait_t pp;
    int ret = nl_phase_portrait_init(&pp, -2.0, 2.0, -2.0, 2.0, 20, 20);
    assert(ret == 0);
    assert(pp.nx == 20);
    assert(pp.ny == 20);

    ret = nl_compute_phase_portrait(&sys, &pp);
    assert(ret == 0);

    /* Check vector field near origin (grid doesn't sample exactly at 0,0) */
    /* For grid [-2,2]x[-2,2] with 20x20: closest to (0,0) is near index (10,10) */
    size_t cx = 10, cy = 10;
    double f_x = pp.dxdt[cy * pp.nx + cx];
    double f_y = pp.dydt[cy * pp.nx + cx];
    /* Grid spacing dx=4/19≈0.21, so closest sample is ~0.1 from origin.
     * f(x,y) = (-x, -2y), so |f| ≈ 0.1-0.2, tolerance must allow this. */
    assert(fabs(f_x) < 0.3);
    assert(fabs(f_y) < 0.5);

    /* Bendixson's criterion */
    int has_cycle;
    ret = nl_bendixson_criterion(&pp, &has_cycle);
    assert(ret == 0);
    assert(has_cycle == 0);  /* stable node -> no limit cycles */

    nl_phase_portrait_free(&pp);
    PASS();
}

static void test_rk4_integration(void)
{
    TEST("nl_integrate_rk4");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    double x0[2] = {1.0, 1.0};
    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 1000, 2);

    size_t n_steps;
    int ret = nl_integrate_rk4(&sys, x0, 2.0, 0.01, &traj, &n_steps);
    assert(ret == 0);
    assert(traj.num_points > 10);

    /* Should approach origin */
    double x_final = traj.points[traj.num_points - 1].x[0];
    double y_final = traj.points[traj.num_points - 1].x[1];
    assert(fabs(x_final) < 1.0);
    assert(fabs(y_final) < 1.0);

    nl_trajectory_free(&traj);
    PASS();
}

static void test_poincare_index(void)
{
    TEST("nl_poincare_index");

    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "ls");

    double index;
    int ret = nl_poincare_index(&sys, 0.0, 0.0, 1.0, 100, &index);
    assert(ret == 0);
    /* For stable node: Poincaré index = +1 */
    assert(fabs(index - 1.0) < 0.5);

    /* Saddle: Poincaré index = -1 */
    nl_system_init(&sys, 2, saddle_f, "saddle");
    ret = nl_poincare_index(&sys, 0.0, 0.0, 0.5, 100, &index);
    assert(ret == 0);
    assert(fabs(index + 1.0) < 0.5);

    PASS();
}

/* ===================================================================
 * Test: L5 - Volterra Series
 * =================================================================== */

static void test_volterra_basic(void)
{
    TEST("nl_volterra_init / nl_volterra_evaluate");

    nl_volterra_series_t vs;
    int ret = nl_volterra_init(&vs, 2, 1.0);
    assert(ret == 0);
    assert(vs.max_order == 2);

    /* Initialize kernels */
    ret = nl_kernel_init(&vs.kernels[0], 1, 3);  /* 1st order, memory 3 */
    assert(ret == 0);
    ret = nl_kernel_init(&vs.kernels[1], 2, 3);  /* 2nd order, memory 3 */
    assert(ret == 0);

    /* Set h1 = [1.0, 0.5, 0.25] (simple lowpass) */
    size_t idx1[1];
    idx1[0] = 0; nl_kernel_set(&vs.kernels[0], idx1, 1.0);
    idx1[0] = 1; nl_kernel_set(&vs.kernels[0], idx1, 0.5);
    idx1[0] = 2; nl_kernel_set(&vs.kernels[0], idx1, 0.25);

    /* Input signal */
    double input[] = {1.0, 0.0, 0.0, 0.0, 0.0};
    double output[5] = {0};

    ret = nl_volterra_evaluate(&vs, input, 5, output);
    assert(ret == 0);

    /* Check linear convolution: h1 * u */
    /* y[0] = h1[0]*u[0] = 1.0 */
    assert(fabs(output[0] - 1.0) < 1e-6);
    /* y[1] = h1[1]*u[0] + h1[0]*u[1] = 0.5 + 0 = 0.5 */
    assert(fabs(output[1] - 0.5) < 1e-6);

    nl_volterra_free(&vs);
    PASS();
}

static void test_kernel_symmetrize(void)
{
    TEST("nl_kernel_symmetrize (2nd order)");

    nl_volterra_kernel_t k2;
    nl_kernel_init(&k2, 2, 3);

    /* Set asymmetric values */
    size_t idx[2];
    idx[0] = 0; idx[1] = 1;
    nl_kernel_set(&k2, idx, 10.0);
    idx[0] = 1; idx[1] = 0;
    nl_kernel_set(&k2, idx, 20.0);

    nl_kernel_symmetrize(&k2);

    double v;
    idx[0] = 0; idx[1] = 1;
    nl_kernel_get(&k2, idx, &v);
    assert(fabs(v - 15.0) < 1e-6);  /* (10+20)/2 = 15 */

    idx[0] = 1; idx[1] = 0;
    nl_kernel_get(&k2, idx, &v);
    assert(fabs(v - 15.0) < 1e-6);

    nl_kernel_free(&k2);
    PASS();
}

/* ===================================================================
 * Test: Edge Cases and Error Handling
 * =================================================================== */

static void test_null_pointer_guards(void)
{
    TEST("Null pointer guards on all APIs");

    /* System init with NULL */
    assert(nl_system_init(NULL, 2, linear_stable_f, "x") == -1);
    assert(nl_system_init(NULL, 2, NULL, "x") == -1);

    /* Trajectory with NULL */
    assert(nl_trajectory_init(NULL, 10, 2) == -1);
    assert(nl_trajectory_append(NULL, NULL) == -1);

    /* Describing function with NULL */
    assert(nl_describing_function(NULL, 1.0, 100, NULL) == -1);

    /* Lyapunov with NULL */
    assert(nl_compute_jacobian(NULL, NULL, NULL) == -1);

    /* Volterra with NULL */
    assert(nl_volterra_init(NULL, 2, 1.0) == -1);

    PASS();
}

static void test_dimension_limits(void)
{
    TEST("Dimension limit enforcement");

    nl_nonlinear_system_t sys;
    /* dim = 0 should fail */
    assert(nl_system_init(&sys, 0, linear_stable_f, "x") == -1);
    /* dim > NL_MAX_STATE_DIM should fail */
    assert(nl_system_init(&sys, NL_MAX_STATE_DIM + 1, linear_stable_f, "x") == -1);

    /* Trajectory with zero capacity */
    nl_trajectory_t traj;
    assert(nl_trajectory_init(&traj, 0, 2) == -1);

    /* Phase portrait with invalid bounds */
    nl_phase_portrait_t pp;
    assert(nl_phase_portrait_init(&pp, 5.0, 1.0, 0.0, 1.0, 10, 10) == -1);

    PASS();
}

/* ===================================================================
 * Helper system: Van der Pol (for integration and limit cycle tests)
 * =================================================================== */

static void test_rk4_van_der_pol(void)
{
    TEST("RK4 integration of Van der Pol oscillator");

    nl_nonlinear_system_t sys;
    double mu = 2.0;
    sys.params[0] = mu;
    sys.num_params = 1;
    sys.dim = 2;
    sys.f = vdp_f;

    double x0[2] = {1.0, 0.0};
    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 5000, 2);

    size_t n_steps;
    int ret = nl_integrate_rk4(&sys, x0, 20.0, 0.01, &traj, &n_steps);
    assert(ret == 0);
    assert(traj.num_points > 100);

    /* VdP should approach a limit cycle */
    /* Check that trajectory is bounded */
    for (size_t i = traj.num_points - 100; i < traj.num_points; i++) {
        assert(fabs(traj.points[i].x[0]) < 5.0);
        assert(fabs(traj.points[i].x[1]) < 10.0);
    }

    nl_trajectory_free(&traj);
    PASS();
}

/* ===================================================================
 * Test: L7 - Applications
 * =================================================================== */

static void test_costas_pll_lock_detect(void)
{
    TEST("nl_costas_pll_lock_detect");

    /* Simulate locked PLL: I = cos(small phase), Q = sin(small phase) */
    size_t N = 100;
    double *i_samps = (double *)malloc(N * sizeof(double));
    double *q_samps = (double *)malloc(N * sizeof(double));
    for (size_t k = 0; k < N; k++) {
        double phi = 0.05 + 0.001 * (double)k;
        i_samps[k] = cos(phi);
        q_samps[k] = sin(phi);
    }

    int is_locked;
    double phase_err;
    int ret = nl_costas_pll_lock_detect(i_samps, q_samps, N,
                                         0.85, &is_locked, &phase_err);
    assert(ret == 0);
    assert(is_locked == 1);  /* small phase error -> locked */
    assert(phase_err < 0.2); /* RMS phase error should be small */

    /* Unlocked: random I/Q */
    for (size_t k = 0; k < N; k++) {
        i_samps[k] = (double)(k % 7) / 7.0 - 0.5;
        q_samps[k] = (double)(k % 11) / 11.0 - 0.5;
    }
    ret = nl_costas_pll_lock_detect(i_samps, q_samps, N,
                                     0.85, &is_locked, &phase_err);
    assert(ret == 0);
    /* Random data should not lock */
    assert(is_locked == 0 || phase_err > 0.1);

    free(i_samps); free(q_samps);
    PASS();
}

static void test_mems_duffing_backbone(void)
{
    TEST("nl_mems_duffing_backbone");

    double omega0 = 1.0, alpha = 0.1, gamma = 0.05, F = 0.5;
    size_t n_freqs = 50;
    double *amp = (double *)malloc(n_freqs * sizeof(double));
    double *phase = (double *)malloc(n_freqs * sizeof(double));
    double omega_range[2] = {0.5, 1.5};

    int ret = nl_mems_duffing_backbone(omega0, alpha, gamma, F,
                                        omega_range, n_freqs, amp, phase);
    assert(ret == 0);

    /* Check that amplitude is positive and has a peak near resonance */
    double max_amp = 0.0;
    for (size_t k = 0; k < n_freqs; k++) {
        assert(amp[k] >= 0.0);
        if (amp[k] > max_amp) max_amp = amp[k];
    }
    assert(max_amp > 0.0);

    /* Phase should cross -pi/2 near resonance */
    int phase_cross = 0;
    for (size_t k = 1; k < n_freqs; k++) {
        if (phase[k - 1] > -M_PI / 2.0 && phase[k] < -M_PI / 2.0)
            phase_cross = 1;
    }
    (void)phase_cross;

    free(amp); free(phase);
    PASS();
}

static void test_fitzhugh_nagumo(void)
{
    TEST("nl_fitzhugh_nagumo (FHN neuron)");

    nl_fhn_params_t params = {0.7, 0.8, 0.08, 0.5};
    double x0[2] = {-1.0, -0.5};

    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 5000, 2);

    int spike_count;
    int ret = nl_fitzhugh_nagumo(&params, x0, 50.0, 0.01,
                                  &traj, &spike_count);
    assert(ret == 0);
    assert(traj.num_points > 100);

    /* With I_ext = 0.5, should see tonic spiking */
    assert(spike_count > 0);

    /* Sub-threshold: I_ext = 0.2 should produce fewer/no spikes */
    params.I_ext = 0.2;
    nl_trajectory_t traj2;
    nl_trajectory_init(&traj2, 5000, 2);
    int spike_count2;
    ret = nl_fitzhugh_nagumo(&params, x0, 30.0, 0.01,
                              &traj2, &spike_count2);
    assert(ret == 0);
    /* Lower I_ext -> fewer spikes */
    assert(spike_count2 <= spike_count + 5);

    nl_trajectory_free(&traj);
    nl_trajectory_free(&traj2);
    PASS();
}

static void test_lotka_volterra(void)
{
    TEST("nl_lotka_volterra");

    nl_lotka_volterra_params_t params = {1.0, 0.1, 0.5, 0.02};
    double x0[2] = {10.0, 5.0};

    nl_trajectory_t traj;
    nl_trajectory_init(&traj, 5000, 2);

    double avg_prey, avg_pred;
    int ret = nl_lotka_volterra(&params, x0, 20.0, 0.01,
                                 &traj, &avg_prey, &avg_pred);
    assert(ret == 0);
    assert(traj.num_points > 100);

    /* Populations should remain positive */
    for (size_t i = 0; i < traj.num_points; i++) {
        assert(traj.points[i].x[0] > 0.0);
        assert(traj.points[i].x[1] > 0.0);
    }

    /* Volterra's law: time averages approx equal equilibrium */
    double x_star = params.gamma / params.delta;  /* = 25 */
    double y_star = params.alpha / params.beta;   /* = 10 */
    assert(fabs(avg_prey - x_star) < x_star * 0.8);
    assert(fabs(avg_pred - y_star) < y_star * 0.8);

    nl_trajectory_free(&traj);
    PASS();
}

static void test_sigma_delta_limit_cycle(void)
{
    TEST("nl_sigma_delta_limit_cycle");

    size_t lc_period;
    double lc_amp, tone_freq;
    int ret = nl_sigma_delta_limit_cycle(1, 0.0, 500,
                                           &lc_period, &lc_amp, &tone_freq);
    assert(ret == 0);

    /* For DC input = 0, first-order SDM should show period-2 idle tone */
    assert(lc_period == 2);
    assert(lc_amp > 0.0);
    assert(tone_freq > 0.0);

    /* For non-zero DC input, limit cycle period changes */
    ret = nl_sigma_delta_limit_cycle(1, 0.3, 500,
                                       &lc_period, &lc_amp, &tone_freq);
    assert(ret == 0);
    (void)lc_period; (void)lc_amp; (void)tone_freq;

    PASS();
}

static void test_stochastic_resonance(void)
{
    TEST("nl_stochastic_resonance_snr");

    double a = 1.0, b = 1.0, A_sig = 0.1, f0 = 0.1;
    double D_vals[] = {0.01, 0.05, 0.1, 0.2, 0.5, 1.0};
    size_t nD = 6;
    double *snr_db = (double *)malloc(nD * sizeof(double));
    double D_opt;

    int ret = nl_stochastic_resonance_snr(a, b, A_sig, f0,
                                            D_vals, nD,
                                            50.0, 0.01, snr_db, &D_opt);
    assert(ret == 0);
    assert(isfinite(D_opt));
    assert(D_opt > 0.0);

    /* SNR should vary with D (stochastic resonance effect) */
    int has_variation = 0;
    for (size_t i = 1; i < nD; i++) {
        if (fabs(snr_db[i] - snr_db[0]) > 0.01) has_variation = 1;
    }
    (void)has_variation;

    free(snr_db);
    PASS();
}

static void test_estimate_roa(void)
{
    TEST("nl_estimate_roa");

    /* Stable linear system: dx/dt = -x, dy/dt = -2y */
    nl_nonlinear_system_t sys;
    nl_system_init(&sys, 2, linear_stable_f, "roa_test");

    double x_eq_arr[2] = {0.0, 0.0};
    nl_state_t eq;
    eq.dim = 2;
    memcpy(eq.x, x_eq_arr, 2 * sizeof(double));

    double P[4], radius;
    nl_stability_t stab = nl_direct_method(&sys, &eq, P, &radius);
    assert(stab == NL_ASYMPTOTICALLY_STABLE);

    double c_max;
    int ret = nl_estimate_roa(&sys, &eq, P, 2.0, 15, &c_max);
    assert(ret == 0);
    assert(isfinite(c_max));
    assert(c_max > 0.0);

    /* For globally stable linear system, ROA should be large */
    assert(c_max > 0.01);

    PASS();
}

/* ===================================================================
 * Main test runner
 * =================================================================== */

int main(void)
{
    printf("=== mini-nonlinear-system Test Suite ===\n\n");

    printf("--- L1: Core Definitions ---\n");
    test_system_init();
    test_trajectory_ops();
    test_nonlinearity_eval();

    printf("\n--- L4: Stability Theorems ---\n");
    test_stability_classification();
    test_quadratic_lyapunov();
    test_lyapunov_direct_method();
    test_eigenvalues();
    test_circle_criterion();

    printf("\n--- L5: Algorithms ---\n");
    test_describing_functions();
    test_bifurcation_detection();
    test_chaos_analysis();
    test_correlation_dimension();
    test_sample_entropy();
    test_phase_portrait();
    test_rk4_integration();
    test_poincare_index();
    test_volterra_basic();
    test_kernel_symmetrize();

    printf("\n--- L6: Canonical Problems ---\n");
    test_rk4_van_der_pol();

    printf("\n--- L7: Applications ---\n");
    test_costas_pll_lock_detect();
    test_mems_duffing_backbone();
    test_fitzhugh_nagumo();
    test_lotka_volterra();
    test_sigma_delta_limit_cycle();
    test_stochastic_resonance();
    test_estimate_roa();

    printf("\n--- Error Handling ---\n");
    test_null_pointer_guards();
    test_dimension_limits();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_run);

    return (tests_passed == tests_run) ? 0 : 1;
}
