/**
 * @file test_kalman.c
 * @brief Test suite for Kalman filter module
 *
 * Tests: allocation, initialization, predict, update, NIS,
 * log-likelihood, steady-state, EKF (nonlinear pendulum),
 * RTS smoother, information form, RMSE, NEES, GPS tracker.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/system_defs.h"
#include "../include/state_space.h"
#include "../include/kalman_filter.h"

#define TOL 1e-6

static int passed = 0, failed = 0;
#define T(n) printf("  TEST: %s ... ", n)
#define P() do { printf("PASS\n"); passed++; } while(0)
#define F(m) do { printf("FAIL: %s\n", m); failed++; } while(0)
#define C(c,m) do { if(!(c)){F(m);return;} } while(0)

static void test_kf_alloc_free(void)
{
    T("kf_alloc/free (n=4,m=2,l=3)");
    kalman_filter_t kf = kf_alloc(4, 2, 3);
    C(kf.x != NULL && kf.P != NULL, "allocation failed");
    C(kf.n == 4 && kf.m == 2 && kf.l == 3, "dimensions wrong");
    C(kf.owns_data == 1, "owns_data not set");
    kf_free(&kf);
    C(kf.x == NULL && kf.n == 0, "free did not reset");
    P();
}

static void test_kf_initialize(void)
{
    T("kf_initialize (1D CV model)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 1.0};
    double P0[4] = {10,0,0,10};
    double F[4]  = {1,1,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0.1,0,0,0.1};
    double R[1]  = {1.0};
    int ret = kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    C(ret == 0, "initialize failed");
    C(fabs(kf.x[0]) < TOL, "x[0] wrong");
    C(fabs(kf.x[1] - 1.0) < TOL, "x[1] wrong");
    kf_free(&kf);
    P();
}

static void test_kf_predict_update(void)
{
    T("kf_predict + kf_update (reduces uncertainty)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 0.0};
    double P0[4] = {100,0,0,100};
    double F[4]  = {1,1,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0.01,0.001,0.001,0.01};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    C(kf_predict(&kf, NULL) == 0, "predict failed");
    double z[1] = {5.0};
    C(kf_update(&kf, z) == 0, "update failed");
    C(fabs(kf.x[0] - 5.0) < 5.0, "no correction toward measurement");
    double trace_before = 200.0;
    double trace_after = kf.P[0] + kf.P[3];
    C(trace_after < trace_before, "uncertainty did not decrease");
    kf_free(&kf);
    P();
}

static void test_kf_step(void)
{
    T("kf_step (combined predict+update)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 0.0};
    double P0[4] = {100,0,0,100};
    double F[4]  = {1,1,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0.01,0.001,0.001,0.01};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    double z[1] = {10.0};
    C(kf_step(&kf, NULL, z) == 0, "step failed");
    C(kf.x[0] > 0.0, "position did not move toward measurement");
    kf_free(&kf);
    P();
}

static void test_kf_innovation(void)
{
    T("kf_innovation (zero innovation at measurement)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {5.0, 1.0};
    double P0[4] = {1,0,0,1};
    double F[4]  = {1,0,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0,0,0,0};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    double z[1] = {5.0}, y[1];
    C(kf_innovation(&kf, z, y) == 0, "innovation failed");
    C(fabs(y[0]) < TOL, "innovation not zero when z=Hx");
    z[0] = 7.0;
    kf_innovation(&kf, z, y);
    C(fabs(y[0] - 2.0) < TOL, "innovation wrong");
    kf_free(&kf);
    P();
}

static void test_kf_nis(void)
{
    T("kf_nis (outlier detection)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 0.0};
    double P0[4] = {1,0,0,1};
    double F[4]  = {1,0,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0,0,0,0};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    double z[1] = {0.0};
    double nis = kf_nis(&kf, z);
    C(nis >= 0.0, "NIS negative");
    z[0] = 100.0;
    double nis2 = kf_nis(&kf, z);
    C(nis2 > nis, "NIS not larger for outlier");
    kf_free(&kf);
    P();
}

static void test_kf_log_likelihood(void)
{
    T("kf_log_likelihood (nominal > outlier)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 0.0};
    double P0[4] = {1,0,0,1};
    double F[4]  = {1,0,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0,0,0,0};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    double z_nom[1] = {0.0}, z_out[1] = {10.0};
    double ll_nom = kf_log_likelihood(&kf, z_nom);
    double ll_out = kf_log_likelihood(&kf, z_out);
    C(ll_nom > ll_out, "outlier more likely than nominal");
    C(isfinite(ll_nom), "log-likelihood not finite");
    kf_free(&kf);
    P();
}

static void test_kf_steady_state(void)
{
    T("kf_steady_state (DARE convergence)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2] = {0.0, 0.0};
    double P0[4] = {100,0,0,100};
    double F[4]  = {1,1,0,1};
    double H[2]  = {1,0};
    double Q[4]  = {0.1,0.05,0.05,0.1};
    double R[1]  = {1.0};
    kf_initialize(&kf, x0, P0, F, H, Q, R, NULL);
    int ret = kf_steady_state(&kf, 1e-6, 500);
    C(ret == 0, "DARE did not converge");
    C(kf.P[0] > 0.0 && kf.P[3] > 0.0, "P diagonal not positive");
    kf_free(&kf);
    P();
}

static void test_kf_rmse(void)
{
    T("kf_rmse (zero error -> zero RMSE)");
    double x_est[6]  = {1,2,3,4,5,6};
    double x_true[6] = {1,2,3,4,5,6};
    double err = kf_rmse(x_est, x_true, 2, 3);
    C(fabs(err) < TOL, "RMSE not zero for identical signals");
    x_est[0] = 3.0;
    err = kf_rmse(x_est, x_true, 2, 3);
    C(fabs(err - sqrt(4.0/6.0)) < TOL, "RMSE wrong");
    P();
}

static void test_kf_nees(void)
{
    T("kf_nees (consistent filter check)");
    int n = 2, M = 3;
    double x_est[6] = {0,0, 0.1,0.1, 0,0};
    double x_true[6] = {0,0, 0,0, 0,0};
    double P[12] = {1,0,0,1, 1,0,0,1, 1,0,0,1};
    double nees[3];
    int ret = kf_nees(x_est, x_true, P, M, n, nees);
    C(ret == 0, "NEES failed");
    C(fabs(nees[0]) < TOL, "NEES[0] not zero for zero error");
    C(fabs(nees[1] - 0.02) < TOL, "NEES[1] wrong");
    P();
}

static void test_kf_gps_tracker(void)
{
    T("kf_configure_gps_tracker (8-state GPS KF)");
    kalman_filter_t kf = kf_alloc(8, 0, 4);
    C(kf.x != NULL, "alloc failed");
    int ret = kf_configure_gps_tracker(&kf, 1.0, 5.0, 0.1, 1.0);
    C(ret == 0, "GPS config failed");
    C(fabs(kf.F[0*8+3] - 1.0) < TOL, "F coupling wrong");
    double q = 1.0;
    C(fabs(kf.Q[0] - q*q/3.0) < TOL, "Q[0] wrong for dt=1");
    double z[4] = {1000.0, 2000.0, 3000.0, 0.0};
    ret = kf_gps_update(&kf, z, 4);
    C(ret == 0, "GPS update failed");
    kf_free(&kf);
    P();
}

static void test_kf_cholesky_2x2(void)
{
    T("kf_cholesky_2x2");
    double A[4] = {4.0, 1.0, 1.0, 3.0};
    double L[4];
    int ret = kf_cholesky_2x2(A, L);
    C(ret == 0, "Cholesky failed on PD matrix");
    C(fabs(L[0] - 2.0) < TOL, "L[0] wrong");
    C(fabs(L[2] - 0.5) < TOL, "L[2] wrong");
    double recon[4];
    recon[0] = L[0]*L[0];
    recon[1] = L[0]*L[2];
    recon[2] = L[2]*L[0];
    recon[3] = L[2]*L[2] + L[3]*L[3];
    C(fabs(recon[0]-A[0])<TOL && fabs(recon[1]-A[1])<TOL, "L*L^T != A");
    P();
}

static void test_kf_process_noise(void)
{
    T("kf_process_noise_discrete (DWNA model)");
    double Q[4];
    int ret = kf_process_noise_discrete(0.1, 2.0, 1, Q);
    C(ret == 0, "DWNA failed");
    double dt=0.1, q=2.0, dt2=0.01, dt3=0.001;
    C(fabs(Q[0] - q*dt3/3.0) < TOL, "Q[0] wrong");
    C(fabs(Q[1] - q*dt2/2.0) < TOL, "Q[1] wrong");
    C(fabs(Q[3] - q*dt) < TOL, "Q[3] wrong");
    P();
}

/* EKF test: nonlinear pendulum model */
static void pendulum_f(const double *x, const double *u, double *nx, void *d)
{ (void)u;(void)d; double gL=9.81/1.0,b=0.1,dt=0.01;
  nx[0]=x[0]+dt*x[1]; nx[1]=x[1]+dt*(-gL*sin(x[0])-b*x[1]); }

static void pendulum_jac_f(const double *x, const double *u, double *F, void *d)
{ (void)u;(void)d; double gL=9.81/1.0,b=0.1,dt=0.01;
  F[0]=1.0;F[1]=dt;F[2]=-dt*gL*cos(x[0]);F[3]=1.0-dt*b; }

static void ph(const double *x, double *zp, void *d)
{ (void)d; zp[0]=x[0]; }

static void pjh(const double *x, double *H, void *d)
{ (void)x;(void)d; H[0]=1.0;H[1]=0.0; }

static void test_ekf_pendulum(void)
{
    T("ekf_predict + ekf_update (nonlinear pendulum)");
    kalman_filter_t kf = kf_alloc(2, 0, 1);
    double x0[2]={0.5,0.0}, P0[4]={0.1,0,0,0.1};
    double F0[4]={1,0.01,-0.0981*cos(0.5),0.999};
    double H[2]={1,0}, Q[4]={0.001,0,0,0.001}, R[1]={0.01};
    kf_initialize(&kf, x0, P0, F0, H, Q, R, NULL);
    int ret = ekf_predict(&kf, NULL, pendulum_f, pendulum_jac_f, NULL);
    C(ret == 0, "EKF predict failed");
    double z[1] = {0.3};
    ret = ekf_update(&kf, z, ph, pjh, NULL);
    C(ret == 0, "EKF update failed");
    C(fabs(kf.x[0]-0.3) < fabs(0.5-0.3)+0.1, "EKF did not correct");
    kf_free(&kf);
    P();
}

static void test_rts_smooth(void)
{
    T("kf_rts_smooth (smoothing improves estimate)");
    int n=2,N=3;
    double xf[6]={0,0,1,1,2,2},Pf[12]={1,0,0,1,1,0,0,1,1,0,0,1};
    double xp[6]={0,0,1,1,2,2},Pp[12]={2,0,0,2,2,0,0,2,2,0,0,2};
    double Fa[8]={1,0,0,1,1,0,0,1};
    double xs[6]={0},Ps[12]={0};
    int ret=kf_rts_smooth(xf,Pf,xp,Pp,Fa,N,n,xs,Ps);
    C(ret==0,"RTS smoother failed");
    C(fabs(xs[4]-2.0)<TOL,"smoothed last wrong");
    P();
}

static void test_info_form(void)
{
    T("kf_to_information_form (Y*P=I)");
    kalman_filter_t kf=kf_alloc(2,0,1);
    double x0[2]={3.0,1.0},P0[4]={4,1,1,3};
    double F[4]={1,0,0,1},H[2]={1,0},Q[4]={0,0,0,0},R[1]={1.0};
    kf_initialize(&kf,x0,P0,F,H,Q,R,NULL);
    double Y[4],y[2];
    C(kf_to_information_form(&kf,Y,y)==0,"info form failed");
    double YP0=Y[0]*P0[0]+Y[1]*P0[2], YP1=Y[0]*P0[1]+Y[1]*P0[3];
    C(fabs(YP0-1.0)<1e-4 && fabs(YP1)<1e-4,"Y*P != I");
    kf_free(&kf);
    P();
}

int main(void)
{
    printf("=== mini-system-analysis Kalman Filter Test Suite ===\n\n");
    test_kf_alloc_free();
    test_kf_initialize();
    test_kf_predict_update();
    test_kf_step();
    test_kf_innovation();
    test_kf_nis();
    test_kf_log_likelihood();
    test_kf_steady_state();
    test_kf_rmse();
    test_kf_nees();
    test_kf_gps_tracker();
    test_kf_cholesky_2x2();
    test_ekf_pendulum();
    test_rts_smooth();
    test_info_form();
    test_kf_process_noise();
    printf("\n=== Results: %d passed, %d failed ===\n", passed, failed);
    return (failed > 0) ? 1 : 0;
}
