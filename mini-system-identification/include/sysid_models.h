/**
 * @file sysid_models.h
 * @brief Parametric model structures for system identification (L2 Concepts, L3 Math Structures)
 *
 * Knowledge coverage:
 *   L2 - Parametric model concept: ARX, ARMAX, OE, BJ, State-Space
 *   L2 - Predictor forms: one-step-ahead predictor ŷ(k|k-1;θ)
 *   L2 - Model order concept and model set hierarchy
 *   L3 - Polynomial representations: A(q), B(q), C(q), D(q), F(q)
 *   L3 - Shift operator q: q⁻¹ y(k) = y(k-1)
 *   L3 - State-space matrix algebra (A, B, C, D, K)
 *   L3 - Transfer function ↔ state-space conversion
 *   L6 - Canonical forms: controllable, observable, diagonal
 *
 * The general linear model (Ljung 1999, Eq 4.15):
 *   A(q) y(k) = [B(q)/F(q)] u(k-nk) + [C(q)/D(q)] e(k)
 *   A(q) = 1 + a₁q⁻¹ + ... + aₙₐq⁻ⁿᵃ
 *   B(q) = b₁q⁻¹ + ... + bₙₐq⁻ⁿᵇ (note: B starts at q⁻¹, delay handled by nk)
 *
 * References:
 *   Ljung (1999) "System Identification: Theory for the User", 2nd ed.
 *   Söderström & Stoica (1989) "System Identification"
 *   Box, Jenkins, Reinsel (2008) "Time Series Analysis", 4th ed.
 */

#ifndef SYSID_MODELS_H
#define SYSID_MODELS_H

#include "sysid_types.h"

/*---------------------------------------------------------------------------
 * L2/L3: Polynomial model structures (input-output form)
 * Each structure stores the polynomial coefficients explicitly.
 *---------------------------------------------------------------------------*/

/** ARX model: A(q) y(k) = B(q) u(k-nk) + e(k)
 *  Predictor: ŷ(k|θ) = B(q) u(k-nk) + [1 - A(q)] y(k)
 */
typedef struct {
    sysid_model_order_t order;
    double *A;          /**< A polynomial [na+1], A[0]=1 */
    double *B;          /**< B polynomial [nb] */
    int     nk;         /**< Input delay in samples */
    double  Ts;         /**< Sampling period */
} sysid_arx_t;

/** ARMAX model: A(q) y(k) = B(q) u(k-nk) + C(q) e(k)
 *  Predictor: C(q) ŷ(k|θ) = B(q) u(k-nk) + [C(q) - A(q)] y(k)
 */
typedef struct {
    sysid_model_order_t order;
    double *A;          /**< A polynomial [na+1], A[0]=1 */
    double *B;          /**< B polynomial [nb] */
    double *C;          /**< C polynomial [nc+1], C[0]=1 */
    int     nk;         /**< Input delay */
    double  Ts;
} sysid_armax_t;

/** Output Error (OE) model: y(k) = [B(q)/F(q)] u(k-nk) + e(k)
 *  Predictor: ŷ(k|θ) = [B(q)/F(q)] u(k-nk)  (simulates without noise model)
 */
typedef struct {
    sysid_model_order_t order;
    double *B;          /**< B polynomial [nb] */
    double *F;          /**< F polynomial [nf+1], F[0]=1 */
    int     nk;
    double  Ts;
} sysid_oe_t;

/** Box-Jenkins (BJ) model: y(k) = [B(q)/F(q)] u(k-nk) + [C(q)/D(q)] e(k)
 *  Most general linear input-output model
 */
typedef struct {
    sysid_model_order_t order;
    double *B;          /**< B polynomial [nb] */
    double *C;          /**< C polynomial [nc+1], C[0]=1 */
    double *D;          /**< D polynomial [nd+1], D[0]=1 */
    double *F;          /**< F polynomial [nf+1], F[0]=1 */
    int     nk;
    double  Ts;
} sysid_bj_t;

/** FIR model: y(k) = B(q) u(k-nk) + e(k)
 *  Special case: non-parametric impulse response estimate
 */
typedef struct {
    int     nb;         /**< FIR length (number of taps) */
    double *B;          /**< Impulse response coefficients [nb] */
    int     nk;
    double  Ts;
} sysid_fir_t;

/** Continuous-time transfer function:
 *  G(s) = K * (bₘsᵐ + ... + b₁s + b₀) / (aₙsⁿ + ... + a₁s + a₀) * exp(-τs)
 *  Used for physical parameter estimation
 */
typedef struct {
    int     num_order;     /**< Numerator order m */
    int     den_order;     /**< Denominator order n */
    double *num;           /**< Numerator coefficients [m+1], descending powers */
    double *den;           /**< Denominator coefficients [n+1], aₙ..a₀ */
    double  K;             /**< Overall gain */
    double  delay;         /**< Transport delay τ [seconds] */
} sysid_tf_ct_t;

/*---------------------------------------------------------------------------
 * L3: State-space model structure
 *   x(k+1) = A x(k) + B u(k) + K e(k)
 *   y(k)   = C x(k) + D u(k) + e(k)
 *
 * Where K is the Kalman gain (innovation form).
 * For output-error models, K = 0.
 *---------------------------------------------------------------------------*/

typedef struct {
    int      nx;           /**< State dimension */
    int      nu;           /**< Input dimension */
    int      ny;           /**< Output dimension */
    double  *A;            /**< State matrix [nx × nx], row-major */
    double  *B;            /**< Input matrix [nx × nu], row-major */
    double  *C;            /**< Output matrix [ny × nx], row-major */
    double  *D;            /**< Feedthrough matrix [ny × nu], row-major */
    double  *K;            /**< Kalman gain [nx × ny], row-major; NULL = no noise model */
    double  *x0;           /**< Initial state estimate [nx] */
    double   Ts;           /**< Sampling period */
    int      is_discrete;  /**< 1 if discrete-time, 0 if continuous-time */
} sysid_ss_t;

/*---------------------------------------------------------------------------
 * L6: Canonical forms for state-space models
 *---------------------------------------------------------------------------*/

/** Canonical form enumeration (Kailath 1980, Ch. 2.4) */
typedef enum {
    SYSID_SS_CANON_CONTROLLABLE = 0,  /**< Controllable canonical form */
    SYSID_SS_CANON_OBSERVABLE   = 1,  /**< Observable canonical form */
    SYSID_SS_CANON_MODAL        = 2,  /**< Modal (diagonal/Jordan) form */
    SYSID_SS_CANON_BALANCED     = 3   /**< Balanced realization */
} sysid_ss_canon_t;

/*---------------------------------------------------------------------------
 * Model allocation and construction API
 *---------------------------------------------------------------------------*/

/** Allocate ARX model with given orders */
int  sysid_arx_alloc(sysid_arx_t *model, int na, int nb, int nk, double Ts);
void sysid_arx_free(sysid_arx_t *model);
/** Set ARX parameters from a flat parameter vector */
void sysid_arx_set_params(sysid_arx_t *model, const sysid_param_vector_t *pv);
/** Extract parameters from ARX model into a flat vector */
void sysid_arx_get_params(const sysid_arx_t *model, sysid_param_vector_t *pv);

int  sysid_armax_alloc(sysid_armax_t *model, int na, int nb, int nc, int nk, double Ts);
void sysid_armax_free(sysid_armax_t *model);
void sysid_armax_set_params(sysid_armax_t *model, const sysid_param_vector_t *pv);
void sysid_armax_get_params(const sysid_armax_t *model, sysid_param_vector_t *pv);

int  sysid_oe_alloc(sysid_oe_t *model, int nb, int nf, int nk, double Ts);
void sysid_oe_free(sysid_oe_t *model);
void sysid_oe_set_params(sysid_oe_t *model, const sysid_param_vector_t *pv);
void sysid_oe_get_params(const sysid_oe_t *model, sysid_param_vector_t *pv);

int  sysid_bj_alloc(sysid_bj_t *model, int nb, int nc, int nd, int nf, int nk, double Ts);
void sysid_bj_free(sysid_bj_t *model);
void sysid_bj_set_params(sysid_bj_t *model, const sysid_param_vector_t *pv);
void sysid_bj_get_params(const sysid_bj_t *model, sysid_param_vector_t *pv);

int  sysid_fir_alloc(sysid_fir_t *model, int nb, int nk, double Ts);
void sysid_fir_free(sysid_fir_t *model);

int  sysid_tf_ct_alloc(sysid_tf_ct_t *model, int num_order, int den_order);
void sysid_tf_ct_free(sysid_tf_ct_t *model);

int  sysid_ss_alloc(sysid_ss_t *model, int nx, int nu, int ny, int is_discrete, double Ts);
void sysid_ss_free(sysid_ss_t *model);

/*---------------------------------------------------------------------------
 * L3: Model simulation and prediction
 *---------------------------------------------------------------------------*/

/** One-step-ahead prediction ŷ(k|k-1; θ)
 *  @param model   Model structure
 *  @param u       Input data array [N × nu], row-major
 *  @param y       Output data array [N × ny], row-major (observed)
 *  @param ypred   Prediction output [N × ny], row-major (computed)
 *  @param N       Number of samples
 *  @return        0 on success, -1 on error
 */
int  sysid_arx_predict(const sysid_arx_t *model, const double *u, const double *y,
                       double *ypred, int N);
int  sysid_armax_predict(const sysid_armax_t *model, const double *u, const double *y,
                         double *ypred, int N);
int  sysid_oe_predict(const sysid_oe_t *model, const double *u, const double *y,
                      double *ypred, int N);
int  sysid_bj_predict(const sysid_bj_t *model, const double *u, const double *y,
                      double *ypred, int N);
int  sysid_fir_predict(const sysid_fir_t *model, const double *u, double *ypred, int N);

/** Pure simulation (k-step-ahead with zero initial conditions)
 *  ŷ(k) = G(q,θ) u(k) — uses past simulated outputs, not measured
 */
int  sysid_arx_simulate(const sysid_arx_t *model, const double *u, double *ysim, int N);
int  sysid_oe_simulate(const sysid_oe_t *model, const double *u, double *ysim, int N);
int  sysid_bj_simulate(const sysid_bj_t *model, const double *u, double *ysim, int N);

/** State-space one-step-ahead prediction (Kalman filter innov. form) */
int  sysid_ss_predict(const sysid_ss_t *model, const double *u, const double *y,
                      double *ypred, int N);
/** State-space k-step simulation (open-loop, no correction) */
int  sysid_ss_simulate(const sysid_ss_t *model, const double *u, double *ysim, int N);

/** Compute model poles (roots of A(q) or eigenvalues of A matrix) */
int  sysid_arx_poles(const sysid_arx_t *model, double complex *poles, int max_poles);
int  sysid_armax_poles(const sysid_armax_t *model, double complex *poles, int max_poles);
int  sysid_ss_poles(const sysid_ss_t *model, double complex *poles, int max_poles);

/** Compute DC gain (steady-state gain) of a model */
double sysid_arx_dcgain(const sysid_arx_t *model);
double sysid_oe_dcgain(const sysid_oe_t *model);
double sysid_ss_dcgain(const sysid_ss_t *model);

/** Compute impulse response of a model (first N samples) */
int  sysid_arx_impulse(const sysid_arx_t *model, double *impulse, int N);
int  sysid_oe_impulse(const sysid_oe_t *model, double *impulse, int N);
int  sysid_ss_impulse(const sysid_ss_t *model, double *impulse, int N, int input_ch);

/** Compute step response of a model (first N samples) */
int  sysid_arx_step(const sysid_arx_t *model, double *step_resp, int N);

/** Compute frequency response G(e^{jωTs}) at given frequencies */
int  sysid_arx_freqresp(const sysid_arx_t *model, const double *omega, int Nf,
                        double complex *G);
int  sysid_oe_freqresp(const sysid_oe_t *model, const double *omega, int Nf,
                       double complex *G);
int  sysid_ss_freqresp(const sysid_ss_t *model, const double *omega, int Nf,
                       double complex *G);

/*---------------------------------------------------------------------------
 * L3: Model conversion (transfer function ↔ state-space)
 *---------------------------------------------------------------------------*/

/** Convert ARX to state-space (observable canonical form) */
int  sysid_arx_to_ss(const sysid_arx_t *arx, sysid_ss_t *ss);

/** Convert state-space to ARX (for SISO systems via observer canonical) */
int  sysid_ss_to_arx(const sysid_ss_t *ss, sysid_arx_t *arx);

/** Convert continuous-time TF to discrete-time (Tustin/bilinear transform) */
int  sysid_tf_ct_to_dt(const sysid_tf_ct_t *ct, sysid_arx_t *dt, double Ts);

/** Convert state-space to canonical form */
int  sysid_ss_canonical(const sysid_ss_t *src, sysid_ss_t *dst, sysid_ss_canon_t form);

#endif /* SYSID_MODELS_H */
