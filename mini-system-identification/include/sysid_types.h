/**
 * @file sysid_types.h
 * @brief Core type definitions for system identification (L1 Definitions)
 *
 * Implements knowledge points:
 *   L1 - System ID problem formulation (Ljung, 1999; Söderström & Stoica, 1989)
 *   L1 - Model parameterization types
 *   L1 - Excitation signal classification
 *   L1 - Estimation bias-variance decomposition (Cramér-Rao lower bound)
 *   L1 - Identification experiment design parameters
 *
 * Reference texts:
 *   Ljung, L. "System Identification: Theory for the User" (2nd ed., 1999)
 *   Söderström, T. & Stoica, P. "System Identification" (1989)
 */

#ifndef SYSID_TYPES_H
#define SYSID_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

/* Maximum model order supported by default buffers */
#define SYSID_MAX_ORDER      128
#define SYSID_MAX_DATALEN    65536
#define SYSID_MAX_OUTPUTS    16
#define SYSID_MAX_INPUTS     16

/*---------------------------------------------------------------------------
 * L1: System identification problem types
 *     - Open-loop vs closed-loop
 *     - Time-domain vs frequency-domain
 *     - Continuous-time vs discrete-time
 *---------------------------------------------------------------------------*/

/** Experiment mode classification */
typedef enum {
    SYSID_OPEN_LOOP           = 0,  /**< System without feedback */
    SYSID_CLOSED_LOOP_DIRECT  = 1,  /**< Direct closed-loop identification */
    SYSID_CLOSED_LOOP_INDIRECT = 2,  /**< Indirect (two-stage) method */
    SYSID_JOINT_IO            = 3   /**< Joint input-output approach */
} sysid_experiment_mode_t;

/** Domain of identification data */
typedef enum {
    SYSID_TIME_DOMAIN     = 0,  /**< Sampled time-domain data */
    SYSID_FREQ_DOMAIN     = 1,  /**< Frequency response data */
    SYSID_CORR_DOMAIN     = 2   /**< Correlation function data */
} sysid_domain_t;

/*---------------------------------------------------------------------------
 * L1: Model structure taxonomy (Ljung 1999, Ch. 4)
 *     - Parametric vs non-parametric
 *     - Linear vs nonlinear
 *     - Input-output vs state-space
 *---------------------------------------------------------------------------*/

/** Linear parametric model structures (Ljung 1999, §4.2) */
typedef enum {
    SYSID_MODEL_FIR       = 0,   /**< Finite Impulse Response: B(q) */
    SYSID_MODEL_ARX       = 1,   /**< AutoRegressive with eXogenous: A(q)y = B(q)u + e */
    SYSID_MODEL_ARMAX     = 2,   /**< ARMA+eXogenous: A(q)y = B(q)u + C(q)e */
    SYSID_MODEL_OE        = 3,   /**< Output Error: y = B(q)/F(q) u + e */
    SYSID_MODEL_BJ        = 4,   /**< Box-Jenkins: y = B(q)/F(q) u + C(q)/D(q) e */
    SYSID_MODEL_SS        = 5,   /**< State-Space: x(k+1)=Ax(k)+Bu(k), y(k)=Cx(k)+Du(k) */
    SYSID_MODEL_TRANSFER  = 6,   /**< Continuous-time transfer function */
    SYSID_MODEL_NARX      = 7,   /**< Nonlinear ARX (L8) */
    SYSID_MODEL_HAMMERSTEIN = 8, /**< Hammerstein model (L8) */
    SYSID_MODEL_WIENER    = 9    /**< Wiener model (L8) */
} sysid_model_type_t;

/** Noise model structures */
typedef enum {
    SYSID_NOISE_WHITE      = 0,  /**< White noise e(k) ~ N(0, sigma^2) */
    SYSID_NOISE_MA         = 1,  /**< Moving Average: C(q)e(k) */
    SYSID_NOISE_ARMA       = 2,  /**< ARMA: C(q)/D(q) e(k) */
    SYSID_NOISE_COLORED    = 3   /**< Generic colored noise */
} sysid_noise_model_t;

/*---------------------------------------------------------------------------
 * L1: Parameterization and degrees of freedom
 *---------------------------------------------------------------------------*/

/** Model order specification */
typedef struct {
    int na;      /**< Order of A(q) polynomial (autoregressive part) */
    int nb;      /**< Order of B(q) polynomial (exogenous input part) */
    int nc;      /**< Order of C(q) polynomial (MA noise part) */
    int nd;      /**< Order of D(q) polynomial (AR noise part) */
    int nf;      /**< Order of F(q) polynomial (output error dynamics) */
    int nk;      /**< Input delay in samples */
    int nx;      /**< State dimension (for state-space models) */
    int ny;      /**< Number of outputs */
    int nu;      /**< Number of inputs */
} sysid_model_order_t;

/** Parameter vector - core representation of estimated model coefficients */
typedef struct {
    double *theta;         /**< Parameter vector (length = nparam) */
    int     nparam;        /**< Total number of parameters */
    int     nparam_free;   /**< Number of free (non-fixed) parameters */
    int    *free_idx;      /**< Indices of free parameters, NULL if all free */
} sysid_param_vector_t;

/** Gradient information for estimation algorithms */
typedef struct {
    double *psi;           /**< Gradient vector psi(k) = dŷ(k|θ)/dθ (length nparam) */
    double *hessian;       /**< Approximate Hessian (nparam × nparam, row-major) */
    int     nparam;
} sysid_gradient_t;

/*---------------------------------------------------------------------------
 * L1: Data structures for identification experiments
 *---------------------------------------------------------------------------*/

/** Time-domain input-output data record */
typedef struct {
    double *u;             /**< Input signal array [N × nu], row-major */
    double *y;             /**< Output signal array [N × ny], row-major */
    int     N;             /**< Number of samples */
    int     nu;            /**< Number of input channels */
    int     ny;            /**< Number of output channels */
    double  Ts;            /**< Sampling period [seconds] */
    int     is_real;       /**< 1 if real data, 0 if complex */
} sysid_data_t;

/** Frequency-domain data (for frequency response fitting) */
typedef struct {
    double   *omega;       /**< Frequency vector [rad/s] (length Nf) */
    double complex *G;     /**< Complex frequency response G(jω) (length Nf) */
    double   *weight;      /**< Frequency weighting (length Nf, can be NULL) */
    int       Nf;          /**< Number of frequency points */
} sysid_freq_data_t;

/** Identification experiment specification */
typedef struct {
    sysid_domain_t       domain;    /**< Time or frequency domain */
    sysid_experiment_mode_t mode;   /**< Open-loop or closed-loop */
    double               Ts;        /**< Sampling period */
    int                  N;         /**< Number of samples to collect */
    int                  N_est;     /**< Estimation data length */
    int                  N_val;     /**< Validation data length */
    double               input_amplitude; /**< Input signal amplitude */
    int                  input_type;      /**< Signal type (PRBS=0, chirp=1, etc.) */
} sysid_experiment_t;

/*---------------------------------------------------------------------------
 * L1: Validation metrics (bias-variance decomposition)
 *---------------------------------------------------------------------------*/

/** Single model fit result */
typedef struct {
    double  NRMSE_fit;       /**< Normalized RMSE fit [%]: 100*(1 - ||y-ŷ||/||y-mean(y)||) */
    double  RMSE;            /**< Root Mean Square Error */
    double  MAE;             /**< Mean Absolute Error */
    double  R_squared;       /**< Coefficient of determination R² */
    double  AIC;             /**< Akaike Information Criterion */
    double  AICc;            /**< Corrected AIC (for small samples) */
    double  BIC;             /**< Bayesian (Schwarz) Information Criterion */
    double  MDL;             /**< Minimum Description Length */
    double  FPE;             /**< Akaike's Final Prediction Error */
    double  residual_mean;   /**< Mean of residuals */
    double  residual_var;    /**< Variance of residuals */
    double  residual_whiteness; /**< Ljung-Box Q-statistic p-value approximation */
    int     residual_indep_test; /**< Cross-correlation independence test passed? (1=yes) */
    int     nparam;           /**< Number of estimated parameters */
    int     N_data;           /**< Number of data points used */
} sysid_fit_result_t;

/** Stability analysis result */
typedef struct {
    int     is_stable;       /**< 1 if all poles inside unit circle */
    int     n_poles;         /**< Number of poles */
    double complex *poles;   /**< Pole locations */
    double  max_pole_mag;    /**< Maximum pole magnitude */
    double  stability_margin; /**< Distance to unit circle (1 - max|pole|) */
} sysid_stability_t;

/*---------------------------------------------------------------------------
 * Core API: memory management and initialization
 *---------------------------------------------------------------------------*/

/** Initialize a model order with default zeros */
void sysid_order_init(sysid_model_order_t *order);

/** Check if model order is valid (non-negative, within bounds) */
int  sysid_order_is_valid(const sysid_model_order_t *order);

/** Compute total number of parameters for a given model type and order */
int  sysid_order_nparam(sysid_model_type_t type, const sysid_model_order_t *order);

/** Allocate and zero a parameter vector */
int  sysid_param_alloc(sysid_param_vector_t *pv, int nparam);

/** Free a parameter vector */
void sysid_param_free(sysid_param_vector_t *pv);

/** Set all parameters to a constant value */
void sysid_param_set_all(sysid_param_vector_t *pv, double value);

/** Copy parameter vector (dst must be pre-allocated, same size) */
void sysid_param_copy(sysid_param_vector_t *dst, const sysid_param_vector_t *src);

/** L2 norm of parameter vector */
double sysid_param_norm2(const sysid_param_vector_t *pv);

/** Inner product of two parameter vectors of same length */
double sysid_param_dot(const sysid_param_vector_t *a, const sysid_param_vector_t *b);

/** Allocate data structure with given dimensions */
int  sysid_data_alloc(sysid_data_t *data, int N, int nu, int ny, double Ts);

/** Free data structure */
void sysid_data_free(sysid_data_t *data);

/** Allocate frequency data */
int  sysid_freq_data_alloc(sysid_freq_data_t *fd, int Nf);

/** Free frequency data */
void sysid_freq_data_free(sysid_freq_data_t *fd);

/** Initialize fit result to zero */
void sysid_fit_result_init(sysid_fit_result_t *result);

/** Allocate and fill stability analysis from pole array */
int  sysid_stability_alloc(sysid_stability_t *stab, int n_poles);

/** Free stability analysis */
void sysid_stability_free(sysid_stability_t *stab);

#endif /* SYSID_TYPES_H */
