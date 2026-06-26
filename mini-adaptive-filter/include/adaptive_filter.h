#ifndef ADAPTIVE_FILTER_H
#define ADAPTIVE_FILTER_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Core Definitions - Enumerations and Type Aliases
 * ========================================================================= */

/** Filter structure topology */
typedef enum {
    AF_STRUCTURE_FIR          = 0,
    AF_STRUCTURE_IIR          = 1,
    AF_STRUCTURE_LATTICE      = 2,
    AF_STRUCTURE_SUBBAND      = 3,
    AF_STRUCTURE_FREQ_DOMAIN  = 4,
    AF_STRUCTURE_VOLTERRA     = 5
} af_structure_t;

/** Adaptation algorithm identifiers */
typedef enum {
    AF_ALGO_LMS              = 0,
    AF_ALGO_NLMS             = 1,
    AF_ALGO_SIGN_ERROR_LMS   = 2,
    AF_ALGO_SIGN_DATA_LMS    = 3,
    AF_ALGO_SIGN_SIGN_LMS    = 4,
    AF_ALGO_LEAKY_LMS        = 5,
    AF_ALGO_BLOCK_LMS        = 6,
    AF_ALGO_RLS              = 100,
    AF_ALGO_QR_RLS           = 101,
    AF_ALGO_SQUARE_ROOT_RLS  = 102,
    AF_ALGO_HOUSEHOLDER_RLS  = 103,
    AF_ALGO_APA              = 200,
    AF_ALGO_KALMAN_LINEAR    = 300,
    AF_ALGO_KALMAN_EXTENDED  = 301,
    AF_ALGO_KALMAN_UNSCENTED = 302,
    AF_ALGO_PARTICLE         = 400
} af_algorithm_t;

/** Adaptation mode for reference signal availability */
typedef enum {
    AF_MODE_TRAINING          = 0,
    AF_MODE_DECISION_DIRECTED = 1,
    AF_MODE_BLIND             = 2
} af_mode_t;

/** Error norm specification */
typedef enum {
    AF_ERROR_L2   = 0,
    AF_ERROR_L1   = 1,
    AF_ERROR_LINF = 2,
    AF_ERROR_RMS  = 3
} af_error_norm_t;

/** Tap-weight initialization strategy */
typedef enum {
    AF_INIT_ZEROS        = 0,
    AF_INIT_RANDOM       = 1,
    AF_INIT_CENTER_TAP   = 2,
    AF_INIT_USER_DEFINED = 3
} af_init_strategy_t;

/* =========================================================================
 * L1: Adaptive Filter Configuration and State Structures
 * ========================================================================= */

/**
 * @brief Real-valued adaptive filter configuration
 *
 * Encapsulates all hyperparameters that govern the adaptive algorithm.
 * A single config can be shared across multiple filter instances.
 *
 * Parameter semantics:
 *   mu: Step-size, 0 < mu < 2/lambda_max for LMS convergence in mean.
 *       Larger mu => faster convergence, higher steady-state error.
 *   lambda: RLS forgetting factor, 0 < lambda <= 1.
 *           lambda=1 = infinite memory; lambda<1 = tracking capability.
 *   delta: RLS regularization (initial diagonal loading).
 *          Prevents singularity: P(0) = delta^{-1} * I.
 *   gamma_leakage: Leakage factor, 0 < gamma <= 1.
 *   epsilon: NLMS stability constant, prevents div-by-zero.
 *   projection_order: APA projection order (1 = NLMS).
 */
typedef struct {
    af_structure_t     structure;
    af_algorithm_t     algorithm;
    af_mode_t          mode;
    af_init_strategy_t init_strategy;
    size_t             order;
    double             mu;
    double             lambda;
    double             delta;
    double             gamma_leakage;
    double             epsilon;
    size_t             block_size;
    size_t             projection_order;
    double             convergence_tol;
    size_t             max_iterations;
} af_config_t;

/**
 * @brief Real-valued adaptive filter runtime state
 *
 * Wiener solution (optimal, not directly computable online):
 *   w_opt = R^{-1} * p,  R = E[x(n)*x^T(n)], p = E[d(n)*x(n)]
 *
 * Instantaneous error:
 *   e(n) = d(n) - y(n) = d(n) - w^T(n) * x(n)
 *
 * Tapped-delay-line operation:
 *   x(n) = [x(n), x(n-1), ..., x(n-M+1)]^T
 *   y(n) = w^T(n) * x(n)
 *   e(n) = d(n) - y(n)
 */
typedef struct {
    const af_config_t *config;
    double  *w;                  /* Tap-weight vector [order] */
    double  *x_buffer;           /* Input delay line [order] */
    double  *d_buffer;           /* Desired signal buffer */
    double   error;              /* e(n) */
    double   output;             /* y(n) */
    double   mse;                /* Running MSE estimate */
    double   emse;               /* Excess MSE estimate */
    double   misadjustment;      /* M = EMSE / J_min */
    size_t   iteration;          /* Iteration counter */
    double   tap_norm_delta;     /* ||w(n)-w(n-1)|| */
    double   gradient_norm;      /* ||gradient|| */
    double  *P;                  /* Inverse correlation [order^2] (RLS) */
    double  *K;                  /* Kalman gain [order] (RLS) */
    double  *pi_vector;          /* Auxiliary vector (APA) */
    double  *g_vector;           /* Gradient vector (block LMS) */
    double  *residual_vector;    /* Residual workspace */
    void    *user_data;          /* Opaque user data */
    void   (*error_cb)(size_t, double, double, void*);
} af_filter_t;

/**
 * @brief Complex-valued adaptive filter configuration
 */
typedef struct {
    af_structure_t     structure;
    af_algorithm_t     algorithm;
    af_mode_t          mode;
    af_init_strategy_t init_strategy;
    size_t             order;
    double             mu;
    double             lambda;
    double             delta;
    double             gamma_leakage;
    double             epsilon;
    size_t             block_size;
    size_t             projection_order;
    double             convergence_tol;
    size_t             max_iterations;
} af_config_complex_t;

/**
 * @brief Complex-valued adaptive filter runtime state
 *
 * Complex LMS update (Widrow, 1975):
 *   w(n+1) = w(n) + mu * conj(e(n)) * x(n)
 *
 * Complex output: y(n) = w^H(n) * x(n) = sum conj(w_k)*x_k
 */
typedef struct {
    const af_config_complex_t *config;
    double  *w_re, *w_im;
    double  *x_re_buffer, *x_im_buffer;
    double   error_re, error_im;
    double   output_re, output_im;
    double   mse, emse, misadjustment;
    size_t   iteration;
    double  *P_re, *P_im;
    double  *K_re, *K_im;
} af_filter_complex_t;

/* =========================================================================
 * L1: Signal Statistics Structure
 * ========================================================================= */

/**
 * @brief Signal statistics for optimal Wiener filter analysis
 *
 * For WSS processes, R is symmetric Toeplitz: R[i][j] = r_xx[|i-j|]
 * Eigenvalue spread chi(R) = lambda_max / lambda_min determines
 * convergence speed of gradient-based algorithms.
 */
typedef struct {
    size_t   length;
    double  *r_xx;
    double  *p_dx;
    double   var_d;
    double   var_noise;
    double   eigenvalue_min;
    double   eigenvalue_max;
    double   eigenvalue_spread;
    double   condition_number;
} af_statistics_t;

/* =========================================================================
 * L1: Performance Metrics Structure
 * ========================================================================= */

/**
 * @brief Adaptive filter performance metrics
 *
 * Learning curve:
 *   MSE(n) = E[|e(n)|^2]
 *   J_min  = sigma_d^2 - p^T * R^{-1} * p
 *   J_excess = lim MSE(n) - J_min
 *   M = J_excess / J_min
 */
typedef struct {
    double   mse_initial;
    double   mse_final;
    double   mse_minimum;
    double   excess_mse;
    double   misadjustment;
    double   mse_time_constant;
    double   convergence_rate;
    size_t   convergence_time;
    double   tracking_error;
    double   snr_improvement;
} af_performance_t;

/* =========================================================================
 * L2: Wiener Filter Structures
 * ========================================================================= */

typedef struct {
    size_t   order;
    double  *r_xx;
    double  *p_dx;
    double   sigma_d2;
} af_wiener_config_t;

/**
 * @brief Wiener filter solution
 *
 * Wiener-Hopf equation: R * w_opt = p
 * Minimum MSE: J_min = sigma_d^2 - p^T * w_opt
 * Orthogonality principle: E[e_opt(n) * x(n-k)] = 0
 */
typedef struct {
    size_t   order;
    double  *w_opt;
    double   j_min;
    double  *r_xx;
    double  *p_dx;
    double   sigma_d2;
    int      well_conditioned;
} af_wiener_filter_t;

/* =========================================================================
 * L3: Mathematical Structures - Matrix and Vector
 * ========================================================================= */

/** Symmetric Toeplitz matrix (compact O(N) storage) */
typedef struct {
    size_t   order;
    double  *first_col;
} af_toeplitz_t;

/** General real matrix (row-major) */
typedef struct {
    size_t   rows;
    size_t   cols;
    size_t   stride;
    double  *data;
} af_matrix_t;

/** Real vector */
typedef struct {
    size_t   length;
    double  *data;
} af_vector_t;

/* =========================================================================
 * L4: Fundamental Laws - LMS Convergence Bounds
 * ========================================================================= */

/**
 * @brief LMS convergence bound constants
 *
 * Convergence in mean:     0 < mu < 2 / lambda_max
 * Convergence in mean-sq:  0 < mu < 2 / (3 * tr(R))
 *
 * Misadjustment (small mu): M approx (mu/2) * tr(R)
 */
typedef struct {
    double   lambda_max;
    double   lambda_min;
    double   trace_R;
    double   mu_max_mean;
    double   mu_max_mean_sq;
    double   theoretical_mu_opt;
} af_convergence_bounds_t;

/* =========================================================================
 * Error Callback Type
 * ========================================================================= */

typedef void (*af_error_callback_t)(size_t iteration, double error,
                                     double mse, void *user_data);


/* =========================================================================
 * Core API Declarations
 * ========================================================================= */

/* Lifecycle */
int  af_init(const af_config_t *config, af_filter_t *filter);
void af_free(af_filter_t *filter);
int  af_init_complex(const af_config_complex_t *config,
                      af_filter_complex_t *filter);
void af_free_complex(af_filter_complex_t *filter);

/* Adaptation */
double af_adapt(af_filter_t *filter, double x, double d);
void   af_adapt_complex(af_filter_complex_t *filter,
                         double x_re, double x_im,
                         double d_re, double d_im,
                         double *out_re, double *out_im);
void   af_adapt_batch(af_filter_t *filter, const double *x, const double *d,
                       size_t n_samples, double *errors);

/* Filtering without adaptation */
double af_filter_output(const af_filter_t *filter, double x);
void   af_reset(af_filter_t *filter);
void   af_reset_complex(af_filter_complex_t *filter);

/* Diagnostics */
void   af_get_weights(const af_filter_t *filter, double *w_out);
void   af_set_weights(af_filter_t *filter, const double *w_in);
double af_get_misadjustment(const af_filter_t *filter);
void   af_set_stepsize(af_filter_t *filter, double new_mu);
void   af_set_error_callback(af_filter_t *filter,
                              af_error_callback_t callback, void *user_data);

/* Wiener filter */
int  af_wiener_solve(const af_wiener_config_t *cfg,
                      af_wiener_filter_t *wiener);
void af_wiener_free(af_wiener_filter_t *wiener);

/* Autocorrelation estimation */
void af_autocorrelation(const double *x, size_t n_samples, size_t max_lag,
                         double *r_xx_out, int biased);

/* Matrix/vector utilities */
int  af_matrix_alloc(size_t rows, size_t cols, af_matrix_t *mat);
void af_matrix_free(af_matrix_t *mat);
void af_matrix_set(af_matrix_t *mat, size_t i, size_t j, double val);
double af_matrix_get(const af_matrix_t *mat, size_t i, size_t j);
void af_matrix_vector_mul(const af_matrix_t *A, const af_vector_t *x,
                           af_vector_t *y);
void af_matrix_transpose(const af_matrix_t *A, af_matrix_t *AT);
int  af_matrix_inverse(const af_matrix_t *A, af_matrix_t *A_inv);
void af_matrix_identity(af_matrix_t *A);
void af_matrix_cholesky(const af_matrix_t *A, af_matrix_t *L);

int  af_vector_alloc(size_t length, af_vector_t *vec);
void af_vector_free(af_vector_t *vec);
double af_vector_dot(const double *a, const double *b, size_t n);
double af_vector_norm(const double *a, size_t n);
void af_vector_scale(double *a, double alpha, size_t n);
void af_vector_add(const double *a, const double *b, double *c, size_t n);
void af_vector_sub(const double *a, const double *b, double *c, size_t n);

/* Toeplitz solver (Levinson-Durbin) */
int  af_toeplitz_solve(const double *r_xx, const double *p_dx,
                        size_t order, double *w_opt);

/* Convergence analysis */
int  af_compute_convergence_bounds(const double *r_xx, size_t order,
                                    af_convergence_bounds_t *bounds);
void af_performance_evaluate(af_filter_t *filter,
                              const double *x, const double *d,
                              size_t n_samples, af_performance_t *perf);

/* LMS family */
void af_lms_direct_update(af_filter_t *filter, const double *x, double d,
                           double *w, double mu, size_t order);
void af_lms_normalized_update(af_filter_t *filter, const double *x, double d,
                               double *w, double mu, double eps, size_t order);
void af_lms_sign_error_update(af_filter_t *filter, const double *x, double d,
                               double *w, double mu, size_t order);
void af_lms_sign_data_update(af_filter_t *filter, const double *x, double d,
                              double *w, double mu, size_t order);
void af_lms_sign_sign_update(af_filter_t *filter, const double *x, double d,
                              double *w, double mu, size_t order);
void af_lms_leaky_update(af_filter_t *filter, const double *x, double d,
                          double *w, double mu, double gamma, size_t order);

/* RLS family */
void af_rls_standard_update(af_filter_t *filter, double x, double d,
                             double lambda, double delta, size_t order);
void af_rls_qr_givens_update(af_filter_t *filter, double x, double d,
                              double lambda, size_t order);
double af_rls_compute_kalman_gain(const af_matrix_t *P, const double *x,
                                   double lambda, size_t order, double *K);
void af_rls_update_inverse_correlation(af_matrix_t *P, const double *K,
                                        const double *x, double lambda,
                                        size_t order);

/* Affine Projection Algorithm */
void af_apa_update(af_filter_t *filter, const double *x, double d,
                    double mu, double epsilon, size_t order, size_t P_order);

/* Kalman filter */
void af_kalman_predict(double *x_hat, double *P, const double *F,
                        const double *Q, size_t n_states);
void af_kalman_update(double *x_hat, double *P, const double *H,
                       const double *R, double z, size_t n_states,
                       size_t n_meas);
void af_kalman_extended_predict(double *x_hat, double *P,
                                 void (*f)(const double*, double*),
                                 void (*F_jac)(const double*, double*),
                                 const double *Q, size_t n_states);
void af_kalman_extended_update(double *x_hat, double *P,
                                void (*h)(const double*, double*),
                                void (*H_jac)(const double*, double*),
                                const double *R, double z,
                                size_t n_states, size_t n_meas);

/* Eigenvalue computation */
int  af_eigenvalue_max(const double *r_xx, size_t order, double *lambda,
                        size_t max_iter, double tol);
int  af_eigenvalue_spread(const double *r_xx, size_t order,
                           double *lambda_min, double *lambda_max,
                           double *spread);

/* =========================================================================
 * Application-level API (L6-L7)
 * ========================================================================= */

/* System identification using LMS */
int  af_app_system_identify_lms(const double *x, const double *d,
                                 size_t n_samples, size_t order,
                                 double mu, double *w_out, double *mse_out);

/* Adaptive noise cancellation using LMS */
int  af_app_noise_cancel_lms(const double *primary, const double *reference,
                              size_t n_samples, size_t order, double mu,
                              double *cleaned_out);

/* Acoustic echo cancellation using NLMS */
int  af_app_echo_cancel_nlms(const double *far_end, const double *mic_signal,
                              size_t n_samples, size_t order, double mu,
                              double *residual_out);

/* Channel equalization using LMS (decision-directed) */
int  af_app_channel_equalize_lms(const double *received,
                                  const double *training,
                                  size_t n_samples, size_t train_len,
                                  size_t order, double mu,
                                  double *equalized_out);

/* Adaptive line enhancement using LMS */
int  af_app_line_enhance_lms(const double *input, size_t n_samples,
                              size_t order, double mu, size_t delay,
                              double *enhanced_out);

/* ECG 60 Hz noise cancellation (L7 biomedical application) */
int  af_app_ecg_noise_cancel(const double *ecg_clean, const double *noise_ref,
                              size_t n_samples, double *cleaned_out);

/* Block LMS gradient accumulation */
void af_lms_block_gradient(af_filter_t *filter,
                            const double *block_x,
                            const double *errors,
                            size_t block_size,
                            size_t order,
                            double *grad_out);

/* RLS system identification */
void af_rls_system_identify(const double *x, const double *d,
                             size_t n_samples, size_t order,
                             double lambda, double delta,
                             double *w_estimated);

/* LMS system identification (internal helper) */
void af_lms_system_identify(const double *x, const double *d,
                             size_t n_samples, size_t order, double mu,
                             double *w_estimated);

/* DC motor EKF helpers (L6 application) */
void af_kalman_dc_motor_f(const double *x, double *x_next);
void af_kalman_dc_motor_F_jac(const double *x, double *F);

#ifdef __cplusplus
}
#endif

#endif /* ADAPTIVE_FILTER_H */
