/**
 * @file nonlinear_core.h
 * @brief Core type definitions and data structures for nonlinear system analysis.
 *
 * Module: mini-nonlinear-system
 * Knowledge Coverage: L1 (Definitions), L2 (Core Concepts)
 *
 * This header defines the fundamental data structures used throughout the
 * nonlinear system analysis module, including state-space representations,
 * system descriptors, trajectory data, and nonlinearity characterizations.
 *
 * Course Mapping:
 *   - MIT 6.003: Signals and Systems (nonlinear extensions)
 *   - Stanford EE263: Introduction to Linear Dynamical Systems
 *   - ETH 227-0427: Signal Processing (nonlinear methods)
 *   - Tsinghua: 信号与系统 (nonlinear systems chapter)
 *
 * Reference:
 *   - Khalil, "Nonlinear Systems", 3rd ed., Prentice Hall, 2002
 *   - Strogatz, "Nonlinear Dynamics and Chaos", Perseus, 1994
 */

#ifndef NONLINEAR_CORE_H
#define NONLINEAR_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* M_PI is not in strict C11; define it for portability */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L1: Core Definitions — Fundamental Data Structures
 * -------------------------------------------------------------------------*/

#define NL_MAX_STATE_DIM   32
#define NL_MAX_PARAMS      16
#define NL_MAX_SERIES_ORDER 8
#define NL_DEFAULT_TOLERANCE 1e-9
#define NL_MAX_ITERATIONS   5000

/**
 * @brief Scalar nonlinearity types (L1: nonlinear element classification).
 */
typedef enum {
    NL_TYPE_IDEAL_RELAY,
    NL_TYPE_DEAD_ZONE,
    NL_TYPE_SATURATION,
    NL_TYPE_BACKLASH,
    NL_TYPE_QUANTIZER,
    NL_TYPE_POLYNOMIAL,
    NL_TYPE_SIGMOID,
    NL_TYPE_RELAY_WITH_HYSTERESIS,
    NL_TYPE_COULOMB_FRICTION,
    NL_TYPE_CUBIC_STIFFNESS,
    NL_TYPE_SECTOR_BOUNDED,
    NL_TYPE_CUSTOM
} nl_type_t;

/**
 * @brief Stability category for an equilibrium point (L1).
 */
typedef enum {
    NL_STABLE_NODE,
    NL_UNSTABLE_NODE,
    NL_STABLE_FOCUS,
    NL_UNSTABLE_FOCUS,
    NL_CENTER,
    NL_SADDLE,
    NL_STABLE_DEGENERATE,
    NL_UNSTABLE_DEGENERATE,
    NL_NONHYPERBOLIC,
    NL_LYAPUNOV_STABLE,
    NL_ASYMPTOTICALLY_STABLE,
    NL_EXPONENTIALLY_STABLE,
    NL_UNSTABLE,
    NL_LIMIT_CYCLE,
    NL_CHAOTIC
} nl_stability_t;

/**
 * @brief Bifurcation type classification (L1, L2).
 */
typedef enum {
    NL_BIF_NONE,
    NL_BIF_SADDLE_NODE,
    NL_BIF_PITCHFORK,
    NL_BIF_TRANSCRITICAL,
    NL_BIF_HOPF,
    NL_BIF_PERIOD_DOUBLING,
    NL_BIF_HOMOCLINIC,
    NL_BIF_HETEROCLINIC,
    NL_BIF_SADDLE_NODE_LC,
    NL_BIF_TORUS
} nl_bifurcation_t;

typedef struct {
    double x[NL_MAX_STATE_DIM];
    size_t dim;
} nl_state_t;

typedef struct {
    nl_state_t *points;
    size_t     num_points;
    size_t     capacity;
    double     t_start;
    double     t_end;
    double     dt;
} nl_trajectory_t;

typedef struct {
    nl_type_t type;
    double    params[NL_MAX_PARAMS];
    size_t    num_params;
    double    lower_bound;
    double    upper_bound;
    double    slope_min;
    double    slope_max;
    double  (*custom_f)(double x, const double *params);
    double  (*custom_df)(double x, const double *params);
} nl_nonlinearity_t;

typedef struct {
    double A[NL_MAX_STATE_DIM][NL_MAX_STATE_DIM];
    double B[NL_MAX_STATE_DIM];
    double C[NL_MAX_STATE_DIM];
    double D;
    size_t n;
} nl_linear_system_t;

typedef struct {
    nl_linear_system_t   linear;
    nl_nonlinearity_t    nonlinearity;
    double               B_nl[NL_MAX_STATE_DIM];
    nl_state_t           state;
    double               t;
} nl_lure_system_t;

typedef struct {
    int   (*f)(const double *x, double t, const double *params,
               double *dxdt, size_t dim);
    int   (*jacobian)(const double *x, double t, const double *params,
                      double *J, size_t dim);
    size_t dim;
    double params[NL_MAX_PARAMS];
    size_t num_params;
    char   name[64];
    int    is_autonomous;
} nl_nonlinear_system_t;

typedef struct {
    double (*V)(const double *x, size_t dim);
    int    (*grad_V)(const double *x, double *grad, size_t dim);
    double (*V_dot)(const double *x, const nl_nonlinear_system_t *sys);
    size_t dim;
    char   form[32];
} nl_lyapunov_function_t;

typedef struct {
    double amplitude;
    double frequency;
    double P;
    double Q;
    double magnitude;
    double phase;
} nl_describing_function_t;

typedef struct {
    nl_bifurcation_t type;
    double           param_value;
    double           eigenvalue_real;
    double           eigenvalue_imag;
    nl_state_t       equilibrium;
    int              is_subcritical;
} nl_bifurcation_point_t;

typedef struct {
    nl_trajectory_t  orbit;
    double           period;
    double           floquet_mult;
    int              is_stable;
    int              is_non_trivial;
} nl_limit_cycle_t;

typedef struct {
    double   lyapunov_exponent;
    int      lyapunov_dimension;
    double   correlation_dim;
    double   ks_entropy;
    int      is_chaotic;
} nl_chaos_metrics_t;

/* ---------------------------------------------------------------------------
 * L3: Mathematical Structures
 * -------------------------------------------------------------------------*/

typedef struct {
    double  matrix[NL_MAX_STATE_DIM][NL_MAX_STATE_DIM];
    double  eigenvalues_real[NL_MAX_STATE_DIM];
    double  eigenvalues_imag[NL_MAX_STATE_DIM];
    double  eigenvectors[NL_MAX_STATE_DIM][NL_MAX_STATE_DIM];
    size_t  dim;
    int     has_complex_pair;
} nl_jacobian_t;

typedef struct {
    double   x_min, x_max, y_min, y_max;
    size_t   nx, ny;
    double  *dxdt;
    double  *dydt;
    nl_trajectory_t *trajectories;
    size_t   num_trajectories;
    double  *nullcline_x;
    double  *nullcline_y;
    size_t   num_nullcline_x;
    size_t   num_nullcline_y;
    nl_state_t *equilibria;
    size_t      num_equilibria;
} nl_phase_portrait_t;

typedef struct {
    double  *data;
    size_t   order;
    size_t  *dims;
    size_t   total_size;
    double   norm_l2;
} nl_volterra_kernel_t;

typedef struct {
    nl_volterra_kernel_t *kernels;
    size_t   max_order;
    double   sampling_rate;
    int      is_causal;
    int      is_symmetric;
} nl_volterra_series_t;

/* ---------------------------------------------------------------------------
 * API: Allocation
 * -------------------------------------------------------------------------*/

int nl_system_init(nl_nonlinear_system_t *sys, size_t dim,
                   int (*f)(const double*, double, const double*,
                            double*, size_t),
                   const char *name);

int nl_trajectory_init(nl_trajectory_t *traj, size_t capacity, size_t dim);

void nl_trajectory_free(nl_trajectory_t *traj);

int nl_trajectory_append(nl_trajectory_t *traj, const nl_state_t *state);

int nl_nonlinearity_init(nl_nonlinearity_t *nl, nl_type_t type,
                         const double *params, size_t nparam);

double nl_nonlinearity_eval(const nl_nonlinearity_t *nl, double x);

int nl_lure_init(nl_lure_system_t *lure, const nl_linear_system_t *lin,
                 const nl_nonlinearity_t *nl, const double *Bnl,
                 const nl_state_t *x0);

int nl_lure_dxdt(const nl_lure_system_t *lure, double *dxdt);

int nl_phase_portrait_init(nl_phase_portrait_t *pp,
                           double xmin, double xmax,
                           double ymin, double ymax,
                           size_t nx, size_t ny);

void nl_phase_portrait_free(nl_phase_portrait_t *pp);

int nl_volterra_init(nl_volterra_series_t *vs, size_t order, double sr);

void nl_volterra_free(nl_volterra_series_t *vs);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_CORE_H */
