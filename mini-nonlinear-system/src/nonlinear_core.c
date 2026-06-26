/**
 * @file nonlinear_core.c
 * @brief Core implementations for nonlinear system module.
 *
 * Implements fundamental operations: system initialization, trajectory
 * management, nonlinearity evaluation, Lur'e system dynamics, phase
 * portrait setup, and Volterra series management.
 *
 * Knowledge points:
 *   L1: nonlinearity types, state/trajectory representation
 *   L2: Lur'e system dynamics, state-space nonlinear systems
 *   L5: numerical infrastructure for all higher-level algorithms
 */

#include "nonlinear_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

/* ===================================================================
 * System Initialization
 * =================================================================== */

int nl_system_init(nl_nonlinear_system_t *sys, size_t dim,
                   int (*f)(const double*, double, const double*,
                            double*, size_t),
                   const char *name)
{
    if (!sys || !f || dim < 1 || dim > NL_MAX_STATE_DIM) {
        errno = EINVAL;
        return -1;
    }
    if (!name) {
        errno = EINVAL;
        return -1;
    }

    memset(sys, 0, sizeof(*sys));
    sys->dim = dim;
    sys->f = f;
    sys->jacobian = NULL;  /* optional, set by user if available */
    sys->num_params = 0;
    sys->is_autonomous = 1;  /* default, override if needed */
    strncpy(sys->name, name, sizeof(sys->name) - 1);
    sys->name[sizeof(sys->name) - 1] = '\0';

    return 0;
}

/* ===================================================================
 * Trajectory Management
 * =================================================================== */

int nl_trajectory_init(nl_trajectory_t *traj, size_t capacity, size_t dim)
{
    if (!traj || capacity < 1 || dim < 1 || dim > NL_MAX_STATE_DIM) {
        errno = EINVAL;
        return -1;
    }

    traj->points = (nl_state_t *)malloc(capacity * sizeof(nl_state_t));
    if (!traj->points) {
        errno = ENOMEM;
        return -1;
    }

    traj->num_points = 0;
    traj->capacity = capacity;
    traj->t_start = 0.0;
    traj->t_end = 0.0;
    traj->dt = 0.0;

    /* Initialize all states with zero dimension temporarily */
    for (size_t i = 0; i < capacity; i++) {
        traj->points[i].dim = dim;
        memset(traj->points[i].x, 0, dim * sizeof(double));
    }

    return 0;
}

void nl_trajectory_free(nl_trajectory_t *traj)
{
    if (traj) {
        free(traj->points);
        traj->points = NULL;
        traj->num_points = 0;
        traj->capacity = 0;
    }
}

int nl_trajectory_append(nl_trajectory_t *traj, const nl_state_t *state)
{
    if (!traj || !state || state->dim < 1
        || state->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL;
        return -1;
    }

    /* Grow if needed (double capacity) */
    if (traj->num_points >= traj->capacity) {
        size_t new_cap = traj->capacity * 2;
        if (new_cap < 4) new_cap = 4;
        nl_state_t *new_pts = (nl_state_t *)realloc(
            traj->points, new_cap * sizeof(nl_state_t));
        if (!new_pts) {
            errno = ENOMEM;
            return -1;
        }
        traj->points = new_pts;
        traj->capacity = new_cap;
        /* Initialize new slots */
        for (size_t i = traj->num_points; i < new_cap; i++) {
            traj->points[i].dim = state->dim;
            memset(traj->points[i].x, 0, state->dim * sizeof(double));
        }
    }

    /* Copy state */
    traj->points[traj->num_points].dim = state->dim;
    memcpy(traj->points[traj->num_points].x, state->x,
           state->dim * sizeof(double));
    traj->num_points++;

    return 0;
}

/* ===================================================================
 * Nonlinearity: Definition and Evaluation
 * =================================================================== */

int nl_nonlinearity_init(nl_nonlinearity_t *nl, nl_type_t type,
                         const double *params, size_t nparam)
{
    if (!nl || nparam > NL_MAX_PARAMS) {
        errno = EINVAL;
        return -1;
    }

    memset(nl, 0, sizeof(*nl));
    nl->type = type;
    nl->num_params = nparam;

    if (params && nparam > 0) {
        memcpy(nl->params, params, nparam * sizeof(double));
    }

    /* Set default sector bounds based on type */
    switch (type) {
    case NL_TYPE_SATURATION:
        nl->slope_min = 0.0;
        nl->slope_max = (nparam >= 1) ? params[0] : 1.0;
        nl->lower_bound = -100.0;
        nl->upper_bound = 100.0;
        break;
    case NL_TYPE_DEAD_ZONE:
        nl->slope_min = 0.0;
        nl->slope_max = 1.0;
        break;
    case NL_TYPE_IDEAL_RELAY:
        nl->slope_min = 0.0;
        nl->slope_max = INFINITY;
        break;
    case NL_TYPE_SECTOR_BOUNDED:
        nl->slope_min = (nparam >= 1) ? params[0] : 0.0;
        nl->slope_max = (nparam >= 2) ? params[1] : 1.0;
        break;
    default:
        nl->slope_min = -INFINITY;
        nl->slope_max = INFINITY;
        break;
    }

    nl->custom_f = NULL;
    nl->custom_df = NULL;

    return 0;
}

double nl_nonlinearity_eval(const nl_nonlinearity_t *nl, double x)
{
    if (!nl) {
        errno = EINVAL;
        return 0.0;
    }

    /* Handle special cases uniformly */
    switch (nl->type) {
    case NL_TYPE_IDEAL_RELAY: {
        /* f(x) = M * sgn(x), param[0] = M */
        double M = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        return (x > 0.0) ? M : ((x < 0.0) ? -M : 0.0);
    }

    case NL_TYPE_DEAD_ZONE: {
        /* f(x) = 0 for |x| ≤ δ, x - δ*sgn(x) for |x| > δ */
        double delta = (nl->num_params >= 1) ? nl->params[0] : 0.1;
        if (fabs(x) <= delta) return 0.0;
        return (x > 0.0) ? (x - delta) : (x + delta);
    }

    case NL_TYPE_SATURATION: {
        /* f(x): slope k for |x| ≤ a/k, then a*sgn(x) */
        double k = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        double a = (nl->num_params >= 2) ? nl->params[1] : 1.0;
        double limit = a / k;
        if (fabs(x) <= limit) return k * x;
        return (x > 0.0) ? a : -a;
    }

    case NL_TYPE_BACKLASH: {
        /* Simplified static backlash model (full model needs state) */
        double gap = (nl->num_params >= 1) ? nl->params[0] : 0.1;
        double slope = (nl->num_params >= 2) ? nl->params[1] : 1.0;
        if (fabs(x) <= gap / 2.0) return 0.0;
        return slope * (x - copysign(gap / 2.0, x));
    }

    case NL_TYPE_QUANTIZER: {
        /* Uniform quantizer with step size Δ, param[0] = Δ */
        double delta = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        if (delta <= 0.0) return x;
        return delta * floor(x / delta + 0.5);
    }

    case NL_TYPE_POLYNOMIAL: {
        /* f(x) = Σ_{i=0}^{n-1} a_i x^i */
        double result = 0.0;
        double x_pow = 1.0;  /* x^0 = 1 */
        for (size_t i = 0; i < nl->num_params; i++) {
            result += nl->params[i] * x_pow;
            x_pow *= x;
        }
        return result;
    }

    case NL_TYPE_SIGMOID: {
        /* f(x) = tanh(α x), param[0] = α (sigmoid steepness) */
        double alpha = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        return tanh(alpha * x);
    }

    case NL_TYPE_RELAY_WITH_HYSTERESIS: {
        /* Simplified: M*sgn(x) for |x| > h, 0 otherwise */
        double M = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        double h = (nl->num_params >= 2) ? nl->params[1] : 0.1;
        if (fabs(x) > h) return (x > 0.0) ? M : -M;
        return 0.0;
    }

    case NL_TYPE_COULOMB_FRICTION: {
        /* f(v) = F_c * sgn(v), param[0] = F_c */
        double Fc = (nl->num_params >= 1) ? nl->params[0] : 0.5;
        if (fabs(x) < 1e-12) return 0.0;
        return Fc * ((x > 0.0) ? 1.0 : -1.0);
    }

    case NL_TYPE_CUBIC_STIFFNESS: {
        /* f(x) = k₁ x + k₃ x³, params = [k₁, k₃] */
        double k1 = (nl->num_params >= 1) ? nl->params[0] : 1.0;
        double k3 = (nl->num_params >= 2) ? nl->params[1] : 0.0;
        return k1 * x + k3 * x * x * x;
    }

    case NL_TYPE_SECTOR_BOUNDED: {
        /* f(x) is bounded between k_min * x and k_max * x */
        /* Return the midpoint as a sensible default */
        double k_mid = (nl->slope_min + nl->slope_max) / 2.0;
        if (!isfinite(k_mid)) k_mid = 0.5;
        return k_mid * x;
    }

    case NL_TYPE_CUSTOM: {
        if (nl->custom_f) {
            return nl->custom_f(x, nl->params);
        }
        errno = ENOSYS;
        return 0.0;
    }

    default:
        errno = EINVAL;
        return 0.0;
    }
}

/* ===================================================================
 * Lur'e System
 * =================================================================== */

int nl_lure_init(nl_lure_system_t *lure, const nl_linear_system_t *lin,
                 const nl_nonlinearity_t *nl, const double *Bnl,
                 const nl_state_t *x0)
{
    if (!lure || !lin || !nl || !Bnl || !x0) {
        errno = EINVAL;
        return -1;
    }
    if (lin->n < 1 || lin->n > NL_MAX_STATE_DIM) {
        errno = EINVAL;
        return -1;
    }

    memset(lure, 0, sizeof(*lure));
    lure->linear = *lin;
    lure->nonlinearity = *nl;
    memcpy(lure->B_nl, Bnl, lin->n * sizeof(double));
    lure->state = *x0;
    lure->t = 0.0;

    return 0;
}

int nl_lure_dxdt(const nl_lure_system_t *lure, double *dxdt)
{
    if (!lure || !dxdt) {
        errno = EINVAL;
        return -1;
    }

    size_t n = lure->linear.n;

    /*
     * Lur'e system dynamics:
     *   y = C^T x
     *   φ = nonlinearity(y)
     *   dx/dt = A x + B_nl φ
     */
    /* Compute output y = C^T x */
    double y = 0.0;
    for (size_t i = 0; i < n; i++) {
        y += lure->linear.C[i] * lure->state.x[i];
    }

    /* Evaluate nonlinearity φ(y) */
    double phi = nl_nonlinearity_eval(&lure->nonlinearity, y);

    /* Compute dx/dt = A x + B_nl φ */
    for (size_t i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (size_t j = 0; j < n; j++) {
            dxdt[i] += lure->linear.A[i][j] * lure->state.x[j];
        }
        dxdt[i] += lure->B_nl[i] * phi;
    }

    return 0;
}

/* ===================================================================
 * Phase Portrait
 * =================================================================== */

int nl_phase_portrait_init(nl_phase_portrait_t *pp,
                           double xmin, double xmax,
                           double ymin, double ymax,
                           size_t nx, size_t ny)
{
    if (!pp || nx < 2 || ny < 2 || xmax <= xmin || ymax <= ymin) {
        errno = EINVAL;
        return -1;
    }

    memset(pp, 0, sizeof(*pp));
    pp->x_min = xmin;
    pp->x_max = xmax;
    pp->y_min = ymin;
    pp->y_max = ymax;
    pp->nx = nx;
    pp->ny = ny;

    /* Allocate vector field arrays */
    size_t total = nx * ny;
    pp->dxdt = (double *)malloc(total * sizeof(double));
    pp->dydt = (double *)malloc(total * sizeof(double));

    if (!pp->dxdt || !pp->dydt) {
        free(pp->dxdt);
        free(pp->dydt);
        pp->dxdt = NULL;
        pp->dydt = NULL;
        errno = ENOMEM;
        return -1;
    }

    /* Initialize to zero */
    memset(pp->dxdt, 0, total * sizeof(double));
    memset(pp->dydt, 0, total * sizeof(double));

    pp->trajectories = NULL;
    pp->num_trajectories = 0;
    pp->nullcline_x = NULL;
    pp->nullcline_y = NULL;
    pp->num_nullcline_x = 0;
    pp->num_nullcline_y = 0;
    pp->equilibria = NULL;
    pp->num_equilibria = 0;

    return 0;
}

void nl_phase_portrait_free(nl_phase_portrait_t *pp)
{
    if (!pp) return;

    free(pp->dxdt);
    free(pp->dydt);
    pp->dxdt = NULL;
    pp->dydt = NULL;

    /* Free trajectories */
    if (pp->trajectories) {
        for (size_t i = 0; i < pp->num_trajectories; i++) {
            nl_trajectory_free(&pp->trajectories[i]);
        }
        free(pp->trajectories);
        pp->trajectories = NULL;
    }

    free(pp->nullcline_x);
    free(pp->nullcline_y);
    pp->nullcline_x = NULL;
    pp->nullcline_y = NULL;

    free(pp->equilibria);
    pp->equilibria = NULL;
    pp->num_equilibria = 0;
}

/* ===================================================================
 * Volterra Series Management
 * =================================================================== */

int nl_volterra_init(nl_volterra_series_t *vs, size_t order, double sr)
{
    if (!vs || order < 1 || order > NL_MAX_SERIES_ORDER || sr <= 0.0) {
        errno = EINVAL;
        return -1;
    }

    memset(vs, 0, sizeof(*vs));
    vs->max_order = order;
    vs->sampling_rate = sr;
    vs->is_causal = 1;
    vs->is_symmetric = 1;

    vs->kernels = (nl_volterra_kernel_t *)calloc(
        order, sizeof(nl_volterra_kernel_t));
    if (!vs->kernels) {
        errno = ENOMEM;
        return -1;
    }

    return 0;
}

void nl_volterra_free(nl_volterra_series_t *vs)
{
    if (!vs) return;

    if (vs->kernels) {
        for (size_t i = 0; i < vs->max_order; i++) {
            free(vs->kernels[i].data);
            free(vs->kernels[i].dims);
        }
        free(vs->kernels);
        vs->kernels = NULL;
    }
    vs->max_order = 0;
}
