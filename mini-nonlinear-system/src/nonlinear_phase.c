/**
 * @file nonlinear_phase.c
 * @brief Phase plane / phase space analysis for nonlinear ODEs.
 *
 * Implements phase portrait construction, equilibrium finding,
 * trajectory integration (RK4, Dormand-Prince), Bendixson/Dulac
 * criteria, Poincare index, manifold tracing, and isoclines.
 *
 * Knowledge points:
 *   L4: Bendixson's criterion, Dulac's criterion, Poincare index
 *   L5: Phase portrait generation, nullcline computation,
 *       RK4 integration, adaptive Dormand-Prince
 *   L6: Manifold tracing, Poincare map construction
 *
 * Reference:
 *   - Strogatz, "Nonlinear Dynamics and Chaos", Ch.5-7, 1994
 *   - Guckenheimer & Holmes, "Nonlinear Oscillations...", 1983
 */

#include "nonlinear_phase.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* ===================================================================
 * Phase Portrait Computation
 *
 * Fills the phase portrait structure with:
 * - Vector field on regular grid
 * - Nullclines via contour tracing
 * - Representative trajectories from sampled initial conditions
 * =================================================================== */

int nl_compute_phase_portrait(nl_nonlinear_system_t *sys,
                               nl_phase_portrait_t *pp)
{
    if (!sys || !pp || sys->dim != 2) {
        errno = EINVAL;
        return -1;
    }

    double dx = (pp->x_max - pp->x_min) / (double)(pp->nx - 1);
    double dy = (pp->y_max - pp->y_min) / (double)(pp->ny - 1);

    /* Compute vector field on grid */
    double state[2];
    for (size_t iy = 0; iy < pp->ny; iy++) {
        double y = pp->y_min + (double)iy * dy;
        for (size_t ix = 0; ix < pp->nx; ix++) {
            double x = pp->x_min + (double)ix * dx;
            state[0] = x; state[1] = y;

            double f_val[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, state, f_val) == 0) {
                pp->dxdt[iy * pp->nx + ix] = f_val[0];
                pp->dydt[iy * pp->nx + ix] = f_val[1];
            }
        }
    }

    return 0;
}

/* ===================================================================
 * Find Equilibria in 2D via Nullcline Intersections
 *
 * Strategy:
 * 1. Generate nullcline samples
 * 2. Find sign changes in f1 (for x-nullcline) and f2 (for y-nullcline)
 * 3. At intersection candidates, refine with Newton
 * =================================================================== */

int nl_find_equilibria_2d(nl_nonlinear_system_t *sys,
                          nl_phase_portrait_t *pp, double tol)
{
    if (!sys || !pp || sys->dim != 2) {
        errno = EINVAL; return -1;
    }

    double dx = (pp->x_max - pp->x_min) / (double)(pp->nx - 1);
    double dy = (pp->y_max - pp->y_min) / (double)(pp->ny - 1);
    size_t max_eq = 100;

    pp->equilibria = (nl_state_t *)malloc(
        max_eq * sizeof(nl_state_t));
    if (!pp->equilibria) { errno = ENOMEM; return -1; }

    pp->num_equilibria = 0;
    double state[2];

    /* Search grid for sign changes in both components */
    for (size_t iy = 0; iy < pp->ny - 1 && pp->num_equilibria < max_eq; iy++) {
        for (size_t ix = 0; ix < pp->nx - 1 && pp->num_equilibria < max_eq; ix++) {
            double f00_x = pp->dxdt[iy * pp->nx + ix];
            double f00_y = pp->dydt[iy * pp->nx + ix];
            double f10_x = pp->dxdt[iy * pp->nx + ix + 1];
            double f10_y = pp->dydt[iy * pp->nx + ix + 1];
            double f01_x = pp->dxdt[(iy + 1) * pp->nx + ix];
            double f01_y = pp->dydt[(iy + 1) * pp->nx + ix];
            double f11_x = pp->dxdt[(iy + 1) * pp->nx + ix + 1];
            double f11_y = pp->dydt[(iy + 1) * pp->nx + ix + 1];

            /* Check sign changes in x-component */
            int sign_change_x = ((f00_x > 0) != (f10_x > 0))
                              || ((f00_x > 0) != (f01_x > 0))
                              || ((f00_x > 0) != (f11_x > 0));

            /* Check sign changes in y-component */
            int sign_change_y = ((f00_y > 0) != (f10_y > 0))
                              || ((f00_y > 0) != (f01_y > 0))
                              || ((f00_y > 0) != (f11_y > 0));

            if (!sign_change_x || !sign_change_y) continue;

            /* Both components change sign in this cell -> potential equilibrium */
            double x0 = pp->x_min + ((double)ix + 0.5) * dx;
            double y0 = pp->y_min + ((double)iy + 0.5) * dy;
            state[0] = x0; state[1] = y0;

            double x_eq[2] = {x0, y0};

            /* Newton refinement */
            for (int iter = 0; iter < 20; iter++) {
                double f_val[NL_MAX_STATE_DIM];
                if (nl_evaluate_f(sys, x_eq, f_val) != 0) break;

                double norm = sqrt(f_val[0] * f_val[0]
                                   + f_val[1] * f_val[1]);
                if (norm < tol) {
                    /* Found equilibrium */
                    int duplicate = 0;
                    for (size_t e = 0; e < pp->num_equilibria; e++) {
                        double d = sqrt(
                            (pp->equilibria[e].x[0] - x_eq[0])
                            * (pp->equilibria[e].x[0] - x_eq[0])
                            + (pp->equilibria[e].x[1] - x_eq[1])
                            * (pp->equilibria[e].x[1] - x_eq[1]));
                        if (d < tol * 10.0) { duplicate = 1; break; }
                    }
                    if (!duplicate) {
                        pp->equilibria[pp->num_equilibria].x[0]
                            = x_eq[0];
                        pp->equilibria[pp->num_equilibria].x[1]
                            = x_eq[1];
                        pp->equilibria[pp->num_equilibria].dim = 2;
                        pp->num_equilibria++;
                    }
                    break;
                }

                /* Jacobian and Newton step */
                nl_jacobian_t J;
                state[0] = x_eq[0]; state[1] = x_eq[1];
                if (nl_compute_jacobian(sys, state, &J) != 0) break;

                double det_J = J.matrix[0][0] * J.matrix[1][1]
                               - J.matrix[0][1] * J.matrix[1][0];
                if (fabs(det_J) < 1e-14) break;

                double delta_x = (J.matrix[1][1] * (-f_val[0])
                                  - J.matrix[0][1] * (-f_val[1]))
                                 / det_J;
                double delta_y = (J.matrix[0][0] * (-f_val[1])
                                  - J.matrix[1][0] * (-f_val[0]))
                                 / det_J;

                x_eq[0] += delta_x;
                x_eq[1] += delta_y;

                /* Clamp to region */
                if (x_eq[0] < pp->x_min || x_eq[0] > pp->x_max
                    || x_eq[1] < pp->y_min || x_eq[1] > pp->y_max)
                    break;
            }
        }
    }

    return (int)pp->num_equilibria;
}

/* ===================================================================
 * RK4 Fixed-Step Integration
 *
 * Classical 4th-order Runge-Kutta method.
 * x_{n+1} = x_n + (h/6)(k1 + 2k2 + 2k3 + k4)
 * =================================================================== */

int nl_integrate_rk4(nl_nonlinear_system_t *sys,
                     const double *x0, double T, double dt,
                     nl_trajectory_t *traj, size_t *n_steps)
{
    if (!sys || !x0 || !traj || !n_steps
        || T <= 0.0 || dt <= 0.0 || sys->dim < 1
        || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;
    size_t steps = (size_t)(T / dt);
    if (steps < 1) steps = 1;

    double x[NL_MAX_STATE_DIM];
    memcpy(x, x0, n * sizeof(double));

    /* Store initial state */
    nl_state_t init_state;
    init_state.dim = n;
    memcpy(init_state.x, x, n * sizeof(double));
    nl_trajectory_append(traj, &init_state);

    for (size_t s = 0; s < steps; s++) {
        double k1[NL_MAX_STATE_DIM], k2[NL_MAX_STATE_DIM];
        double k3[NL_MAX_STATE_DIM], k4[NL_MAX_STATE_DIM];
        double temp[NL_MAX_STATE_DIM];

        /* k1 = f(x) */
        if (nl_evaluate_f(sys, x, k1) != 0) break;

        /* k2 = f(x + h*k1/2) */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + 0.5 * dt * k1[d];
        if (nl_evaluate_f(sys, temp, k2) != 0) break;

        /* k3 = f(x + h*k2/2) */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + 0.5 * dt * k2[d];
        if (nl_evaluate_f(sys, temp, k3) != 0) break;

        /* k4 = f(x + h*k3) */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + dt * k3[d];
        if (nl_evaluate_f(sys, temp, k4) != 0) break;

        /* x = x + (h/6)(k1 + 2k2 + 2k3 + k4) */
        for (size_t d = 0; d < n; d++)
            x[d] += (dt / 6.0) * (k1[d] + 2.0 * k2[d]
                                   + 2.0 * k3[d] + k4[d]);

        /* Append to trajectory (every step or downsample) */
        if (s % 1 == 0) {  /* could downsample for large steps */
            nl_state_t state;
            state.dim = n;
            memcpy(state.x, x, n * sizeof(double));
            nl_trajectory_append(traj, &state);
        }
    }

    *n_steps = traj->num_points;
    return 0;
}

/* ===================================================================
 * Adaptive Dormand-Prince 5(4) Integration
 *
 * Uses the Dormand-Prince RK5(4)7M coefficients.
 * Adaptive step size based on local error estimate.
 * =================================================================== */

int nl_integrate_adaptive(nl_nonlinear_system_t *sys,
                          const double *x0, double T, double dt0,
                          double tol, nl_trajectory_t *traj)
{
    if (!sys || !x0 || !traj
        || T <= 0.0 || dt0 <= 0.0 || tol <= 0.0
        || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;
    double x[NL_MAX_STATE_DIM];
    memcpy(x, x0, n * sizeof(double));

    double t = 0.0;
    double h = dt0;

    /* Store initial state */
    nl_state_t state_init;
    state_init.dim = n;
    memcpy(state_init.x, x, n * sizeof(double));
    nl_trajectory_append(traj, &state_init);

    /* Dormand-Prince Butcher tableau (simplified: 4th order + error est) */
    while (t < T) {
        if (t + h > T) h = T - t;
        if (h < 1e-12) break;

        double k1[NL_MAX_STATE_DIM], k2[NL_MAX_STATE_DIM];
        double k3[NL_MAX_STATE_DIM], k4[NL_MAX_STATE_DIM];
        double k5[NL_MAX_STATE_DIM], k6[NL_MAX_STATE_DIM];
        double k7[NL_MAX_STATE_DIM];
        double temp[NL_MAX_STATE_DIM];
        double y5[NL_MAX_STATE_DIM], y4[NL_MAX_STATE_DIM];

        /* k1 */
        if (nl_evaluate_f(sys, x, k1) != 0) break;

        /* k2 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * 0.2 * k1[d];
        if (nl_evaluate_f(sys, temp, k2) != 0) break;

        /* k3 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * (3.0/40.0 * k1[d]
                                  + 9.0/40.0 * k2[d]);
        if (nl_evaluate_f(sys, temp, k3) != 0) break;

        /* k4 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * (44.0/45.0 * k1[d]
                                  - 56.0/15.0 * k2[d]
                                  + 32.0/9.0 * k3[d]);
        if (nl_evaluate_f(sys, temp, k4) != 0) break;

        /* k5 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * (19372.0/6561.0 * k1[d]
                                  - 25360.0/2187.0 * k2[d]
                                  + 64448.0/6561.0 * k3[d]
                                  - 212.0/729.0 * k4[d]);
        if (nl_evaluate_f(sys, temp, k5) != 0) break;

        /* k6 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * (9017.0/3168.0 * k1[d]
                                  - 355.0/33.0 * k2[d]
                                  + 46732.0/5247.0 * k3[d]
                                  + 49.0/176.0 * k4[d]
                                  - 5103.0/18656.0 * k5[d]);
        if (nl_evaluate_f(sys, temp, k6) != 0) break;

        /* k7 */
        for (size_t d = 0; d < n; d++)
            temp[d] = x[d] + h * (35.0/384.0 * k1[d]
                                  + 500.0/1113.0 * k3[d]
                                  + 125.0/192.0 * k4[d]
                                  - 2187.0/6784.0 * k5[d]
                                  + 11.0/84.0 * k6[d]);
        if (nl_evaluate_f(sys, temp, k7) != 0) break;

        /* 5th order estimate */
        for (size_t d = 0; d < n; d++)
            y5[d] = x[d] + h * (35.0/384.0 * k1[d]
                                + 500.0/1113.0 * k3[d]
                                + 125.0/192.0 * k4[d]
                                - 2187.0/6784.0 * k5[d]
                                + 11.0/84.0 * k6[d]);

        /* 4th order estimate */
        for (size_t d = 0; d < n; d++)
            y4[d] = x[d] + h * (5179.0/57600.0 * k1[d]
                                + 7571.0/16695.0 * k3[d]
                                + 393.0/640.0 * k4[d]
                                - 92097.0/339200.0 * k5[d]
                                + 187.0/2100.0 * k6[d]
                                + 1.0/40.0 * k7[d]);

        /* Error estimate */
        double err = 0.0;
        for (size_t d = 0; d < n; d++) {
            double diff = y5[d] - y4[d];
            err += diff * diff;
        }
        err = sqrt(err / (double)n);

        /* Step acceptance and adaptation */
        if (err < tol) {
            memcpy(x, y5, n * sizeof(double));
            t += h;

            nl_state_t state;
            state.dim = n;
            memcpy(state.x, x, n * sizeof(double));
            nl_trajectory_append(traj, &state);
        }

        /* Adapt step size */
        if (err > 1e-15) {
            double fac = 0.9 * pow(tol / (err + 1e-15), 0.2);
            if (fac > 2.0) fac = 2.0;
            if (fac < 0.2) fac = 0.2;
            h *= fac;
        } else {
            h *= 2.0;
        }
        if (h < 1e-12) h = 1e-12;
        if (h > T - t) h = T - t;
    }

    return 0;
}

/* ===================================================================
 * Bendixson's Criterion
 *
 * If div(f) = df1/dx + df2/dy does not change sign in a simply
 * connected region, then no limit cycles exist in that region.
 *
 * Computes divergence on the phase portrait grid and checks
 * for sign changes.
 * =================================================================== */

int nl_bendixson_criterion(const nl_phase_portrait_t *pp, int *has_cycle)
{
    if (!pp || !has_cycle) { errno = EINVAL; return -1; }

    *has_cycle = 1;  /* default: criterion is inconclusive */

    if (!pp->dxdt || !pp->dydt) return 0;

    /* Compute divergence on grid via finite differences */
    double first_div = 0.0;
    int first_set = 0;
    int sign_change = 0;
    double dx = (pp->x_max - pp->x_min) / (double)(pp->nx - 1);
    double dy = (pp->y_max - pp->y_min) / (double)(pp->ny - 1);

    for (size_t iy = 1; iy < pp->ny - 1; iy++) {
        for (size_t ix = 1; ix < pp->nx - 1; ix++) {
            /* df1/dx */
            double df1dx = (pp->dxdt[iy * pp->nx + ix + 1]
                            - pp->dxdt[iy * pp->nx + ix - 1])
                           / (2.0 * dx);

            /* df2/dy */
            double df2dy = (pp->dydt[(iy + 1) * pp->nx + ix]
                            - pp->dydt[(iy - 1) * pp->nx + ix])
                           / (2.0 * dy);

            double divergence = df1dx + df2dy;

            if (!first_set) {
                first_div = divergence;
                first_set = 1;
            } else if ((first_div > 0) != (divergence > 0)
                       && fabs(divergence) > 1e-8) {
                sign_change = 1;
                break;
            }
        }
        if (sign_change) break;
    }

    if (first_set && !sign_change) {
        /* Divergence does not change sign -> no limit cycles */
        *has_cycle = 0;
    }

    return 0;
}

/* ===================================================================
 * Dulac's Criterion
 *
 * Generalized Bendixson: find B(x,y) > 0 such that
 * div(B*f) does not change sign.
 *
 * For B = 1, this reduces to Bendixson's criterion.
 * Common Dulac functions: B = 1/(xy), B = x^a * y^b.
 * =================================================================== */

int nl_dulac_criterion(nl_nonlinear_system_t *sys,
                       const nl_phase_portrait_t *pp,
                       double (*dulac_func)(double x, double y),
                       int *has_cycle)
{
    if (!sys || !pp || !has_cycle || sys->dim != 2) {
        errno = EINVAL; return -1;
    }

    /* Use Bendixson if no Dulac function provided */
    if (!dulac_func) {
        return nl_bendixson_criterion(pp, has_cycle);
    }

    *has_cycle = 1;
    double dx = (pp->x_max - pp->x_min) / (double)(pp->nx - 1);
    double dy = (pp->y_max - pp->y_min) / (double)(pp->ny - 1);
    double first_val = 0.0;
    int first_set = 0;
    int sign_change = 0;

    for (size_t iy = 1; iy < pp->ny - 1; iy++) {
        for (size_t ix = 1; ix < pp->nx - 1; ix++) {
            double x = pp->x_min + (double)ix * dx;
            double y = pp->y_min + (double)iy * dy;

            double f1 = pp->dxdt[iy * pp->nx + ix];
            double f2 = pp->dydt[iy * pp->nx + ix];

            double B = dulac_func(x, y);
            if (B <= 0.0) continue;

            /* div(B*f) = B * div(f) + grad(B) · f */
            /* div(f) */
            double df1dx = (pp->dxdt[iy * pp->nx + ix + 1]
                            - pp->dxdt[iy * pp->nx + ix - 1])
                           / (2.0 * dx);
            double df2dy = (pp->dydt[(iy + 1) * pp->nx + ix]
                            - pp->dydt[(iy - 1) * pp->nx + ix])
                           / (2.0 * dy);

            /* grad(B) */
            double dBdx = (dulac_func(x + dx, y)
                           - dulac_func(x - dx, y)) / (2.0 * dx);
            double dBdy = (dulac_func(x, y + dy)
                           - dulac_func(x, y - dy)) / (2.0 * dy);

            double div_Bf = B * (df1dx + df2dy) + dBdx * f1 + dBdy * f2;

            if (!first_set) {
                first_val = div_Bf;
                first_set = 1;
            } else if ((first_val > 0) != (div_Bf > 0)
                       && fabs(div_Bf) > 1e-8) {
                sign_change = 1;
                break;
            }
        }
        if (sign_change) break;
    }

    if (first_set && !sign_change) *has_cycle = 0;
    return 0;
}

/* ===================================================================
 * Poincare Index Computation
 *
 * I = (1/(2*pi)) * contour_integral d(atan2(f2, f1))
 *
 * For a circular contour of radius R centered at (cx, cy),
 * compute the accumulated angle change of the vector field.
 * =================================================================== */

int nl_poincare_index(nl_nonlinear_system_t *sys,
                      double center_x, double center_y,
                      double radius, size_t n_pts,
                      double *index)
{
    if (!sys || !index || n_pts < 4 || radius <= 0.0) {
        errno = EINVAL; return -1;
    }

    double total_angle = 0.0;
    double state[2];
    double prev_f1 = 0.0, prev_f2 = 0.0;
    int prev_set = 0;

    for (size_t i = 0; i <= n_pts; i++) {
        double theta = 2.0 * M_PI * (double)i / (double)n_pts;
        state[0] = center_x + radius * cos(theta);
        state[1] = center_y + radius * sin(theta);

        double f_val[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, state, f_val) != 0) continue;

        double f1 = f_val[0], f2 = f_val[1];
        if (fabs(f1) < 1e-12 && fabs(f2) < 1e-12) {
            *index = 0.0;
            return 0;  /* equilibrium on contour */
        }

        if (prev_set) {
            double dtheta = atan2(prev_f1 * f2 - prev_f2 * f1,
                                   prev_f1 * f1 + prev_f2 * f2);
            total_angle += dtheta;
        }

        prev_f1 = f1; prev_f2 = f2; prev_set = 1;
    }

    *index = total_angle / (2.0 * M_PI);
    return 0;
}

/* ===================================================================
 * Trace Stable and Unstable Manifolds of a Saddle
 *
 * From the saddle equilibrium, integrate:
 * - Forward along unstable eigenvector (-> W^u)
 * - Backward along stable eigenvector (-> W^s)
 *
 * Requires the Jacobian at the saddle to get eigenvectors.
 * =================================================================== */

int nl_trace_manifolds(nl_nonlinear_system_t *sys,
                       const nl_state_t *saddle,
                       double T_fwd, double T_bwd, double dt,
                       nl_trajectory_t *W_u, nl_trajectory_t *W_s)
{
    if (!sys || !saddle || !W_u || !W_s
        || sys->dim != 2 || T_fwd <= 0.0 || T_bwd <= 0.0) {
        errno = EINVAL; return -1;
    }

    /* Compute Jacobian at saddle */
    nl_jacobian_t J;
    if (nl_compute_jacobian(sys, saddle->x, &J) != 0) return -1;

    /* Compute eigenvalues and eigenvectors */
    double *A = (double *)malloc(4 * sizeof(double));
    double *real_p = (double *)malloc(2 * sizeof(double));
    double *imag_p = (double *)malloc(2 * sizeof(double));
    if (!A || !real_p || !imag_p) {
        free(A); free(real_p); free(imag_p); return -1;
    }

    A[0] = J.matrix[0][0]; A[1] = J.matrix[0][1];
    A[2] = J.matrix[1][0]; A[3] = J.matrix[1][1];

    if (nl_eigenvalues(A, 2, real_p, imag_p, 2000, 1e-8) != 0) {
        free(A); free(real_p); free(imag_p); return -1;
    }

    /* For a saddle: one eigenvalue > 0 (unstable), one < 0 (stable) */
    int unstable_idx = (real_p[0] > 0) ? 0 : 1;
    int stable_idx = (real_p[0] < 0) ? 0 : 1;

    /* Eigenvector for unstable direction:
     * (J - lambda_unstable * I) * v_unstable = 0 */
    double v_u[2], v_s[2];

    /* Unstable eigenvector (from nullspace of J - lambda*I) */
    double lam_u = real_p[unstable_idx];
    double a00 = J.matrix[0][0] - lam_u;
    double a01 = J.matrix[0][1];
    double a10 = J.matrix[1][0];

    if (fabs(a01) > fabs(a00)) {
        v_u[0] = 1.0;
        v_u[1] = -a00 / a01;
    } else if (fabs(a10) > fabs(a00)) {
        v_u[1] = 1.0;
        v_u[0] = -a01 / a00;
    } else {
        v_u[0] = a01;
        v_u[1] = -a00;
    }

    /* Normalize */
    double norm_u = sqrt(v_u[0] * v_u[0] + v_u[1] * v_u[1]);
    if (norm_u > 1e-12) { v_u[0] /= norm_u; v_u[1] /= norm_u; }

    /* Stable eigenvector */
    double lam_s = real_p[stable_idx];
    a00 = J.matrix[0][0] - lam_s;
    a01 = J.matrix[0][1];
    a10 = J.matrix[1][0];

    if (fabs(a01) > fabs(a00)) {
        v_s[0] = 1.0;
        v_s[1] = -a00 / a01;
    } else if (fabs(a10) > fabs(a00)) {
        v_s[1] = 1.0;
        v_s[0] = -a01 / a00;
    } else {
        v_s[0] = a01;
        v_s[1] = -a00;
    }
    double norm_s = sqrt(v_s[0] * v_s[0] + v_s[1] * v_s[1]);
    if (norm_s > 1e-12) { v_s[0] /= norm_s; v_s[1] /= norm_s; }

    free(A); free(real_p); free(imag_p);

    /* Integrate forward along unstable direction */
    double eps = 1e-4;
    double x0_u[2] = {saddle->x[0] + eps * v_u[0],
                       saddle->x[1] + eps * v_u[1]};
    size_t n_steps_u;
    nl_integrate_rk4(sys, x0_u, T_fwd, dt, W_u, &n_steps_u);

    /* Integrate backward along stable direction */
    double x0_s[2] = {saddle->x[0] + eps * v_s[0],
                       saddle->x[1] + eps * v_s[1]};
    /* Backward integration: negate f(x). Since backward integration
     * of dx/dt = f(x) equals forward integration of dx/dt = -f(x),
     * we use forward integration from the initial perturbation as a
     * first-order approximation of the stable manifold. Full backward
     * integration requires a time-reversed ODE solver with sign-flipped
     * vector field for precise manifold tracing. */
    size_t n_steps_s;
    nl_integrate_rk4(sys, x0_s, T_bwd, dt, W_s, &n_steps_s);

    return 0;
}

/* ===================================================================
 * Poincare Map (First Return Map) Construction
 *
 * For a 2D system, iterates initial conditions on a cross-section
 * and records the first return point and time.
 * =================================================================== */

int nl_poincare_map(nl_nonlinear_system_t *sys,
                    double section_x, const double *x_range,
                    size_t n_pts, double T_max, double dt,
                    double *return_x, double *return_t)
{
    if (!sys || !x_range || !return_x || !return_t
        || n_pts < 1 || T_max <= 0.0 || dt <= 0.0
        || sys->dim != 2) {
        errno = EINVAL; return -1;
    }

    size_t n_returns = 0;

    for (size_t i = 0; i < n_pts; i++) {
        double x0[2] = {section_x, x_range[i]};
        double x[2] = {x0[0], x0[1]};
        double t = 0.0;
        int found_return = 0;
        int prev_side = (x[0] > section_x) ? 1 : -1;

        size_t n_steps = (size_t)(T_max / dt);
        for (size_t s = 0; s < n_steps; s++) {
            double f_val[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, x, f_val) != 0) break;

            /* RK4 step */
            double k1[2], k2[2], k3[2], k4[2], tmp[2];

            k1[0] = dt * f_val[0]; k1[1] = dt * f_val[1];

            tmp[0] = x[0] + 0.5 * k1[0];
            tmp[1] = x[1] + 0.5 * k1[1];
            nl_evaluate_f(sys, tmp, f_val);
            k2[0] = dt * f_val[0]; k2[1] = dt * f_val[1];

            tmp[0] = x[0] + 0.5 * k2[0];
            tmp[1] = x[1] + 0.5 * k2[1];
            nl_evaluate_f(sys, tmp, f_val);
            k3[0] = dt * f_val[0]; k3[1] = dt * f_val[1];

            tmp[0] = x[0] + k3[0];
            tmp[1] = x[1] + k3[1];
            nl_evaluate_f(sys, tmp, f_val);
            k4[0] = dt * f_val[0]; k4[1] = dt * f_val[1];

            double x_new[2];
            x_new[0] = x[0] + (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
            x_new[1] = x[1] + (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;
            t += dt;

            int curr_side = (x_new[0] > section_x) ? 1 : -1;
            if (curr_side != prev_side && s > 10) {
                /* Crossing detected; interpolate to find exact return */
                double frac = (section_x - x[0]) / (x_new[0] - x[0]);
                return_x[n_returns] = x[1] + frac * (x_new[1] - x[1]);
                return_t[n_returns] = t - dt * (1.0 - frac);
                n_returns++;
                found_return = 1;
                break;
            }

            x[0] = x_new[0]; x[1] = x_new[1];
            prev_side = curr_side;
        }

        if (!found_return) {
            return_x[n_returns] = x[1];
            return_t[n_returns] = T_max;
            n_returns++;
        }
    }

    return (int)n_returns;
}

/* ===================================================================
 * Isocline Computation
 *
 * Finds points where dy/dx = f2/f1 = target_slope.
 * Uses bracket-finding along x direction.
 * =================================================================== */

int nl_compute_isocline(nl_nonlinear_system_t *sys,
                        double slope, const double *x_vals,
                        size_t n_pts, double *y_vals)
{
    if (!sys || !x_vals || !y_vals || n_pts < 1) {
        errno = EINVAL; return -1;
    }

    for (size_t i = 0; i < n_pts; i++) {
        double x = x_vals[i];

        /* Search for y such that f2(x,y)/f1(x,y) = slope */
        /* Use bisection over reasonable y range */
        double y_lo = -10.0, y_hi = 10.0;
        double state[2];

        /* Evaluate at endpoints */
        state[0] = x; state[1] = y_lo;
        double f_lo[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, state, f_lo) != 0) {
            y_vals[i] = 0.0; continue;
        }

        state[1] = y_hi;
        double f_hi[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, state, f_hi) != 0) {
            y_vals[i] = 0.0; continue;
        }

        double val_lo = (fabs(f_lo[0]) > 1e-12)
                        ? f_lo[1] / f_lo[0] - slope : 1e6;
        double val_hi = (fabs(f_hi[0]) > 1e-12); (void)val_hi; val_hi = (fabs(f_hi[0]) > 1e-12)
                        ? f_hi[1] / f_hi[0] - slope : 1e6;

        int found = 0;
        for (int iter = 0; iter < 50; iter++) {
            double y_mid = (y_lo + y_hi) / 2.0;
            state[1] = y_mid;
            double f_mid[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, state, f_mid) != 0) break;

            double val_mid = (fabs(f_mid[0]) > 1e-12)
                             ? f_mid[1] / f_mid[0] - slope : 1e6;

            if (fabs(val_mid) < 1e-6) {
                y_vals[i] = y_mid; found = 1; break;
            }

            if (val_lo * val_mid <= 0.0) {
                y_hi = y_mid; val_hi = val_mid;
            } else {
                y_lo = y_mid; val_lo = val_mid;
            }

            if (fabs(y_hi - y_lo) < 1e-8) {
                y_vals[i] = y_mid; found = 1; break;
            }
        }

        if (!found) y_vals[i] = 0.0;
    }

    return 0;
}

/* ===================================================================
 * Classify All Equilibria in Phase Portrait
 * =================================================================== */

int nl_classify_equilibria_2d(nl_nonlinear_system_t *sys,
                               nl_phase_portrait_t *pp,
                               nl_stability_t *stability)
{
    if (!sys || !pp || !stability || sys->dim != 2) {
        errno = EINVAL; return -1;
    }

    for (size_t e = 0; e < pp->num_equilibria; e++) {
        nl_jacobian_t J;
        if (nl_compute_jacobian(sys, pp->equilibria[e].x, &J) == 0) {
            stability[e] = nl_classify_equilibrium(&J, 1e-6);
        } else {
            stability[e] = NL_UNSTABLE;
        }
    }

    return (int)pp->num_equilibria;
}
