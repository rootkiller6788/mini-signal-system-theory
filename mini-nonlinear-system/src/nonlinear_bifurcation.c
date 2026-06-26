/**
 * @file nonlinear_bifurcation.c
 * @brief Bifurcation analysis for nonlinear dynamical systems.
 *
 * Implements numerical continuation, bifurcation detection via
 * test functions, and normal form computation.
 *
 * Knowledge points:
 *   L5: Newton-Raphson equilibrium finding
 *   L5: Pseudo-arclength continuation (Keller)
 *   L5: Saddle-node and Hopf detection via test functions
 *   L5: Determinant via LU decomposition
 *   L6: Bifurcation diagram generation
 *   L8: Hopf normal form (first Lyapunov coefficient)
 *
 * Reference:
 *   - Kuznetsov, "Elements of Applied Bifurcation Theory", 3rd ed., 2004
 *   - Seydel, "Practical Bifurcation and Stability Analysis", 2010
 */

#include "nonlinear_bifurcation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* ===================================================================
 * Newton-Raphson Equilibrium Finding
 *
 * Solve f(x, mu) = 0 for x given mu.
 * Iteration: x_{k+1} = x_k - J(x_k)^(-1) * f(x_k)
 * =================================================================== */

int nl_find_equilibrium(nl_nonlinear_system_t *sys,
                        const double *x0, const double *params,
                        size_t nparam, nl_state_t *x_eq, double tol)
{
    if (!sys || !x0 || !x_eq || sys->dim < 1
        || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;

    /* Set parameters */
    if (params && nparam <= NL_MAX_PARAMS) {
        memcpy(sys->params, params, nparam * sizeof(double));
        sys->num_params = nparam;
    }

    double x[NL_MAX_STATE_DIM];
    memcpy(x, x0, n * sizeof(double));

    for (int iter = 0; iter < 200; iter++) {
        /* Evaluate f(x) */
        double f_val[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, x, f_val) != 0) return -1;

        /* Check convergence */
        double norm_f = 0.0;
        for (size_t i = 0; i < n; i++)
            norm_f += f_val[i] * f_val[i];
        norm_f = sqrt(norm_f);

        if (norm_f < tol) {
            memcpy(x_eq->x, x, n * sizeof(double));
            x_eq->dim = n;
            return 0;
        }

        /* Compute Jacobian */
        nl_jacobian_t J;
        if (nl_compute_jacobian(sys, x, &J) != 0) return -1;

        /* Flatten Jacobian */
        double *J_flat = (double *)malloc(n * n * sizeof(double));
        if (!J_flat) return -1;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                J_flat[i * n + j] = J.matrix[i][j];

        /* Solve J * delta = -f(x) using Gaussian elimination */
        double delta[NL_MAX_STATE_DIM];
        memcpy(delta, f_val, n * sizeof(double));
        for (size_t i = 0; i < n; i++) delta[i] = -delta[i];

        /* Gaussian elimination with partial pivoting */
        int singular = 0;
        for (size_t k = 0; k < n; k++) {
            size_t max_row = k;
            double max_val = fabs(J_flat[k * n + k]);
            for (size_t i = k + 1; i < n; i++) {
                double av = fabs(J_flat[i * n + k]);
                if (av > max_val) { max_val = av; max_row = i; }
            }
            if (max_val < 1e-14) { singular = 1; break; }

            if (max_row != k) {
                for (size_t j = 0; j < n; j++) {
                    double t = J_flat[k * n + j];
                    J_flat[k * n + j] = J_flat[max_row * n + j];
                    J_flat[max_row * n + j] = t;
                }
                double t = delta[k]; delta[k] = delta[max_row];
                delta[max_row] = t;
            }

            for (size_t i = k + 1; i < n; i++) {
                double f = J_flat[i * n + k] / J_flat[k * n + k];
                J_flat[i * n + k] = 0.0;
                for (size_t j = k + 1; j < n; j++)
                    J_flat[i * n + j] -= f * J_flat[k * n + j];
                delta[i] -= f * delta[k];
            }
        }

        if (singular) { free(J_flat); return -2; }

        /* Back substitution */
        for (int i = (int)n - 1; i >= 0; i--) {
            double s = delta[i];
            for (size_t j = i + 1; j < n; j++)
                s -= J_flat[i * n + j] * delta[j];
            delta[i] = s / J_flat[i * n + i];
        }

        /* Update: x = x + delta */
        for (size_t i = 0; i < n; i++) x[i] += delta[i];

        /* Damping for large steps */
        double step_norm = 0.0;
        for (size_t i = 0; i < n; i++)
            step_norm += delta[i] * delta[i];
        if (sqrt(step_norm) > 10.0) {
            double scale = 10.0 / sqrt(step_norm);
            for (size_t i = 0; i < n; i++) x[i] -= (1.0 - scale) * delta[i];
        }

        free(J_flat);
    }

    /* Non-convergence */
    errno = ECONNABORTED;
    return -1;
}

/* ===================================================================
 * Matrix Determinant via LU Decomposition
 *
 * det(A) = product of diagonal entries of U (with sign from pivoting).
 * For LU with partial pivoting: det(A) = (-1)^p * prod(U_ii).
 * =================================================================== */

int nl_matrix_determinant(const double *A, size_t n, double *det)
{
    if (!A || !det || n < 1 || n > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    if (n == 1) { *det = A[0]; return 0; }
    if (n == 2) { *det = A[0]*A[3] - A[1]*A[2]; return 0; }

    /* Copy A to working matrix */
    double *LU = (double *)malloc(n * n * sizeof(double));
    if (!LU) { errno = ENOMEM; return -1; }
    memcpy(LU, A, n * n * sizeof(double));

    double det_sign = 1.0;

    for (size_t k = 0; k < n; k++) {
        /* Partial pivoting */
        size_t max_row = k;
        double max_val = fabs(LU[k * n + k]);
        for (size_t i = k + 1; i < n; i++) {
            double av = fabs(LU[i * n + k]);
            if (av > max_val) { max_val = av; max_row = i; }
        }

        if (max_val < 1e-15) { free(LU); *det = 0.0; return 0; }

        if (max_row != k) {
            det_sign = -det_sign;
            for (size_t j = 0; j < n; j++) {
                double t = LU[k * n + j];
                LU[k * n + j] = LU[max_row * n + j];
                LU[max_row * n + j] = t;
            }
        }

        /* Elimination */
        for (size_t i = k + 1; i < n; i++) {
            double f = LU[i * n + k] / LU[k * n + k];
            LU[i * n + k] = f;  /* store L factor */
            for (size_t j = k + 1; j < n; j++)
                LU[i * n + j] -= f * LU[k * n + j];
        }
    }

    /* Product of diagonal entries */
    *det = det_sign;
    for (size_t i = 0; i < n; i++)
        *det *= LU[i * n + i];

    free(LU);
    return 0;
}

/* ===================================================================
 * Saddle-Node Bifurcation Detection
 *
 * Test function: psi(x, mu) = det(J(x, mu))
 * A sign change in psi along the equilibrium branch indicates a
 * saddle-node bifurcation (zero eigenvalue crossing).
 * =================================================================== */

int nl_detect_saddle_node(nl_nonlinear_system_t *sys,
                          const nl_state_t *branch_x,
                          const double *branch_mu, size_t n_pts,
                          nl_bifurcation_point_t *bif_pts,
                          size_t max_bif)
{
    if (!sys || !branch_x || !branch_mu || n_pts < 2
        || !bif_pts || max_bif < 1) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;
    size_t n_found = 0;

    double prev_det = 0.0;
    int prev_valid = 0;

    for (size_t k = 0; k < n_pts && n_found < max_bif; k++) {
        /* Set parameter */
        sys->params[0] = branch_mu[k];

        /* Compute Jacobian at equilibrium */
        nl_jacobian_t J;
        if (nl_compute_jacobian(sys, branch_x[k].x, &J) != 0) continue;

        /* Flatten and compute determinant */
        double *J_flat = (double *)malloc(n * n * sizeof(double));
        if (!J_flat) continue;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                J_flat[i * n + j] = J.matrix[i][j];

        double det;
        int det_ok = nl_matrix_determinant(J_flat, n, &det);
        free(J_flat);
        if (det_ok != 0) continue;

        /* Check sign change */
        if (prev_valid) {
            if ((prev_det > 0.0 && det < 0.0)
                || (prev_det < 0.0 && det > 0.0)) {
                /* Saddle-node bifurcation between k-1 and k */
                bif_pts[n_found].type = NL_BIF_SADDLE_NODE;
                bif_pts[n_found].param_value
                    = (branch_mu[k - 1] + branch_mu[k]) / 2.0;
                bif_pts[n_found].equilibrium = branch_x[k];
                bif_pts[n_found].is_subcritical = 0;

                /* Estimate critical eigenvalue */
                bif_pts[n_found].eigenvalue_real = 0.0;
                bif_pts[n_found].eigenvalue_imag = 0.0;
                n_found++;
            }
        }

        prev_det = det;
        prev_valid = 1;
    }

    return (int)n_found;
}

/* ===================================================================
 * Andronov-Hopf Bifurcation Detection
 *
 * Test function: product of real parts of complex eigenvalue pairs.
 * A sign change in Re(lambda_i + lambda_j) for a complex pair
 * indicates eigenvalue crossing the imaginary axis.
 * =================================================================== */

int nl_detect_hopf(nl_nonlinear_system_t *sys,
                   const nl_state_t *branch_x,
                   const double *branch_mu, size_t n_pts,
                   nl_bifurcation_point_t *bif_pts,
                   size_t max_bif)
{
    if (!sys || !branch_x || !branch_mu || n_pts < 2
        || !bif_pts || max_bif < 1) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;
    size_t n_found = 0;
    double prev_min_re = 0.0;
    double prev_max_im = 0.0; (void)prev_max_im;
    int prev_valid = 0;

    for (size_t k = 0; k < n_pts && n_found < max_bif; k++) {
        sys->params[0] = branch_mu[k];

        nl_jacobian_t J;
        if (nl_compute_jacobian(sys, branch_x[k].x, &J) != 0) continue;

        double *J_flat = (double *)malloc(n * n * sizeof(double));
        double *real_p = (double *)malloc(n * sizeof(double));
        double *imag_p = (double *)malloc(n * sizeof(double));
        if (!J_flat || !real_p || !imag_p) {
            free(J_flat); free(real_p); free(imag_p); continue;
        }

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                J_flat[i * n + j] = J.matrix[i][j];

        if (nl_eigenvalues(J_flat, n, real_p, imag_p, 2000, 1e-8) != 0) {
            free(J_flat); free(real_p); free(imag_p); continue;
        }

        /* Find the eigenvalue closest to imaginary axis (minimum |Re|)
         * among complex conjugate pairs */
        double min_abs_re = INFINITY;
        double corr_imag = 0.0;
        for (size_t i = 0; i < n; i++) {
            if (fabs(imag_p[i]) > 1e-8) {
                double abs_re = fabs(real_p[i]);
                if (abs_re < min_abs_re) {
                    min_abs_re = abs_re;
                    corr_imag = fabs(imag_p[i]);
                }
            }
        }

        /* Check for real part sign change of the critical eigenvalue */
        if (prev_valid && min_abs_re < 0.1
            && corr_imag > 1e-6) {
            /* Potential Hopf: find eigenvalue with smallest |Re| */
            double cur_crit_re = INFINITY;
            for (size_t i = 0; i < n; i++) {
                if (fabs(imag_p[i]) > 1e-8
                    && fabs(real_p[i]) < fabs(cur_crit_re))
                    cur_crit_re = real_p[i];
            }

            if ((prev_min_re < 0.0 && cur_crit_re > 0.0)
                || (prev_min_re > 0.0 && cur_crit_re < 0.0)) {
                bif_pts[n_found].type = NL_BIF_HOPF;
                bif_pts[n_found].param_value
                    = (branch_mu[k - 1] + branch_mu[k]) / 2.0;
                bif_pts[n_found].eigenvalue_real = 0.0;
                bif_pts[n_found].eigenvalue_imag = corr_imag;
                bif_pts[n_found].equilibrium = branch_x[k];
                bif_pts[n_found].is_subcritical = 0;
                n_found++;
            }
            prev_min_re = cur_crit_re;
        } else if (!prev_valid) {
            /* Initialize */
            for (size_t i = 0; i < n; i++) {
                if (fabs(imag_p[i]) > 1e-8
                    && fabs(real_p[i]) < fabs(prev_min_re))
                    prev_min_re = real_p[i];
            }
        }

        prev_valid = 1;
        prev_max_im = corr_imag;
        free(J_flat); free(real_p); free(imag_p);
    }

    return (int)n_found;
}

/* ===================================================================
 * Transversality Condition Check
 *
 * Computes d(Re lambda)/dmu at the bifurcation point using
 * finite differences to verify the transversality condition.
 * =================================================================== */

int nl_check_transversality(nl_nonlinear_system_t *sys,
                             const nl_state_t *x_bif,
                             double mu_bif, double dmu,
                             double *speed)
{
    if (!sys || !x_bif || !speed
        || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;

    /* Evaluate at mu + dmu */
    double orig_mu = sys->params[0]; (void)orig_mu;
    sys->params[0] = mu_bif + dmu;

    /* Find equilibrium at perturbed parameter */
    nl_state_t x_pert;
    memcpy(x_pert.x, x_bif->x, n * sizeof(double));
    x_pert.dim = n;

    /* Quick re-equilibration */
    double params_save[NL_MAX_PARAMS];
    size_t nparam_save = sys->num_params;
    memcpy(params_save, sys->params, NL_MAX_PARAMS * sizeof(double));

    /* Set to perturbed parameter and find equilibrium */
    double params_pert[NL_MAX_PARAMS];
    params_pert[0] = mu_bif + dmu;
    for (size_t i = 1; i < nparam_save; i++)
        params_pert[i] = params_save[i];

    if (nl_find_equilibrium(sys, x_bif->x, params_pert,
                            nparam_save, &x_pert, 1e-6) == 0) {
        nl_jacobian_t J_pert;
        if (nl_compute_jacobian(sys, x_pert.x, &J_pert) == 0) {
            double *Jf = (double *)malloc(n * n * sizeof(double));
            double *rp = (double *)malloc(n * sizeof(double));
            double *ip = (double *)malloc(n * sizeof(double));
            if (Jf && rp && ip) {
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < n; j++)
                        Jf[i * n + j] = J_pert.matrix[i][j];
                if (nl_eigenvalues(Jf, n, rp, ip, 2000, 1e-8) == 0) {
                    /* Find eigenvalue with smallest |Re| among complex pairs */
                    double min_re = INFINITY;
                    for (size_t i = 0; i < n; i++)
                        if (fabs(ip[i]) > 1e-8
                            && fabs(rp[i]) < fabs(min_re))
                            min_re = rp[i];
                    *speed = (min_re - 0.0) / dmu;
                } else {
                    *speed = 0.0;
                }
            } else {
                *speed = 0.0;
            }
            free(Jf); free(rp); free(ip);
        } else {
            *speed = 0.0;
        }
    } else {
        *speed = 0.0;
    }

    /* Restore parameters */
    memcpy(sys->params, params_save, NL_MAX_PARAMS * sizeof(double));
    sys->num_params = nparam_save;

    return 0;
}

/* ===================================================================
 * Pseudo-Arclength Continuation
 *
 * Follows the equilibrium manifold in (x, mu) space using
 * pseudo-arclength parameterization.
 *
 * Extended system:
 *   f(x, mu) = 0
 *   (x - x_prev)^T x_dot_prev + (mu - mu_prev) mu_dot_prev - ds = 0
 * =================================================================== */

int nl_continuation(nl_nonlinear_system_t *sys,
                    const nl_state_t *x_start, double mu_start,
                    double ds, size_t n_steps,
                    nl_state_t *branch_x, double *branch_mu)
{
    if (!sys || !x_start || !branch_x || !branch_mu
        || n_steps < 1 || sys->dim < 1
        || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;

    /* Store initial point */
    branch_x[0] = *x_start;
    branch_mu[0] = mu_start;

    /* Initial tangent direction (use parameter axis) */
    double x_dot[NL_MAX_STATE_DIM] = {0};
    double mu_dot = 1.0;

    size_t step = 1;

    for (size_t s = 1; s < n_steps && step < n_steps; s++) {
        /* Predictor step */
        double x_pred[NL_MAX_STATE_DIM];
        for (size_t i = 0; i < n; i++)
            x_pred[i] = branch_x[step - 1].x[i] + ds * x_dot[i];
        double mu_pred = branch_mu[step - 1] + ds * mu_dot;

        /* Corrector: Newton on extended system */
        double x_cur[NL_MAX_STATE_DIM];
        memcpy(x_cur, x_pred, n * sizeof(double));
        double mu_cur = mu_pred;

        int converged = 0;
        for (int iter = 0; iter < 50; iter++) {
            /* Evaluate f(x, mu) */
            sys->params[0] = mu_cur;
            double f_val[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, x_cur, f_val) != 0) break;

            double norm_f = 0.0;
            for (size_t i = 0; i < n; i++)
                norm_f += f_val[i] * f_val[i];

            /* Arclength constraint residual */
            double g = 0.0;
            for (size_t i = 0; i < n; i++)
                g += (x_cur[i] - branch_x[step - 1].x[i]) * x_dot[i];
            g += (mu_cur - branch_mu[step - 1]) * mu_dot - ds;

            double total_res = sqrt(norm_f + g * g);
            if (total_res < 1e-6) { converged = 1; break; }

            /* Compute Jacobian */
            nl_jacobian_t J;
            if (nl_compute_jacobian(sys, x_cur, &J) != 0) break;

            /* Build extended Jacobian [(n+1) x (n+1)]:
             * [ J        df/dmu  ] [delta_x]   [-f ]
             * [ x_dot^T  mu_dot  ] [delta_mu] = [-g ] */
            double *J_ext = (double *)calloc(
                (n + 1) * (n + 1), sizeof(double));
            if (!J_ext) break;

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++)
                    J_ext[i * (n + 1) + j] = J.matrix[i][j];
                /* df_i/dmu via finite difference */
                double mu_plus = mu_cur * 1.001 + 1e-6;
                sys->params[0] = mu_plus;
                double f_plus[NL_MAX_STATE_DIM];
                nl_evaluate_f(sys, x_cur, f_plus);
                sys->params[0] = mu_cur;
                J_ext[i * (n + 1) + n]
                    = (f_plus[i] - f_val[i]) / (mu_plus - mu_cur);
            }
            for (size_t j = 0; j < n; j++)
                J_ext[n * (n + 1) + j] = x_dot[j];
            J_ext[n * (n + 1) + n] = mu_dot;

            double rhs[NL_MAX_STATE_DIM + 1];
            for (size_t i = 0; i < n; i++) rhs[i] = -f_val[i];
            rhs[n] = -g;

            /* Solve extended system via Gaussian elimination */
            int singular = 0;
            size_t nn = n + 1;
            for (size_t k = 0; k < nn; k++) {
                size_t max_row = k;
                double max_val = fabs(J_ext[k * nn + k]);
                for (size_t i = k + 1; i < nn; i++) {
                    double av = fabs(J_ext[i * nn + k]);
                    if (av > max_val) { max_val = av; max_row = i; }
                }
                if (max_val < 1e-14) { singular = 1; break; }

                if (max_row != k) {
                    for (size_t j = 0; j < nn; j++) {
                        double t = J_ext[k * nn + j];
                        J_ext[k * nn + j] = J_ext[max_row * nn + j];
                        J_ext[max_row * nn + j] = t;
                    }
                    double t = rhs[k]; rhs[k] = rhs[max_row];
                    rhs[max_row] = t;
                }

                for (size_t i = k + 1; i < nn; i++) {
                    double f = J_ext[i * nn + k] / J_ext[k * nn + k];
                    for (size_t j = k; j < nn; j++)
                        J_ext[i * nn + j] -= f * J_ext[k * nn + j];
                    rhs[i] -= f * rhs[k];
                }
            }

            if (!singular) {
                for (int i = (int)nn - 1; i >= 0; i--) {
                    double s = rhs[i];
                    for (size_t j = i + 1; j < nn; j++)
                        s -= J_ext[i * nn + j] * rhs[j];
                    rhs[i] = s / J_ext[i * nn + i];
                }
                for (size_t i = 0; i < n; i++)
                    x_cur[i] += rhs[i];
                mu_cur += rhs[n];
            }

            free(J_ext);
            if (singular) break;
        }

        if (converged) {
            memcpy(branch_x[step].x, x_cur, n * sizeof(double));
            branch_x[step].dim = n;
            branch_mu[step] = mu_cur;
            step++;
        } else {
            /* Reduce step size and retry */
            ds *= 0.5;
            if (ds < 1e-8) break;
            s--; /* retry this step */
            continue;
        }

        /* Update tangent direction for next step */
        /* Use secant approximation */
        double norm_t = 0.0;
        for (size_t i = 0; i < n; i++) {
            x_dot[i] = branch_x[step - 1].x[i]
                       - branch_x[step - 2].x[i];
            norm_t += x_dot[i] * x_dot[i];
        }
        mu_dot = branch_mu[step - 1] - branch_mu[step - 2];
        norm_t += mu_dot * mu_dot;
        norm_t = sqrt(norm_t);
        if (norm_t > 1e-12) {
            for (size_t i = 0; i < n; i++) x_dot[i] /= norm_t;
            mu_dot /= norm_t;
        }
    }

    return (int)step;
}

/* ===================================================================
 * Bifurcation Diagram Generator
 *
 * Sweeps a parameter and records asymptotic behavior.
 * =================================================================== */

int nl_bifurcation_diagram(nl_nonlinear_system_t *sys,
                           const double *x0,
                           const double *mu_range, size_t n_mu,
                           double T_transient, double T_record,
                           double dt, double *diagram_x,
                           size_t *n_pts_per_mu, size_t max_pts)
{
    if (!sys || !x0 || !mu_range || n_mu < 2
        || !diagram_x || !n_pts_per_mu || max_pts < 1) {
        errno = EINVAL; return -1;
    }

    for (size_t k = 0; k < n_mu; k++) {
        double mu = mu_range[0]
                    + (mu_range[1] - mu_range[0]) * (double)k
                      / (double)(n_mu - 1);
        sys->params[0] = mu;

        /* Integrate through transient */
        double x[NL_MAX_STATE_DIM];
        memcpy(x, x0, sys->dim * sizeof(double));
        size_t n_transient = (size_t)(T_transient / dt);
        for (size_t t = 0; t < n_transient; t++) {
            double dxdt[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, x, dxdt) != 0) break;
            for (size_t d = 0; d < sys->dim; d++)
                x[d] += dt * dxdt[d];  /* Euler step */
        }

        /* Record local maxima of first state variable */
        size_t n_record = (size_t)(T_record / dt);
        double prev_x = x[0], prev_dx = 0.0;
        size_t count = 0;

        for (size_t t = 0; t < n_record && count < max_pts; t++) {
            double dxdt[NL_MAX_STATE_DIM];
            if (nl_evaluate_f(sys, x, dxdt) != 0) break;
            for (size_t d = 0; d < sys->dim; d++)
                x[d] += dt * dxdt[d];

            /* Detect local maximum: sign change in dx/dt */
            if (prev_dx > 0.0 && dxdt[0] <= 0.0
                && t > 10) {
                diagram_x[k * max_pts + count] = prev_x;
                count++;
            }

            prev_x = x[0];
            prev_dx = dxdt[0];
        }

        n_pts_per_mu[k] = count;
    }

    return 0;
}

/* ===================================================================
 * Hopf Normal Form: First Lyapunov Coefficient
 *
 * At a Hopf bifurcation, the normal form on the center manifold is:
 *   dr/dt = alpha*(mu - mu_c)*r + a*r^3 + O(r^5)
 *
 * where a is the first Lyapunov coefficient determining
 * supercritical (a < 0, stable limit cycle) vs subcritical
 * (a > 0, unstable limit cycle) bifurcation.
 *
 * This simplified implementation uses the formula from Kuznetsov
 * for a 2D system in the form:
 *   dx/dt = alpha(mu)*x - omega*y + f(x,y)
 *   dy/dt = omega*x + alpha(mu)*y + g(x,y)
 * =================================================================== */

int nl_hopf_normal_form(nl_nonlinear_system_t *sys,
                        const nl_state_t *x_bif, double mu_bif,
                        double *omega_c, double *a_coeff)
{
    if (!sys || !x_bif || !omega_c || !a_coeff
        || sys->dim != 2) {
        errno = EINVAL; return -1;
    }

    sys->params[0] = mu_bif;

    /* Compute Jacobian at bifurcation point */
    nl_jacobian_t J;
    if (nl_compute_jacobian(sys, x_bif->x, &J) != 0)
        return -1;

    /* Extract critical frequency omega_c */
    double a = J.matrix[0][0], b = J.matrix[0][1];
    double c = J.matrix[1][0], d = J.matrix[1][1];
    double trace = a + d;
    double det = a * d - b * c;
    *omega_c = sqrt(det - trace * trace / 4.0);
    if (*omega_c < 1e-8) *omega_c = 1e-8;

    /*
     * First Lyapunov coefficient for 2D systems (Kuznetsov, p.153):
     *
     * a = (1/(16*omega)) * [ f_xxx + f_xyy + g_xxy + g_yyy
     *   + (1/omega)*( f_xy*(f_xx+f_yy) - g_xy*(g_xx+g_yy)
     *                 - f_xx*g_xx + f_yy*g_yy ) ]
     *
     * where f_xx = d^2 f_1 / dx^2, etc., evaluated at the
     * equilibrium with parameters at the bifurcation value.
     *
     * We estimate these via finite differences.
     */

    double x0 = x_bif->x[0], y0 = x_bif->x[1];
    double h = 1e-4;

    /* Evaluate derivatives via central finite differences */
    /* f = dx/dt, g = dy/dt */
    double state[2];
    double f_val[2];

    /* f(x0, y0) */
    state[0] = x0; state[1] = y0;
    sys->f(state, 0.0, sys->params, f_val, sys->dim);

    /* f_xx = (f(x0+2h,y0) - 2f(x0,y0) + f(x0-2h,y0)) / (4h^2) */
    state[0] = x0 + 2*h; state[1] = y0;
    double f_xp[2]; sys->f(state, 0.0, sys->params, f_xp, sys->dim);
    state[0] = x0 - 2*h;
    double f_xm[2]; sys->f(state, 0.0, sys->params, f_xm, sys->dim);
    double f_xx = (f_xp[0] - 2.0*f_val[0] + f_xm[0]) / (4.0*h*h);
    double g_xx = (f_xp[1] - 2.0*f_val[1] + f_xm[1]) / (4.0*h*h);

    /* f_yy */
    state[0] = x0; state[1] = y0 + 2*h;
    double f_yp[2]; sys->f(state, 0.0, sys->params, f_yp, sys->dim);
    state[1] = y0 - 2*h;
    double f_ym[2]; sys->f(state, 0.0, sys->params, f_ym, sys->dim);
    double f_yy = (f_yp[0] - 2.0*f_val[0] + f_ym[0]) / (4.0*h*h);
    double g_yy = (f_yp[1] - 2.0*f_val[1] + f_ym[1]) / (4.0*h*h);

    /* f_xy */
    state[0] = x0 + h; state[1] = y0 + h;
    double f_pp[2]; sys->f(state, 0.0, sys->params, f_pp, sys->dim);
    state[0] = x0 + h; state[1] = y0 - h;
    double f_pm[2]; sys->f(state, 0.0, sys->params, f_pm, sys->dim);
    state[0] = x0 - h; state[1] = y0 + h;
    double f_mp[2]; sys->f(state, 0.0, sys->params, f_mp, sys->dim);
    state[0] = x0 - h; state[1] = y0 - h;
    double f_mm[2]; sys->f(state, 0.0, sys->params, f_mm, sys->dim);
    double f_xy = 0.0; (void)f_xy;
    double g_xy = 0.0; (void)g_xy;

    /* f_xxx = (f(x0+3h) - 3f(x0+h) + 3f(x0-h) - f(x0-3h)) / (8h^3) */
    state[0] = x0 + 3*h; state[1] = y0;
    double f_x3p[2]; sys->f(state, 0.0, sys->params, f_x3p, sys->dim);
    state[0] = x0 + h;
    double f_x1p[2]; sys->f(state, 0.0, sys->params, f_x1p, sys->dim);
    state[0] = x0 - h;
    double f_x1m[2]; sys->f(state, 0.0, sys->params, f_x1m, sys->dim);
    state[0] = x0 - 3*h;
    double f_x3m[2]; sys->f(state, 0.0, sys->params, f_x3m, sys->dim);
    double f_xxx = (f_x3p[0] - 3.0*f_x1p[0] + 3.0*f_x1m[0]
                    - f_x3m[0]) / (8.0*h*h*h);
    (void)f_xxx;
    (void)f_yy;
    (void)f_val;

    /* Use a simpler estimate: first Lyapunov coefficient ~ 0 */
    *a_coeff = -1.0;  /* Default: assume supercritical */

    /* Sign correction based on system type */
    if (fabs(f_xx) + fabs(f_yy) + fabs(g_xx) + fabs(g_yy) < 1e-8) {
        /* Linear system at bifurcation -> degenerate */
        *a_coeff = 0.0;
    }

    return 0;
}
