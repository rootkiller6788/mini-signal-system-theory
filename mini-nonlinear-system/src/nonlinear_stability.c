/**
 * @file nonlinear_stability.c
 * @brief Lyapunov stability analysis implementations.
 *
 * Implements Lyapunov direct and indirect methods, eigenvalue
 * computation via QR algorithm, algebraic Lyapunov equation solver,
 * the circle criterion, Popov criterion, and related utilities.
 *
 * Knowledge points:
 *   L3: Jacobian via finite differences, QR eigensolver
 *   L4: Lyapunov direct/indirect methods, Circle/Popov criteria
 *   L5: Gershgorin for positive definiteness check
 */

#include "nonlinear_stability.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* Vector and Matrix Helpers */

/** Gaussian elimination with partial pivoting.
 *  A is row-major n*n, b is length n, modified in place.
 *  Solution returned in b. Returns 0 on success, -1 if singular. */
static int gauss_solve(double *A, double *b, size_t n) {
    for (size_t k = 0; k < n; k++) {
        size_t max_row = k;
        double max_val = fabs(A[k * n + k]);
        for (size_t i = k + 1; i < n; i++) {
            double av = fabs(A[i * n + k]);
            if (av > max_val) { max_val = av; max_row = i; }
        }
        if (max_val < 1e-15) return -1;
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k]; b[k] = b[max_row]; b[max_row] = tmp;
        }
        for (size_t i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / A[k * n + k];
            A[i * n + k] = 0.0;
            for (size_t j = k + 1; j < n; j++)
                A[i * n + j] -= factor * A[k * n + j];
            b[i] -= factor * b[k];
        }
    }
    for (int i = (int)n - 1; i >= 0; i--) {
        double sum = b[i];
        for (size_t j = i + 1; j < n; j++)
            sum -= A[i * n + j] * b[j];
        if (fabs(A[i * n + i]) < 1e-15) return -1;
        b[i] = sum / A[i * n + i];
    }
    return 0;
}

/* System Evaluation */

int nl_evaluate_f(nl_nonlinear_system_t *sys,
                  const double *x, double *dxdt) {
    if (!sys || !x || !dxdt || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }
    if (!sys->f) { errno = ENOSYS; return -1; }
    return sys->f(x, 0.0, sys->params, dxdt, sys->dim);
}

/* Quadratic Lyapunov Function: V(x) = x^T P x */

double nl_quadratic_lyapunov(const double *P, const double *x, size_t dim) {
    if (!P || !x || dim < 1 || dim > NL_MAX_STATE_DIM) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < dim; j++)
            row_sum += P[i * dim + j] * x[j];
        sum += x[i] * row_sum;
    }
    return sum;
}

double nl_quadratic_lyapunov_dot(const double *P, const double *x,
                                  nl_nonlinear_system_t *sys,
                                  size_t dim) {
    if (!P || !x || !sys || dim < 1 || dim > NL_MAX_STATE_DIM) return 0.0;
    double f_x[NL_MAX_STATE_DIM];
    if (nl_evaluate_f(sys, x, f_x) != 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < dim; j++)
            row_sum += P[i * dim + j] * f_x[j];
        sum += x[i] * row_sum;
    }
    return 2.0 * sum;
}

/* Jacobian via Central Finite Differences */

int nl_compute_jacobian(nl_nonlinear_system_t *sys,
                        const double *x, nl_jacobian_t *J) {
    if (!sys || !x || !J || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }
    size_t n = sys->dim;
    J->dim = n;
    if (sys->jacobian)
        return sys->jacobian(x, 0.0, sys->params, (double *)J->matrix, n);
    double h_base = 1.5e-8;
    double x_pert[NL_MAX_STATE_DIM];
    memcpy(x_pert, x, n * sizeof(double));
    for (size_t j = 0; j < n; j++) {
        double h = h_base * fmax(1.0, fabs(x[j]));
        if (h < 1e-12) h = 1e-8;
        x_pert[j] = x[j] + h;
        double f_plus[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, x_pert, f_plus) != 0) {
            x_pert[j] = x[j]; return -1;
        }
        x_pert[j] = x[j] - h;
        double f_minus[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, x_pert, f_minus) != 0) {
            x_pert[j] = x[j]; return -1;
        }
        x_pert[j] = x[j];
        for (size_t i = 0; i < n; i++)
            J->matrix[i][j] = (f_plus[i] - f_minus[i]) / (2.0 * h);
    }
    return 0;
}

/* Francis QR Double-Shift Eigensolver */

static void to_hessenberg(double *A, size_t n) {
    for (size_t k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++)
            sigma += A[i * n + k] * A[i * n + k];
        if (sigma < 1e-30) continue;
        double alpha = sqrt(sigma);
        if (A[(k + 1) * n + k] > 0) alpha = -alpha;
        double r = sqrt(2.0 * (sigma - A[(k + 1) * n + k] * alpha));
        double v[NL_MAX_STATE_DIM] = {0};
        v[k + 1] = (A[(k + 1) * n + k] - alpha) / r;
        for (size_t i = k + 2; i < n; i++)
            v[i] = A[i * n + k] / r;
        for (size_t j = k; j < n; j++) {
            double dot = 0.0;
            for (size_t i = k + 1; i < n; i++)
                dot += v[i] * A[i * n + j];
            for (size_t i = k + 1; i < n; i++)
                A[i * n + j] -= 2.0 * v[i] * dot;
        }
        for (size_t i = 0; i < n; i++) {
            double dot = 0.0;
            for (size_t j = k + 1; j < n; j++)
                dot += A[i * n + j] * v[j];
            for (size_t j = k + 1; j < n; j++)
                A[i * n + j] -= 2.0 * v[j] * dot;
        }
    }
}

static void extract_eigenvalues(const double *T, size_t n,
                                 double *real_part, double *imag_part) {
    size_t i = 0;
    while (i < n) {
        if (i == n - 1 || fabs(T[(i + 1) * n + i]) < 1e-12) {
            real_part[i] = T[i * n + i];
            imag_part[i] = 0.0;
            i++;
        } else {
            double a = T[i * n + i], b = T[i * n + i + 1];
            double c = T[(i + 1) * n + i], d = T[(i + 1) * n + i + 1];
            double trace = a + d, det = a * d - b * c;
            double disc = trace * trace - 4.0 * det;
            if (disc >= 0.0) {
                double sd = sqrt(disc);
                real_part[i] = (trace + sd) / 2.0;
                imag_part[i] = 0.0;
                real_part[i + 1] = (trace - sd) / 2.0;
                imag_part[i + 1] = 0.0;
            } else {
                real_part[i] = trace / 2.0;
                imag_part[i] = sqrt(-disc) / 2.0;
                real_part[i + 1] = trace / 2.0;
                imag_part[i + 1] = -sqrt(-disc) / 2.0;
            }
            i += 2;
        }
    }
}

static void qr_double_shift(double *H, size_t n, size_t m) {
    double a = H[(m-2)*n+(m-2)], b = H[(m-2)*n+(m-1)];
    double c = H[(m-1)*n+(m-2)], d = H[(m-1)*n+(m-1)];
    double trace = a + d, det = a * d - b * c;
    double x = H[0]*H[0] + H[n]*H[n+1] - trace*H[0] + det;
    double y = H[n] * (H[0] + H[n+1] - trace);
    double z = H[n] * H[2*n+1];
    for (size_t k = 0; k < m - 2; k++) {
        double r = sqrt(x*x + y*y + z*z);
        if (r < 1e-15) { x = 1.0; y = 0.0; z = 0.0; r = 1.0; }
        double v1 = x + copysign(r, x);
        double v2 = y, v3 = z;
        double nrm = sqrt(v1*v1 + v2*v2 + v3*v3);
        v1 /= nrm; v2 /= nrm; v3 /= nrm;
        size_t start_col = (k > 0) ? (k - 1) : 0;
        size_t end_col = (k < m - 3) ? (k + 3) : n;
        for (size_t j = start_col; j < end_col; j++) {
            double dot = v1*H[k*n+j] + v2*H[(k+1)*n+j];
            if (k + 2 < n) dot += v3*H[(k+2)*n+j];
            H[k*n+j] -= 2.0*v1*dot;
            H[(k+1)*n+j] -= 2.0*v2*dot;
            if (k + 2 < n) H[(k+2)*n+j] -= 2.0*v3*dot;
        }
        size_t row_end = (k + 3 < m) ? (k + 3) : m;
        for (size_t i = 0; i < row_end; i++) {
            double dot = v1*H[i*n+k] + v2*H[i*n+k+1];
            if (k + 2 < n) dot += v3*H[i*n+k+2];
            H[i*n+k] -= 2.0*v1*dot;
            H[i*n+k+1] -= 2.0*v2*dot;
            if (k + 2 < n) H[i*n+k+2] -= 2.0*v3*dot;
        }
        if (k < m - 3) {
            x = H[(k+1)*n+k];
            y = H[(k+2)*n+k];
            z = H[(k+3)*n+k];
        }
    }
}

int nl_eigenvalues(double *A, size_t n, double *real_part,
                   double *imag_part, int max_iter, double tol) {
    if (!A || !real_part || !imag_part || n < 1 || n > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    /* Direct formula for 2x2: λ = (τ ± sqrt(τ^2 - 4Δ)) / 2 */
    if (n == 2) {
        double trace = A[0] + A[3];
        double det = A[0] * A[3] - A[1] * A[2];
        double disc = trace * trace - 4.0 * det;
        if (disc >= 0.0) {
            double sd = sqrt(disc);
            real_part[0] = (trace + sd) / 2.0;
            imag_part[0] = 0.0;
            real_part[1] = (trace - sd) / 2.0;
            imag_part[1] = 0.0;
        } else {
            real_part[0] = trace / 2.0;
            imag_part[0] = sqrt(-disc) / 2.0;
            real_part[1] = trace / 2.0;
            imag_part[1] = -sqrt(-disc) / 2.0;
        }
        return 0;
    }

    size_t n2 = n * n;
    double *H = (double *)malloc(n2 * sizeof(double));
    if (!H) { errno = ENOMEM; return -1; }
    memcpy(H, A, n2 * sizeof(double));
    to_hessenberg(H, n);
    int converged = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        size_t m = n;
        while (m > 2) {
            if (fabs(H[(m-1)*n+(m-2)]) < tol) m--;
            else if (fabs(H[(m-2)*n+(m-3)]) < tol) m--;
            else break;
        }
        /* For m==2: 2x2 block is converged (real Schur form keeps
         * subdiagonal nonzero for complex eigenvalues) */
        if (m <= 2) {
            converged = 1; break;
        }
        qr_double_shift(H, n, m);
    }
    extract_eigenvalues(H, n, real_part, imag_part);
    free(H);
    return converged ? 0 : -1;
}

/* Stability Classification */

nl_stability_t nl_classify_equilibrium(const nl_jacobian_t *J, double tol) {
    if (!J || J->dim < 1 || J->dim > NL_MAX_STATE_DIM) return NL_UNSTABLE;
    size_t n = J->dim;
    double *A_flat = (double *)malloc(n * n * sizeof(double));
    double *real_parts = (double *)malloc(n * sizeof(double));
    double *imag_parts = (double *)malloc(n * sizeof(double));
    if (!A_flat || !real_parts || !imag_parts) {
        free(A_flat); free(real_parts); free(imag_parts);
        return NL_UNSTABLE;
    }
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            A_flat[i * n + j] = J->matrix[i][j];
    if (nl_eigenvalues(A_flat, n, real_parts, imag_parts, 2000, tol) != 0) {
        free(A_flat); free(real_parts); free(imag_parts);
        return NL_UNSTABLE;
    }
    int all_neg = 1, any_pos = 0, any_zero_real = 0, has_complex = 0;
    int n_pos = 0, n_neg = 0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(imag_parts[i]) > tol) has_complex = 1;
        if (real_parts[i] > tol) { any_pos = 1; all_neg = 0; n_pos++; }
        else if (real_parts[i] < -tol) n_neg++;
        else { any_zero_real = 1; all_neg = 0; }
    }
    nl_stability_t result;
    if (all_neg)
        result = has_complex ? NL_STABLE_FOCUS : NL_STABLE_NODE;
    else if (any_zero_real && !any_pos)
        result = NL_NONHYPERBOLIC;
    else if (any_pos && n_neg > 0)
        result = NL_SADDLE;
    else if (any_pos)
        result = has_complex ? NL_UNSTABLE_FOCUS : NL_UNSTABLE_NODE;
    else
        result = NL_UNSTABLE;
    free(A_flat); free(real_parts); free(imag_parts);
    return result;
}

/* Lyapunov Equation: A^T P + P A = -Q via Kronecker form */

int nl_solve_lyapunov_equation(const double *A, const double *Q,
                                double *P, size_t n) {
    if (!A || !Q || !P || n < 1 || n > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }
    size_t n2 = n * n, n4 = n2 * n2;
    double *M = (double *)calloc(n4, sizeof(double));
    double *b = (double *)malloc(n2 * sizeof(double));
    if (!M || !b) { free(M); free(b); errno = ENOMEM; return -1; }
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t row = i * n + j;
            b[row] = -Q[i * n + j];
            for (size_t k = 0; k < n; k++) {
                for (size_t l = 0; l < n; l++) {
                    size_t col = k * n + l;
                    double term = 0.0;
                    if (i == k) term += A[l * n + j];
                    if (j == l) term += A[k * n + i];
                    M[row * n2 + col] = term;
                }
            }
        }
    }
    memcpy(P, b, n2 * sizeof(double));
    int ret = gauss_solve(M, P, n2);
    free(M); free(b);
    return ret;
}

/* Lyapunov Direct Method */

nl_stability_t nl_direct_method(nl_nonlinear_system_t *sys,
                                 const nl_state_t *equilibrium,
                                 double *P, double *radius) {
    if (!sys || !equilibrium || !P || !radius
        || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM)
        return NL_UNSTABLE;
    size_t n = sys->dim;
    nl_jacobian_t J;
    if (nl_compute_jacobian(sys, equilibrium->x, &J) != 0) return NL_UNSTABLE;
    double *A_flat = (double *)malloc(n * n * sizeof(double));
    double *real_p = (double *)malloc(n * sizeof(double));
    double *imag_p = (double *)malloc(n * sizeof(double));
    if (!A_flat || !real_p || !imag_p) {
        free(A_flat); free(real_p); free(imag_p); return NL_UNSTABLE;
    }
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            A_flat[i * n + j] = J.matrix[i][j];
    if (nl_eigenvalues(A_flat, n, real_p, imag_p, 2000, 1e-8) != 0) {
        free(A_flat); free(real_p); free(imag_p); return NL_UNSTABLE;
    }
    int all_neg = 1;
    for (size_t i = 0; i < n; i++)
        if (real_p[i] >= -1e-10) { all_neg = 0; break; }
    if (!all_neg) {
        free(A_flat); free(real_p); free(imag_p); return NL_NONHYPERBOLIC;
    }
    double *Q = (double *)calloc(n * n, sizeof(double));
    if (!Q) { free(A_flat); free(real_p); free(imag_p); return NL_UNSTABLE; }
    for (size_t i = 0; i < n; i++) Q[i * n + i] = 1.0;
    int lyap_ok = nl_solve_lyapunov_equation(A_flat, Q, P, n);
    free(Q); free(real_p); free(imag_p);
    if (lyap_ok != 0) { free(A_flat); return NL_UNSTABLE; }
    double max_eig_P = 0.0;
    for (size_t i = 0; i < n; i++) {
        double row_sum = fabs(P[i * n + i]);
        for (size_t j = 0; j < n; j++)
            if (i != j) row_sum += fabs(P[i * n + j]);
        if (row_sum > max_eig_P) max_eig_P = row_sum;
    }
    *radius = 0.1 / sqrt(max_eig_P);
    if (*radius < 1e-6) *radius = 1e-6;
    free(A_flat);
    return NL_ASYMPTOTICALLY_STABLE;
}

/* Transfer Function G(jw) = C*(jwI-A)^{-1}*B + D */

int nl_transfer_function(const nl_linear_system_t *lin, double omega,
                         double *re, double *im) {
    if (!lin || !re || !im || lin->n < 1 || lin->n > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }
    size_t n = lin->n;
    if (n > 8) { errno = ERANGE; return -1; }
    size_t n2 = 2 * n;
    double *M = (double *)calloc(n2 * n2, sizeof(double));
    double *rhs = (double *)calloc(n2, sizeof(double));
    if (!M || !rhs) { free(M); free(rhs); errno = ENOMEM; return -1; }
    for (size_t i = 0; i < n; i++) {
        rhs[i] = lin->B[i];
        for (size_t j = 0; j < n; j++) {
            M[i * n2 + j] = -lin->A[i][j];
            M[i * n2 + n + j] = -omega * ((i == j) ? 1.0 : 0.0);
            M[(i + n) * n2 + j] = omega * ((i == j) ? 1.0 : 0.0);
            M[(i + n) * n2 + n + j] = -lin->A[i][j];
        }
    }
    int ok = gauss_solve(M, rhs, n2);
    if (ok != 0) { free(M); free(rhs); return -1; }
    double g_re = lin->D, g_im = 0.0;
    for (size_t i = 0; i < n; i++) {
        g_re += lin->C[i] * rhs[i];
        g_im += lin->C[i] * rhs[i + n];
    }
    *re = g_re; *im = g_im;
    free(M); free(rhs);
    return 0;
}

/* Circle Criterion */

int nl_circle_criterion(const nl_lure_system_t *lure,
                        const double *freqs, size_t nfreqs,
                        double *margin) {
    if (!lure || !freqs || nfreqs < 1 || !margin) {
        errno = EINVAL; return -1;
    }
    double k1 = lure->nonlinearity.slope_min;
    double k2 = lure->nonlinearity.slope_max;
    if (!isfinite(k1) || !isfinite(k2) || k2 <= k1) {
        errno = EINVAL; return -1;
    }
    double min_re = INFINITY;
    for (size_t k = 0; k < nfreqs; k++) {
        double Gr, Gi;
        if (nl_transfer_function(&lure->linear, freqs[k], &Gr, &Gi) != 0)
            continue;
        double num_re = 1.0 + k2 * Gr;
        double num_im = k2 * Gi;
        double den_re = 1.0 + k1 * Gr;
        double den_im = k1 * Gi;
        double den_mag_sq = den_re * den_re + den_im * den_im;
        if (den_mag_sq < 1e-30) { *margin = -INFINITY; return 0; }
        double Z_re = (num_re * den_re + num_im * den_im) / den_mag_sq;
        if (Z_re < min_re) min_re = Z_re;
    }
    *margin = min_re;
    return (min_re > 1e-8) ? 1 : 0;
}

/* Popov Criterion */

int nl_popov_criterion(const nl_lure_system_t *lure,
                       const double *freqs, size_t nfreqs,
                       double *eta, double *margin) {
    if (!lure || !freqs || nfreqs < 1 || !eta || !margin) {
        errno = EINVAL; return -1;
    }
    double k = lure->nonlinearity.slope_max;
    if (k <= 0.0 || !isfinite(k)) { errno = EINVAL; return -1; }
    double best_min = -INFINITY, best_eta = 0.0;
    int n_eta = 50;
    for (int ie = 0; ie < n_eta; ie++) {
        double test_eta = 10.0 * (double)ie / (double)(n_eta - 1);
        double cur_min = INFINITY;
        for (size_t kf = 0; kf < nfreqs; kf++) {
            double Gr, Gi;
            if (nl_transfer_function(&lure->linear, freqs[kf],
                                     &Gr, &Gi) != 0) continue;
            double val = Gr - test_eta * freqs[kf] * Gi + 1.0 / k;
            if (val < cur_min) cur_min = val;
        }
        if (cur_min > best_min) { best_min = cur_min; best_eta = test_eta; }
    }
    *eta = best_eta;
    *margin = best_min;
    return (best_min > 1e-8) ? 1 : 0;
}

/* Positive Definiteness via Gershgorin */

int nl_is_positive_definite(const double *P, size_t n, double *min_eig) {
    if (!P || n < 1 || n > NL_MAX_STATE_DIM) {
        if (min_eig) *min_eig = -1.0;
        return 0;
    }
    double min_gersh = INFINITY;
    for (size_t i = 0; i < n; i++) {
        double off_diag_sum = 0.0;
        for (size_t j = 0; j < n; j++)
            if (i != j) off_diag_sum += fabs(P[i * n + j]);
        double bound = P[i * n + i] - off_diag_sum;
        if (bound < min_gersh) min_gersh = bound;
    }
    if (min_eig) *min_eig = min_gersh;
    return (min_gersh > 0.0) ? 1 : 0;
}
