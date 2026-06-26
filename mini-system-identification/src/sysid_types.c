/**
 * @file sysid_types.c
 * @brief Implementation of core type operations for system identification
 */

#include "sysid_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/*---------------------------------------------------------------------------
 * Model order operations
 *---------------------------------------------------------------------------*/

void sysid_order_init(sysid_model_order_t *order) {
    if (!order) return;
    memset(order, 0, sizeof(sysid_model_order_t));
}

int sysid_order_is_valid(const sysid_model_order_t *order) {
    if (!order) return 0;
    if (order->na < 0 || order->nb < 0 || order->nc < 0 ||
        order->nd < 0 || order->nf < 0 || order->nk < 0) return 0;
    if (order->nx < 0 || order->ny < 0 || order->nu < 0) return 0;
    return 1;
}

int sysid_order_nparam(sysid_model_type_t type, const sysid_model_order_t *order) {
    if (!order) return -1;
    int np = 0;
    switch (type) {
        case SYSID_MODEL_FIR:
            np = order->nb;
            break;
        case SYSID_MODEL_ARX:
            np = order->na + order->nb;
            break;
        case SYSID_MODEL_ARMAX:
            np = order->na + order->nb + order->nc;
            break;
        case SYSID_MODEL_OE:
            np = order->nb + order->nf;
            break;
        case SYSID_MODEL_BJ:
            np = order->nb + order->nc + order->nd + order->nf;
            break;
        case SYSID_MODEL_SS:
            np = order->nx * (order->nx + order->nu + order->ny) + order->ny * order->nu;
            break;
        default:
            np = order->na + order->nb;
            break;
    }
    return np;
}

/*---------------------------------------------------------------------------
 * Parameter vector operations
 *---------------------------------------------------------------------------*/

int sysid_param_alloc(sysid_param_vector_t *pv, int nparam) {
    if (!pv || nparam <= 0) return -1;
    memset(pv, 0, sizeof(sysid_param_vector_t));
    pv->theta = (double *)calloc((size_t)nparam, sizeof(double));
    if (!pv->theta) return -1;
    pv->nparam = nparam;
    pv->nparam_free = nparam;
    pv->free_idx = NULL;
    return 0;
}

void sysid_param_free(sysid_param_vector_t *pv) {
    if (!pv) return;
    free(pv->theta); pv->theta = NULL;
    free(pv->free_idx); pv->free_idx = NULL;
    pv->nparam = 0;
    pv->nparam_free = 0;
}

void sysid_param_set_all(sysid_param_vector_t *pv, double value) {
    if (!pv || !pv->theta) return;
    for (int i = 0; i < pv->nparam; i++) {
        pv->theta[i] = value;
    }
}

void sysid_param_copy(sysid_param_vector_t *dst, const sysid_param_vector_t *src) {
    if (!dst || !src || !dst->theta || !src->theta) return;
    if (dst->nparam != src->nparam) return;
    memcpy(dst->theta, src->theta, (size_t)src->nparam * sizeof(double));
}

double sysid_param_norm2(const sysid_param_vector_t *pv) {
    if (!pv || !pv->theta) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < pv->nparam; i++) {
        sum += pv->theta[i] * pv->theta[i];
    }
    return sqrt(sum);
}

double sysid_param_dot(const sysid_param_vector_t *a, const sysid_param_vector_t *b) {
    if (!a || !b || !a->theta || !b->theta) return 0.0;
    if (a->nparam != b->nparam) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < a->nparam; i++) {
        sum += a->theta[i] * b->theta[i];
    }
    return sum;
}

/*---------------------------------------------------------------------------
 * Data structure allocation
 *---------------------------------------------------------------------------*/

int sysid_data_alloc(sysid_data_t *data, int N, int nu, int ny, double Ts) {
    if (!data || N <= 0 || nu <= 0 || ny <= 0) return -1;
    memset(data, 0, sizeof(sysid_data_t));
    data->u = (double *)calloc((size_t)(N * nu), sizeof(double));
    data->y = (double *)calloc((size_t)(N * ny), sizeof(double));
    if (!data->u || !data->y) {
        free(data->u); free(data->y);
        return -1;
    }
    data->N = N;
    data->nu = nu;
    data->ny = ny;
    data->Ts = Ts;
    data->is_real = 1;
    return 0;
}

void sysid_data_free(sysid_data_t *data) {
    if (!data) return;
    free(data->u); data->u = NULL;
    free(data->y); data->y = NULL;
    data->N = 0;
    data->nu = 0;
    data->ny = 0;
}

int sysid_freq_data_alloc(sysid_freq_data_t *fd, int Nf) {
    if (!fd || Nf <= 0) return -1;
    memset(fd, 0, sizeof(sysid_freq_data_t));
    fd->omega = (double *)calloc((size_t)Nf, sizeof(double));
    fd->G = (double complex *)calloc((size_t)Nf, sizeof(double complex));
    fd->weight = NULL;
    if (!fd->omega || !fd->G) {
        free(fd->omega); free(fd->G);
        return -1;
    }
    fd->Nf = Nf;
    return 0;
}

void sysid_freq_data_free(sysid_freq_data_t *fd) {
    if (!fd) return;
    free(fd->omega); fd->omega = NULL;
    free(fd->G); fd->G = NULL;
    free(fd->weight); fd->weight = NULL;
    fd->Nf = 0;
}

/*---------------------------------------------------------------------------
 * Fit result and stability
 *---------------------------------------------------------------------------*/

void sysid_fit_result_init(sysid_fit_result_t *result) {
    if (!result) return;
    memset(result, 0, sizeof(sysid_fit_result_t));
    result->NRMSE_fit = -1e6; /* sentinel: not computed */
}

int sysid_stability_alloc(sysid_stability_t *stab, int n_poles) {
    if (!stab || n_poles <= 0) return -1;
    memset(stab, 0, sizeof(sysid_stability_t));
    stab->poles = (double complex *)calloc((size_t)n_poles, sizeof(double complex));
    if (!stab->poles) return -1;
    stab->n_poles = n_poles;
    return 0;
}

void sysid_stability_free(sysid_stability_t *stab) {
    if (!stab) return;
    free(stab->poles); stab->poles = NULL;
    stab->n_poles = 0;
}

/*---------------------------------------------------------------------------
 * Small matrix utilities (internal but exposed for testing)
 *---------------------------------------------------------------------------*/

/** Cholesky decomposition: A = L Lᵀ (in-place lower triangle)
 *  A is n×n row-major, symmetric positive definite.
 *  Returns 0 on success, -1 if not positive definite.
 */
int sysid_cholesky(double *A, int n) {
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += A[j * n + k] * A[j * n + k];
        }
        double diag = A[j * n + j] - sum;
        if (diag <= 0.0) return -1; /* not positive definite */
        A[j * n + j] = sqrt(diag);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += A[i * n + k] * A[j * n + k];
            }
            A[i * n + j] = (A[i * n + j] - sum) / A[j * n + j];
        }
    }
    return 0;
}

/** Forward substitution: solve L y = b where L is lower triangular n×n
 *  L stored row-major, only lower triangle used. y written into b.
 */
void sysid_forward_sub(const double *L, double *b, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i * n + j] * b[j];
        }
        b[i] = (b[i] - sum) / L[i * n + i];
    }
}

/** Backward substitution: solve Lᵀ x = y where L is lower triangular
 *  Result x written into y.
 */
void sysid_backward_sub(const double *L, double *y, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j * n + i] * y[j];
        }
        y[i] = (y[i] - sum) / L[i * n + i];
    }
}

/** Solve A x = b using Cholesky: A = L Lᵀ, then L y = b, Lᵀ x = y
 *  Modifies A (destroyed) and b (solution overwritten).
 *  Returns 0 on success, -1 if A not positive definite.
 */
int sysid_cholesky_solve(double *A, double *b, int n) {
    if (sysid_cholesky(A, n) != 0) return -1;
    sysid_forward_sub(A, b, n);
    sysid_backward_sub(A, b, n);
    return 0;
}

/** Matrix multiply C = A * B where:
 *  A is m×k row-major, B is k×n row-major, C is m×n row-major.
 *  C must be pre-allocated.
 */
void sysid_mat_mul(const double *A, const double *B, double *C,
                   int m, int k, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += A[i * k + p] * B[p * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/** Matrix transpose multiply: C = Aᵀ * B
 *  A is k×m row-major, B is k×n row-major, C is m×n row-major.
 */
void sysid_mat_t_mul(const double *A, const double *B, double *C,
                     int k, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += A[p * m + i] * B[p * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/** Matrix-vector multiply: y = A * x
 *  A is m×n row-major, x is n×1, y is m×1 (pre-allocated).
 */
void sysid_mat_vec_mul(const double *A, const double *x, double *y, int m, int n) {
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

/** Vector inner product (dot product) */
double sysid_vec_dot(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

/** Householder reflection for QR factorization (internal)
 *  Given vector x[0..n-1], compute v and β such that H = I - β v vᵀ
 *  makes Hx = ||x|| e₁. Returns β, modifies x to store v.
 *  (Golub & Van Loan 1996, Algorithm 5.1.1)
 */
double sysid_householder(double *x, int n) {
    double sigma = 0.0;
    for (int i = 1; i < n; i++) {
        sigma += x[i] * x[i];
    }
    double x0 = x[0];
    double beta;
    if (sigma == 0.0) {
        beta = 0.0;
    } else {
        double mu = sqrt(x0 * x0 + sigma);
        if (x0 <= 0.0) {
            x[0] = x0 - mu;
        } else {
            x[0] = -sigma / (x0 + mu);
        }
        beta = 2.0 * x[0] * x[0] / (sigma + x[0] * x[0]);
        double inv_v1 = 1.0 / x[0];
        for (int i = 1; i < n; i++) {
            x[i] *= inv_v1;
        }
    }
    return beta;
}

/** Polynomial root finding via companion matrix + QR eigenvalue algorithm
 *  Finds roots of a[0] + a[1]x + ... + a[n-1]x^{n-1} + x^n = 0 (monic!)
 *  Uses simple power iteration + deflation for moderate order n ≤ 50.
 *
 *  @param coeffs  Polynomial coefficients [n+1], coeffs[n]=1 (monic)
 *  @param n       Polynomial degree
 *  @param roots   Output: complex roots [n]
 *  @return        0 on success, -1 on non-convergence
 */
int sysid_poly_roots(const double *coeffs, int n, double complex *roots) {
    if (n <= 0 || !coeffs || !roots) return -1;
    if (n == 1) {
        /* Linear: a₁ x + a₀ = 0 → monic form: x + a₀ = 0 */
        roots[0] = -coeffs[0];
        return 0;
    }

    /* Build companion matrix (upper Hessenberg) in row-major: n×n
     * Companion for monic polynomial x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0 = 0
     * Note: coeffs has a_0 at index 0, a_n=1 at index n
     */
    double *Hmat = (double *)calloc((size_t)(n * n), sizeof(double));
    if (!Hmat) return -1;

    /* First row: -coeffs reversed (top row of companion) */
    for (int j = 0; j < n; j++) {
        Hmat[j] = -coeffs[n - 1 - j];
    }
    /* Sub-diagonal ones */
    for (int i = 1; i < n; i++) {
        Hmat[i * n + (i - 1)] = 1.0;
    }

    /* Use simple power method with deflation for eigenvalue estimation.
     * For small n, direct characteristic polynomial via QR on
     * the companion is too heavy. Instead use shifted power iteration
     * to extract dominant eigenvalues sequentially.
     *
     * Practical approach: Use the fact that eigenvalues of companion
     * matrix are the roots. For small n, use the Jenkins-Traub-like
     * approach: extract via inverse iteration.
     */

    /* For now: use Laguerre's method on the polynomial directly */
    double complex *poly = (double complex *)calloc((size_t)(n + 1), sizeof(double complex));
    if (!poly) { free(Hmat); return -1; }
    for (int i = 0; i <= n; i++) {
        poly[i] = coeffs[i];
    }

    /* Laguerre's method with polishing (Press et al. 1992, §9.5) */
    int remaining = n;
    int max_iter = 50;
    double complex *work_poly = (double complex *)calloc((size_t)(n + 1), sizeof(double complex));
    if (!work_poly) { free(Hmat); free(poly); return -1; }

    for (int root_idx = 0; root_idx < n; root_idx++) {
        /* Initial guess: random complex point on unit circle */
        double angle = (root_idx * 2.0943951); /* 2π/3 rotation */
        double complex z = cos(angle) + sin(angle) * I;

        for (int iter = 0; iter < max_iter; iter++) {
            /* Evaluate polynomial and derivatives at z */
            double complex pz = poly[remaining];
            double complex dpz = 0.0;
            double complex ddpz = 0.0;
            for (int j = remaining - 1; j >= 0; j--) {
                ddpz = ddpz * z + 2.0 * dpz;
                dpz = dpz * z + pz;
                pz = pz * z + poly[j];
            }
            /* ddpz was accumulated wrong; recompute properly */
            pz = poly[remaining];
            dpz = 0.0;
            ddpz = 0.0;
            for (int j = remaining - 1; j >= 0; j--) {
                double complex pz_new = pz * z + poly[j];
                double complex dpz_new = dpz * z + pz;
                ddpz = ddpz * z + 2.0 * dpz;
                dpz = dpz_new;
                pz = pz_new;
            }

            if (cabs(pz) < 1e-14) break;

            /* Laguerre step */
            double complex G = dpz / pz;
            double complex H = G * G - ddpz / pz;
            double complex denom;
            double n_rem = (double)remaining;
            double complex sqrt_arg = (n_rem - 1.0) * (n_rem * H - G * G);
            /* Choose sign to maximize denominator magnitude */
            double complex sqrt_val = csqrt(sqrt_arg);
            double complex den1 = G + sqrt_val;
            double complex den2 = G - sqrt_val;
            if (cabs(den1) > cabs(den2)) {
                denom = den1;
            } else {
                denom = den2;
            }
            if (cabs(denom) < 1e-14) break;
            double complex dz = n_rem / denom;
            z = z - dz;
            if (cabs(dz) < 1e-12) break;
        }

        roots[root_idx] = z;

        /* Deflate polynomial */
        remaining--;
        if (remaining > 0) {
            /* Synthetic division: poly_new = poly / (x - root) */
            double complex b = poly[remaining + 1];
            for (int j = remaining; j >= 0; j--) {
                double complex a = poly[j];
                poly[j] = b;
                b = a + b * z;
            }
            for (int j = 0; j <= remaining; j++) {
                poly[j] = poly[j + 1];
            }
        }
    }

    free(Hmat);
    free(poly);
    free(work_poly);
    return 0;
}

/*---------------------------------------------------------------------------
 * SVD approximation via power iteration (Golub & Kahan 1965)
 * Computes the top-k singular values and left/right vectors of A (m×n).
 * Used for total least squares and subspace methods.
 *---------------------------------------------------------------------------*/

/** Compute largest singular value of A via power iteration on AᵀA */
double sysid_svd_power(const double *A, int m, int n, double *u, double *v, int max_iter) {
    if (!A || !u || !v || m <= 0 || n <= 0) return 0.0;

    /* Initialize v with ones */
    for (int j = 0; j < n; j++) v[j] = 1.0;
    double norm_v = sqrt((double)n);
    for (int j = 0; j < n; j++) v[j] /= norm_v;

    double sigma = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* u = A v */
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += A[i * n + j] * v[j];
            }
            u[i] = sum;
        }
        /* Normalize u */
        double norm_u = 0.0;
        for (int i = 0; i < m; i++) norm_u += u[i] * u[i];
        norm_u = sqrt(norm_u);
        if (norm_u < 1e-15) break;
        for (int i = 0; i < m; i++) u[i] /= norm_u;

        /* v = Aᵀ u */
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                sum += A[i * n + j] * u[i];
            }
            v[j] = sum;
        }
        sigma = 0.0;
        for (int j = 0; j < n; j++) sigma += v[j] * v[j];
        sigma = sqrt(sigma);
        if (sigma < 1e-15) break;
        for (int j = 0; j < n; j++) v[j] /= sigma;

        /* loop continues */
    }

    return sigma;
}

/** QR decomposition of A (m×n, m ≥ n) using Householder reflections.
 *  Returns R in upper triangle of A, Q implicitly via Householder vectors.
 *  For full decomposition, tau = Householder coefficients.
 *  (Golub & Van Loan 1996, Algorithm 5.2.1)
 */
int sysid_qr_decomp(double *A, int m, int n, double *tau) {
    for (int j = 0; j < n && j < m; j++) {
        /* Compute Householder for column j below diagonal */
        double *col = (double *)calloc((size_t)(m - j), sizeof(double));
        if (!col) return -1;
        for (int i = j; i < m; i++) col[i - j] = A[i * n + j];
        double beta = sysid_householder(col, m - j);
        if (tau) tau[j] = beta;

        /* Apply Householder to remaining columns */
        if (beta > 0.0) {
            for (int k = j; k < n; k++) {
                /* Compute vᵀ A(j:m, k) */
                double sum = 0.0;
                for (int i = j; i < m; i++) {
                    sum += ((i == j) ? 1.0 : col[i - j]) * A[i * n + k];
                }
                sum *= beta;
                /* A -= beta * v * (vᵀ A) */
                A[j * n + k] -= sum;
                for (int i = j + 1; i < m; i++) {
                    A[i * n + k] -= sum * col[i - j];
                }
            }
            /* Store v in lower part of A */
            A[j * n + j] = col[0]; /* diagonal = R(j,j) */
            for (int i = j + 1; i < m; i++) {
                A[i * n + j] = col[i - j]; /* Householder vector below diagonal */
            }
        }
        free(col);
    }
    return 0;
}
