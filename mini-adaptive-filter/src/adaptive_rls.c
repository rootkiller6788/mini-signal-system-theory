/**
 * @file    adaptive_rls.c
 * @brief   RLS (Recursive Least Squares) family implementations
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): RLS, QR-RLS (Givens), Kalman gain computation,
 *                     inverse correlation update (Riccati equation)
 *   L2 (Core Concepts): Exponentially weighted least squares,
 *                       deterministic normal equations, matrix inversion lemma
 *   L4 (Fundamental Laws): RLS convergence properties,
 *                          Woodbury matrix identity, Riccati equation
 *   L6 (Canonical Problems): Fast-adapting system identification,
 *                            non-stationary tracking with forgetting factor
 *
 * Reference:
 *   Haykin (2014), "Adaptive Filter Theory", 5th ed., Pearson: Chapters 9-10
 *   Sayed (2008), "Adaptive Filters", Wiley-IEEE Press: Chapters 14-17
 *   Diniz (2013), "Adaptive Filtering", 4th ed., Springer: Chapters 5-6
 *
 * Course Mapping:
 *   Stanford EE264 - Adaptive Filtering (RLS algorithms)
 *   Berkeley EE225D - Adaptive Signal Processing
 *   MIT 6.450 - Digital Communications (RLS equalization)
 *   ETH 227-0436 - Communications (adaptive receiver design)
 */

#include "adaptive_filter.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =========================================================================
 * L5: Standard RLS Algorithm
 * ========================================================================= */

/**
 * @brief af_rls_standard_update: Standard RLS weight update
 *
 * The RLS algorithm minimizes the exponentially weighted sum of
 * squared errors:
 *   J(n) = sum_{i=0}^{n} lambda^{n-i} * |d(i) - w^T * x(i)|^2
 *
 * where lambda is the forgetting factor (0 < lambda <= 1).
 * lambda = 1 gives equal weight to all past data (stationary case).
 * lambda < 1 gives more weight to recent data (tracking capability).
 *
 * Recursive updates (derived via the matrix inversion lemma):
 *
 * Step 1 - Compute Kalman gain:
 *   pi(n) = P(n-1) * x(n)                        [intermediate vector]
 *   gamma(n) = lambda + x^T(n) * pi(n)            [scalar]
 *   K(n) = pi(n) / gamma(n)                       [Kalman gain]
 *
 * Step 2 - Compute a priori error:
 *   xi(n) = d(n) - w^T(n-1) * x(n)
 *
 * Step 3 - Update weights:
 *   w(n) = w(n-1) + K(n) * xi(n)
 *
 * Step 4 - Update inverse correlation matrix (Riccati equation):
 *   P(n) = lambda^{-1} * [P(n-1) - K(n) * pi^T(n)]
 *
 * Initialization:
 *   w(0) = 0
 *   P(0) = delta^{-1} * I   (delta = regularization, small for high SNR,
 *                             large for low SNR)
 *
 * Complexity: O(M^2) per iteration (dominated by matrix-vector products)
 *
 * Convergence properties:
 *   - RLS converges in approximately 2*M iterations (order of M)
 *   - Rate is independent of eigenvalue spread (unlike LMS)
 *   - Steady-state misadjustment with forgetting:
 *     M_RLS ¡Ö (1-lambda) * M / 2  (for lambda close to 1)
 *
 * Numerical stability note:
 *   The direct form (this implementation) can suffer from numerical
 *   instability due to loss of symmetry and positive definiteness of P.
 *   For long-running applications, QR-RLS or square-root RLS is preferred.
 *
 * Reference: Haykin (2014) ¡ì9.2, Algorithm 9.1
 */
void af_rls_standard_update(af_filter_t *filter, double x, double d,
                             double lambda, double delta, size_t order) {
    (void)x; (void)delta;
    size_t i, j;
    double *P, *K, *w, *pi, *xbuf;
    double gamma, xi, inv_gamma;

    if (!filter || !filter->P || !filter->K || order == 0) return;

    P = filter->P;
    K = filter->K;
    w = filter->w;
    xbuf = filter->x_buffer;

    /* pi is stored in residual_vector workspace */
    pi = filter->residual_vector;
    if (!pi) return;

    /* Ensure lambda is in valid range */
    if (lambda <= 0.0) lambda = 0.99;
    if (lambda > 1.0) lambda = 1.0;

    /* Step 1: Compute pi = P * x */
    for (i = 0; i < order; i++) {
        pi[i] = 0.0;
        for (j = 0; j < order; j++) {
            pi[i] += P[i * order + j] * xbuf[j];
        }
    }

    /* Compute gamma = lambda + x^T * pi */
    gamma = lambda;
    for (i = 0; i < order; i++) {
        gamma += xbuf[i] * pi[i];
    }

    /* Check for numerical stability */
    if (fabs(gamma) < 1e-15) {
        gamma = (gamma >= 0) ? 1e-15 : -1e-15;
    }

    /* Step 2: Kalman gain K = pi / gamma */
    inv_gamma = 1.0 / gamma;
    for (i = 0; i < order; i++) {
        K[i] = pi[i] * inv_gamma;
    }

    /* Step 3: A priori error xi = d - w^T * x */
    xi = d;
    for (i = 0; i < order; i++) {
        xi -= w[i] * xbuf[i];
    }

    /* Step 4: Weight update w = w + K * xi */
    for (i = 0; i < order; i++) {
        w[i] += K[i] * xi;
    }

    /* Step 5: P update = lambda^{-1} * (P - K * pi^T) */
    {
        double inv_lambda = 1.0 / lambda;
        for (i = 0; i < order; i++) {
            for (j = 0; j < order; j++) {
                P[i * order + j] = inv_lambda * (P[i * order + j] - K[i] * pi[j]);
            }
        }
    }
}

/* =========================================================================
 * L5: QR-RLS using Givens Rotations
 * ========================================================================= */

/**
 * @brief af_rls_qr_givens_update: QR-RLS using Givens rotations
 *
 * QR-RLS maintains the Cholesky factor R_chol of the inverse
 * correlation matrix P (i.e., P^{-1}^{1/2}) and updates it via
 * Givens rotations, providing significantly better numerical
 * stability than the standard RLS.
 *
 * The key idea: instead of updating P directly, we maintain:
 *   Phi^{1/2}(n) * w(n) = theta(n)
 * where Phi(n) = sum lambda^{n-i} * x(i) * x^T(i) + delta * lambda^n * I
 * is the data correlation matrix (inverse of P, up to scaling).
 *
 * A Givens rotation matrix G = [c  s; -s  c] with c^2 + s^2 = 1
 * is applied to zero out elements below the diagonal of the
 * augmented data matrix, maintaining the upper-triangular structure.
 *
 * Givens rotation (zeroing element y using x as reference):
 *   r = hypot(x, y) = sqrt(x^2 + y^2)
 *   c = x / r
 *   s = y / r
 *   Apply: [c  s; -s  c] * [x] = [r]
 *                           [y]   [0]
 *
 * This implementation applies Givens rotations to update the
 * upper-triangular Cholesky factor stored in filter->P.
 *
 * Note: This is a simplified QR-RLS that maintains the Cholesky
 * factor in the same P matrix. In a production implementation,
 * one would maintain a separate R matrix (upper triangular).
 *
 * Complexity: O(M^2) per iteration but with better numerical
 * properties than standard RLS.
 *
 * Reference: Haykin (2014) ¡ì10.4, "QR-Decomposition-Based RLS Filters"
 *            Gentleman & Kung (1981), "Matrix triangularization by systolic arrays"
 */
void af_rls_qr_givens_update(af_filter_t *filter, double x, double d,
                              double lambda, size_t order) {
    (void)x;
    size_t i, j;
    double sqrt_lambda, c, s, r;

    if (!filter || !filter->P || order == 0) return;

    if (lambda <= 0.0) lambda = 0.99;
    if (lambda > 1.0) lambda = 1.0;
    sqrt_lambda = sqrt(lambda);

    double *R_chol = filter->P;  /* Use P to store Cholesky factor R */
    double *u = filter->residual_vector;  /* Workspace */
    double *d_rot = filter->K;   /* Workspace for rotated desired vector */

    if (!u || !d_rot) return;

    /* Step 1: Form the new row of the data matrix:
     * [sqrt(lambda)*R(n-1) | sqrt(lambda)*p(n-1)]
     * [x^T(n)              | d(n)              ]
     *
     * Then apply Givens rotations to restore upper-triangular form.
     */

    /* Copy current x into u (the new row to be rotated in) */
    for (i = 0; i < order; i++) {
        u[i] = filter->x_buffer[i];
    }
    /* d_rot stores the transformed desired signal */
    *d_rot = d;

    /* Scale existing R by sqrt(lambda) for exponential weighting */
    for (i = 0; i < order; i++) {
        for (j = i; j < order; j++) {
            R_chol[i * order + j] *= sqrt_lambda;
        }
        /* Also scale the right-hand side (stored implicitly) */
        if (filter->w) {
            filter->w[i] *= sqrt_lambda;
        }
    }

    /* Step 2: Apply Givens rotations to zero out elements of u
     * against the rows of R_chol */
    for (i = 0; i < order; i++) {
        if (fabs(u[i]) < 1e-20) continue;

        /* Compute Givens rotation to zero u[i] using R_chol[i][i] */
        double r_ii = R_chol[i * order + i];
        r = hypot(r_ii, u[i]);

        if (r < 1e-20) continue;

        c = r_ii / r;
        s = u[i] / r;

        /* Apply rotation to row i of R_chol and the new row u */
        R_chol[i * order + i] = r;

        for (j = i + 1; j < order; j++) {
            double r_ij = R_chol[i * order + j];
            double u_j = u[j];
            R_chol[i * order + j] = c * r_ij + s * u_j;
            u[j] = -s * r_ij + c * u_j;
        }

        /* Apply rotation to the right-hand side (desired signal) */
        {
            double d_i = (filter->w) ? filter->w[i] : 0.0;
            double d_new = *d_rot;
            if (filter->w) filter->w[i] = c * d_i + s * d_new;
            *d_rot = -s * d_i + c * d_new;
        }
    }

    /* After all rotations, the weight vector can be solved by
     * back-substitution: R * w = p (where p is the rotated desired vector).
     * For simplicity in this implementation, the updated w is stored
     * in filter->w during the rotation process above.
     *
     * Full back-substitution would be needed for a more rigorous
     * implementation that explicitly solves R * w_new = p_new.
     */
}

/* =========================================================================
 * L5: Kalman Gain Computation (RLS)
 * ========================================================================= */

/**
 * @brief af_rls_compute_kalman_gain: Compute RLS Kalman gain vector
 *
 * The Kalman gain K(n) is the key quantity in RLS filtering:
 *
 *   K(n) = P(n-1) * x(n) / (lambda + x^T(n) * P(n-1) * x(n))
 *
 * This is equivalent to:
 *   K = (P * x) / (lambda + x^T * P * x)
 *
 * The denominator gamma(n) = lambda + x^T * P * x is the
 * "conversion factor" that relates a priori and a posteriori errors:
 *   e_post(n) = gamma(n) * xi(n)
 *
 * All eigenvalues of K(n) * x^T(n) are in [0, 1], ensuring
 * stability of the Riccati update for P.
 *
 * @param P      Inverse correlation matrix [order*order] (row-major, as in af_matrix_t)
 * @param x      Input vector [order]
 * @param lambda Forgetting factor
 * @param order  Filter order
 * @param K      Output: Kalman gain vector [order]
 * @return       gamma = lambda + x^T * P * x (denominator)
 */
double af_rls_compute_kalman_gain(const af_matrix_t *P, const double *x,
                                   double lambda, size_t order, double *K) {
    size_t i, j;
    double gamma;

    if (!P || !x || !K || order == 0) return 1.0;

    /* Compute pi = P * x */
    double *pi = (double *)malloc(order * sizeof(double));
    if (!pi) return 1.0;

    for (i = 0; i < order; i++) {
        pi[i] = 0.0;
        for (j = 0; j < order; j++) {
            pi[i] += P->data[i * P->stride + j] * x[j];
        }
    }

    /* gamma = lambda + x^T * pi */
    gamma = lambda;
    for (i = 0; i < order; i++) {
        gamma += x[i] * pi[i];
    }

    if (fabs(gamma) < 1e-15) gamma = (gamma >= 0) ? 1e-15 : -1e-15;

    /* K = pi / gamma */
    double inv_gamma = 1.0 / gamma;
    for (i = 0; i < order; i++) {
        K[i] = pi[i] * inv_gamma;
    }

    free(pi);
    return gamma;
}

/* =========================================================================
 * L5: Inverse Correlation Matrix Update (RLS Riccati Equation)
 * ========================================================================= */

/**
 * @brief af_rls_update_inverse_correlation: RLS Riccati update for P
 *
 * Updates the inverse correlation matrix using the Riccati equation:
 *
 *   P(n) = lambda^{-1} * (I - K(n) * x^T(n)) * P(n-1)
 *
 * or equivalently:
 *   P = lambda^{-1} * (P - K * (x^T * P))
 *
 * This update is the matrix form of the discrete Riccati equation.
 * The term K * x^T is a rank-1 matrix, making this a rank-1 update
 * to the inverse correlation matrix.
 *
 * Woodbury matrix identity (used in deriving RLS):
 *   (A + U*C*V)^{-1} = A^{-1} - A^{-1}*U*(C^{-1} + V*A^{-1}*U)^{-1}*V*A^{-1}
 *
 * In RLS: A = lambda*P^{-1}(n-1), U = V^T = x(n), C = 1
 * giving the rank-1 update formula used in RLS.
 *
 * Numerical note:
 *   Direct implementation of this formula can lead to loss of
 *   symmetry due to floating-point errors. In practice, only the
 *   upper triangular part of P is updated and the lower triangle
 *   is set by symmetry.
 *
 * @param P      Inverse correlation matrix [order*order] (in/out)
 * @param K      Kalman gain vector [order]
 * @param x      Input vector [order]
 * @param lambda Forgetting factor
 * @param order  Filter order
 */
void af_rls_update_inverse_correlation(af_matrix_t *P, const double *K,
                                        const double *x, double lambda,
                                        size_t order) {
    size_t i, j;
    double *temp;
    double inv_lambda;

    if (!P || !K || !x || order == 0) return;

    if (lambda <= 0.0) lambda = 0.99;
    inv_lambda = 1.0 / lambda;

    /* Compute temp = x^T * P (row vector) */
    temp = (double *)malloc(order * sizeof(double));
    if (!temp) return;

    for (j = 0; j < order; j++) {
        temp[j] = 0.0;
        for (i = 0; i < order; i++) {
            temp[j] += x[i] * P->data[i * P->stride + j];
        }
    }

    /* P = lambda^{-1} * (P - K * temp)
     * where (K * temp)[i][j] = K[i] * temp[j] (rank-1 outer product) */
    for (i = 0; i < order; i++) {
        for (j = 0; j < order; j++) {
            P->data[i * P->stride + j] =
                inv_lambda * (P->data[i * P->stride + j] - K[i] * temp[j]);
        }
    }

    free(temp);
}

/* =========================================================================
 * L6: Fast-Adapting System Identification with RLS
 * ========================================================================= */

/**
 * @brief af_rls_system_identify: Identify unknown system using RLS
 *
 * Performs system identification using RLS, which converges in
 * approximately 2*M iterations regardless of eigenvalue spread.
 * This is significantly faster than LMS for colored input signals.
 *
 * The forgetting factor lambda controls the trade-off between
 * steady-state error and tracking capability:
 *   lambda close to 1: Low misadjustment, slow tracking
 *   lambda < 1:        Higher misadjustment, fast tracking
 *
 * Typical lambda values:
 *   Stationary system:              lambda = 1.0
 *   Slowly varying:                 lambda = 0.995 - 0.999
 *   Moderately non-stationary:      lambda = 0.98 - 0.995
 *   Fast time-varying:              lambda = 0.95 - 0.98
 *
 * @param x            Input signal [n_samples]
 * @param d            Observed output [n_samples]
 * @param n_samples    Training samples
 * @param order        System order
 * @param lambda       Forgetting factor
 * @param delta        Regularization
 * @param w_estimated  Output: estimated impulse response [order]
 */
void af_rls_system_identify(const double *x, const double *d,
                             size_t n_samples, size_t order,
                             double lambda, double delta,
                             double *w_estimated) {
    size_t n, i;
    af_config_t cfg;
    af_filter_t filter;

    if (!x || !d || !w_estimated || n_samples == 0 || order == 0) return;

    /* Initialize RLS filter */
    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_RLS;
    cfg.order = order;
    cfg.mu = 0.01;          /* Not used by RLS but must be > 0 */
    cfg.lambda = lambda;
    cfg.delta = delta;
    cfg.epsilon = 1e-6;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return;

    /* Run RLS over all samples */
    for (n = 0; n < n_samples; n++) {
        af_rls_standard_update(&filter, x[n], d[n], lambda, delta, order);
    }

    /* Extract estimated weights */
    for (i = 0; i < order; i++) w_estimated[i] = filter.w[i];

    af_free(&filter);
}

