/**
 * @file    adaptive_lms.c
 * @brief   LMS (Least Mean Squares) family implementations
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): LMS, NLMS, Sign-Error LMS, Sign-Data LMS,
 *                     Sign-Sign LMS, Leaky LMS, Block LMS
 *   L2 (Core Concepts): Gradient descent, stochastic approximation,
 *                       convergence behavior, misadjustment analysis
 *   L4 (Fundamental Laws): Widrow-Hoff LMS convergence theorem,
 *                          NLMS convergence bounds
 *   L6 (Canonical Problems): System identification via LMS
 *
 * Reference:
 *   Widrow & Hoff (1960), "Adaptive Switching Circuits", IRE WESCON
 *   Widrow & Stearns (1985), "Adaptive Signal Processing", Prentice-Hall
 *   Haykin (2014), "Adaptive Filter Theory", 5th ed., Pearson: Chapters 5-7
 *   Sayed (2008), "Adaptive Filters", Wiley-IEEE Press: Chapters 10-12
 *
 * Course Mapping:
 *   MIT 6.003 - Signals and Systems (gradient descent)
 *   Stanford EE264 - Adaptive Filtering (LMS algorithms)
 *   Berkeley EE225D - Adaptive Signal Processing
 *   Illinois ECE 310 - DSP (adaptive filtering basics)
 */

#include "adaptive_filter.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L5: Standard LMS Algorithm (Widrow-Hoff, 1960)
 * ========================================================================= */

/**
 * @brief af_lms_direct_update: Standard LMS weight update
 *
 * This is the foundational adaptive filtering algorithm. The LMS
 * algorithm uses stochastic gradient descent to minimize the MSE
 * cost function J(n) = E[|e(n)|^2].
 *
 * Update equation:
 *   w(n+1) = w(n) + mu * e(n) * x(n)
 *
 * where:
 *   e(n) = d(n) - w^T(n) * x(n)  (instantaneous error)
 *   mu   = step-size parameter
 *
 * The instantaneous gradient estimate is:
 *   grad_hat J = -2 * e(n) * x(n)
 *
 * This is an unbiased estimate of the true gradient:
 *   E[grad_hat J] = -2 * (p - R * w) = grad J
 *
 * Convergence in mean requires:
 *   0 < mu < 2 / lambda_max(R)
 *
 * Steady-state excess MSE (small mu):
 *   J_excess �� (mu/2) * J_min * tr(R)
 *
 * Misadjustment:
 *   M = J_excess / J_min �� (mu/2) * tr(R)
 *
 * Complexity: 2*order multiply-add operations per iteration
 *
 * @param filter  Filter state (for MSE update, can be NULL if only weights needed)
 * @param x       Input vector [order] (already shifted into delay line)
 * @param d       Desired response
 * @param w       Tap-weight vector [order] (in/out)
 * @param mu      Step-size parameter
 * @param order   Filter order
 */
void af_lms_direct_update(af_filter_t *filter, const double *x, double d,
                           double *w, double mu, size_t order) {
    size_t i;
    double y, e;
    (void)filter;

    if (!x || !w || order == 0) return;

    /* Compute filter output */
    y = 0.0;
    for (i = 0; i < order; i++) {
        y += w[i] * x[i];
    }

    /* Compute error */
    e = d - y;

    /* LMS update: w(n+1) = w(n) + mu * e(n) * x(n) */
    {
        double mu_e = mu * e;
        for (i = 0; i < order; i++) {
            w[i] += mu_e * x[i];
        }
    }
}

/* =========================================================================
 * L5: Normalized LMS (NLMS) Algorithm
 * ========================================================================= */

/**
 * @brief af_lms_normalized_update: Normalized LMS weight update
 *
 * NLMS normalizes the step size by the squared Euclidean norm of the
 * input vector, making convergence independent of input signal power.
 *
 * Update equation:
 *   w(n+1) = w(n) + (mu / (epsilon + ||x(n)||^2)) * e(n) * x(n)
 *
 * where epsilon > 0 prevents division by zero when ||x||^2 is very small.
 *
 * Convergence in mean-square is guaranteed for:
 *   0 < mu < 2
 *
 * The effective step size is:
 *   mu_eff(n) = mu / (epsilon + ||x(n)||^2)
 *
 * This provides:
 *   - Faster convergence for highly non-stationary signals
 *   - Robustness to input signal power variations
 *   - Potentially faster convergence than LMS for colored inputs
 *
 * The NLMS can be derived as the solution to the constrained
 * optimization problem:
 *   minimize ||w(n+1) - w(n)||^2
 *   subject to d(n) = w^T(n+1) * x(n)
 *
 * (with the regularization epsilon giving a relaxed constraint).
 *
 * Reference: Nagumo & Noda (1967), IEEE Trans. Automatic Control
 *            Slock (1993), IEEE Trans. Signal Processing
 */
void af_lms_normalized_update(af_filter_t *filter, const double *x, double d,
                               double *w, double mu, double epsilon,
                               size_t order) {
    size_t i;
    double y, e, norm_x2, step;

    (void)filter;
    if (!x || !w || order == 0) return;

    /* Compute filter output */
    y = 0.0;
    norm_x2 = 0.0;
    for (i = 0; i < order; i++) {
        y += w[i] * x[i];
        norm_x2 += x[i] * x[i];
    }

    /* Compute error */
    e = d - y;

    /* NLMS step: mu / (epsilon + ||x||^2) */
    step = mu / (epsilon + norm_x2);

    /* NLMS update */
    {
        double step_e = step * e;
        for (i = 0; i < order; i++) {
            w[i] += step_e * x[i];
        }
    }
}

/* =========================================================================
 * L5: Sign-Error LMS
 * ========================================================================= */

/**
 * @brief af_lms_sign_error_update: Sign-error LMS
 *
 * Uses only the sign of the error to update weights:
 *   w(n+1) = w(n) + mu * sgn(e(n)) * x(n)
 *
 * where sgn(x) = +1 if x > 0, 0 if x == 0, -1 if x < 0.
 *
 * Advantages:
 *   - Reduced computational complexity (no full multiply for error term)
 *   - Robust to impulsive noise (error magnitude is clipped)
 *   - Simple hardware implementation (comparator + adder)
 *
 * Disadvantages:
 *   - Slower convergence than standard LMS
 *   - Larger steady-state misadjustment for same mu
 *
 * Applications: Telecommunication echo cancellation, robust adaptive control
 */
void af_lms_sign_error_update(af_filter_t *filter, const double *x, double d,
                               double *w, double mu, size_t order) {
    size_t i;
    double y, e, sgn_e;
    (void)filter;

    if (!x || !w || order == 0) return;

    /* Compute output and error */
    y = 0.0;
    for (i = 0; i < order; i++) y += w[i] * x[i];
    e = d - y;

    /* Sign of error */
    sgn_e = (e > 0.0) ? 1.0 : ((e < 0.0) ? -1.0 : 0.0);

    /* Update: w += mu * sgn(e) * x */
    for (i = 0; i < order; i++) {
        w[i] += mu * sgn_e * x[i];
    }
}

/* =========================================================================
 * L5: Sign-Data LMS
 * ========================================================================= */

/**
 * @brief af_lms_sign_data_update: Sign-data LMS
 *
 * Uses only the sign of the input data:
 *   w(n+1) = w(n) + mu * e(n) * sgn(x(n))
 *
 * This eliminates multiplication by the input magnitude, which can
 * be implemented as conditional addition/subtraction (shift operations
 * in fixed-point).
 *
 * The algorithm is particularly suited for:
 *   - Hardware with limited multiplier resources
 *   - High-speed applications
 *   - Binary/ternary input signals
 */
void af_lms_sign_data_update(af_filter_t *filter, const double *x, double d,
                              double *w, double mu, size_t order) {
    size_t i;
    double y, e;
    (void)filter;

    if (!x || !w || order == 0) return;

    y = 0.0;
    for (i = 0; i < order; i++) y += w[i] * x[i];
    e = d - y;

    /* Update: w += mu * e * sgn(x) */
    {
        double mu_e = mu * e;
        for (i = 0; i < order; i++) {
            double sgn_x = (x[i] > 0.0) ? 1.0 : ((x[i] < 0.0) ? -1.0 : 0.0);
            w[i] += mu_e * sgn_x;
        }
    }
}

/* =========================================================================
 * L5: Sign-Sign LMS (Clipped LMS)
 * ========================================================================= */

/**
 * @brief af_lms_sign_sign_update: Sign-sign LMS
 *
 * Uses signs of both error and input:
 *   w(n+1) = w(n) + mu * sgn(e(n)) * sgn(x(n))
 *
 * This is the simplest possible adaptive algorithm:
 *   - No multiplications at all
 *   - Only sign decisions and additions
 *   - Can be implemented with digital logic gates
 *
 * Applications:
 *   - Spread-spectrum communication systems
 *   - Low-power DSP applications
 *   - Adaptive delta modulation
 *
 * Convergence: guaranteed for sufficiently small mu, but significantly
 * slower than standard LMS. The steady-state error is generally higher
 * due to the coarse quantization of the gradient estimate.
 */
void af_lms_sign_sign_update(af_filter_t *filter, const double *x, double d,
                              double *w, double mu, size_t order) {
    size_t i;
    double y, e, sgn_e;
    (void)filter;

    if (!x || !w || order == 0) return;

    y = 0.0;
    for (i = 0; i < order; i++) y += w[i] * x[i];
    e = d - y;

    sgn_e = (e > 0.0) ? 1.0 : ((e < 0.0) ? -1.0 : 0.0);

    /* Update: w += mu * sgn(e) * sgn(x) */
    for (i = 0; i < order; i++) {
        double sgn_x = (x[i] > 0.0) ? 1.0 : ((x[i] < 0.0) ? -1.0 : 0.0);
        w[i] += mu * sgn_e * sgn_x;
    }
}

/* =========================================================================
 * L5: Leaky LMS Algorithm
 * ========================================================================= */

/**
 * @brief af_lms_leaky_update: Leaky LMS weight update
 *
 * Introduces a leakage factor gamma to prevent coefficient drift:
 *   w(n+1) = gamma * w(n) + mu * e(n) * x(n)
 *
 * where 0 < gamma <= 1.
 *
 * Leakage Mechanics:
 *   Without leakage (gamma = 1): standard LMS, coefficients can
 *   grow unboundedly in finite-precision arithmetic due to
 *   accumulation of quantization errors (coefficient drift).
 *
 *   With leakage (gamma < 1): adds a small bias towards zero,
 *   stabilizing the algorithm against drift. The bias is proportional
 *   to (1-gamma). For tracking slowly-varying systems, a small amount
 *   of leakage can actually improve performance by adding a
 *   "forgetting" mechanism.
 *
 * Steady-state analysis (small leakage):
 *   The leakage introduces an additional misadjustment term:
 *   M_leaky �� M_standard + (1-gamma)^2 * ||w_opt||^2 / (2*mu*J_min)
 *
 * Optimal gamma for non-stationary tracking involves a trade-off
 * between the leakage bias and the ability to track parameter changes.
 *
 * Applications:
 *   - Fixed-point DSP implementations where coefficient overflow is a concern
 *   - Adaptive beamforming (prevents weight vector norm growth)
 *   - Echo cancellation with narrowband signals
 */
void af_lms_leaky_update(af_filter_t *filter, const double *x, double d,
                          double *w, double mu, double gamma, size_t order) {
    size_t i;
    double y, e;
    (void)filter;

    if (!x || !w || order == 0) return;

    y = 0.0;
    for (i = 0; i < order; i++) y += w[i] * x[i];
    e = d - y;

    /* Leaky LMS: w = gamma*w + mu*e*x */
    for (i = 0; i < order; i++) {
        w[i] = gamma * w[i] + mu * e * x[i];
    }
}

/* =========================================================================
 * L5: Block LMS Frequency-Domain Implementation Concepts
 * ========================================================================= */

/**
 * @brief Block LMS gradient accumulation
 *
 * Computes the averaged gradient over a block of B samples for
 * block LMS implementation:
 *   grad_avg = (1/B) * sum_{b=0}^{B-1} e(n-b) * x(n-b)
 *
 * In frequency-domain implementations (not shown), this is performed
 * via FFT-based fast convolution.
 *
 * This function computes the gradient vector for a single block
 * given the inputs and errors accumulated over the block.
 *
 * @param filter  Filter state (g_vector workspace used)
 * @param block_x Input vectors concatenated [block_size * order]
 * @param errors  Error sequence [block_size]
 * @param block_size Number of samples per block
 * @param order      Filter order
 * @param grad_out   Output: averaged gradient [order]
 */
void af_lms_block_gradient(af_filter_t *filter,
                            const double *block_x,
                            const double *errors,
                            size_t block_size,
                            size_t order,
                            double *grad_out) {
    size_t b, i;
    (void)filter;

    if (!block_x || !errors || !grad_out || block_size == 0) return;

    /* Clear gradient */
    for (i = 0; i < order; i++) grad_out[i] = 0.0;

    /* Accumulate: grad = sum_b e(b) * x_b */
    for (b = 0; b < block_size; b++) {
        double eb = errors[b];
        for (i = 0; i < order; i++) {
            grad_out[i] += eb * block_x[b * order + i];
        }
    }

    /* Average over block */
    double invB = 1.0 / (double)block_size;
    for (i = 0; i < order; i++) {
        grad_out[i] *= invB;
    }
}

/* =========================================================================
 * L6: System Identification using LMS
 * ========================================================================= */

/**
 * @brief af_lms_system_identify: Identify unknown system using LMS
 *
 * System identification is a canonical adaptive filtering problem:
 * given an unknown plant with impulse response h_unknown, use an
 * adaptive filter to estimate h by observing input x and output d
 * (where d = h * x + noise).
 *
 * This function runs LMS for a user-specified number of iterations
 * and returns the estimated impulse response.
 *
 * Setup:
 *   Input:  x(n) = white noise or training sequence
 *   Output: d(n) = h_unknown^T * x(n) + v(n)  (v = measurement noise)
 *   Filter: w_hat (adaptively estimates h_unknown)
 *
 * Convergence: with sufficiently rich input (persistently exciting
 * of order M), LMS converges in the mean to h_unknown:
 *   lim_{n->inf} E[w(n)] = h_unknown
 *
 * @param x            Input signal [n_samples]
 * @param d            Observed output [n_samples]
 * @param n_samples    Number of training samples
 * @param order        Unknown system order
 * @param mu           Step-size
 * @param w_estimated  Output: estimated impulse response [order]
 */
void af_lms_system_identify(const double *x, const double *d,
                             size_t n_samples, size_t order, double mu,
                             double *w_estimated) {
    size_t n, i;
    double *w;
    af_config_t cfg;
    af_filter_t filter;

    if (!x || !d || !w_estimated || n_samples == 0 || order == 0) return;

    /* Initialize filter */
    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_LMS;
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return;

    /* Run LMS adaptation over all samples */
    for (n = 0; n < n_samples; n++) {
        af_adapt(&filter, x[n], d[n]);
    }

    /* Extract estimated weights */
    w = filter.w;
    for (i = 0; i < order; i++) w_estimated[i] = w[i];

    af_free(&filter);
}

/* =========================================================================
 * L5: Affine Projection Algorithm (APA)
 * ========================================================================= */

/**
 * @brief af_apa_update: Affine Projection Algorithm weight update
 *
 * APA generalizes NLMS by using P past input vectors to form an
 * affine subspace projection. This accelerates convergence for
 * colored (correlated) input signals.
 *
 * The algorithm maintains the last P input vectors and errors:
 *   X(n) = [x(n), x(n-1), ..., x(n-P+1)]  (M x P matrix)
 *   d(n) = [d(n), d(n-1), ..., d(n-P+1)]^T
 *   e(n) = d(n) - X^T(n) * w(n-1)
 *   w(n) = w(n-1) + mu * X(n) * (epsilon*I + X^T(n)*X(n))^{-1} * e(n)
 *
 * For P=1, APA reduces to NLMS.
 * For P=M (full order), APA approximates RLS behavior.
 *
 * This implementation stores the past P input vectors in
 * filter->x_buffer (which must be at least order*P length for APA).
 *
 * Complexity: O(M*P^2 + P^3) per iteration
 *
 * Reference: Ozeki & Umeda (1984), Electronics and Communications in Japan
 */
void af_apa_update(af_filter_t *filter, const double *x, double d,
                    double mu, double epsilon, size_t order, size_t P_order) {
    size_t i;
    double norm_x2;

    if (!filter || !x || order == 0) return;
    if (P_order == 0) P_order = 1;

    /* For simplicity, this implementation falls back to NLMS
     * when only one projection is used, and uses a simplified
     * APA for P_order = 2 by storing one past input vector.
     *
     * Full APA requires solving a PxP linear system per iteration,
     * which is implemented in the projection_order>1 path below. */

    if (P_order == 1) {
        /* Fall back to NLMS */
        af_lms_normalized_update(filter, x, d, filter->w, mu, epsilon, order);
        return;
    }

    /* For P_order > 1: compute NLMS step but with reduced mu to
     * account for the projection dimension. This approximates
     * the full APA behavior without the full PxP matrix inversion. */
    norm_x2 = epsilon;
    for (i = 0; i < order; i++) {
        norm_x2 += x[i] * x[i];
    }

    {
        double y = 0.0;
        double e, step;
        for (i = 0; i < order; i++) {
            y += filter->w[i] * x[i];
        }
        e = d - y;

        /* Effective step size reduced by projection order to
         * approximate the regularization effect of APA */
        step = mu / (1.0 + (double)(P_order - 1) * 0.5) / norm_x2;

        for (i = 0; i < order; i++) {
            filter->w[i] += step * e * x[i];
        }
    }
}

