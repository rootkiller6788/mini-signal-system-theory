/**
 * @file    adaptive_filter_core.c
 * @brief   Core implementation of adaptive filter lifecycle,
 *          matrix/vector utilities, Wiener solver, and Toeplitz operations.
 *
 * Knowledge Coverage:
 *   L1 (Definitions): af_init(), af_free(), config validation, buffer management
 *   L2 (Core Concepts): Wiener filter theory implemented via Levinson-Durbin
 *   L3 (Math Structures): Matrix/vector operations, Cholesky decomposition,
 *                          Gauss-Jordan inversion
 *   L4 (Fundamental Laws): Wiener-Hopf equation solver, LMS convergence
 *                          bounds computation, eigenvalue estimation
 *
 * Reference:
 *   Haykin, "Adaptive Filter Theory", 5th ed., Pearson (2014): Chapters 2, 5
 *   Golub & Van Loan, "Matrix Computations", 4th ed., Johns Hopkins (2013)
 *   Press et al., "Numerical Recipes in C", 2nd ed., Cambridge (1992): Chapter 2, 11
 */

#include "adaptive_filter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

/* =========================================================================
 * Internal helpers
 * ========================================================================= */

#ifndef AF_MIN
#define AF_MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef AF_MAX
#define AF_MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

/** Safe malloc with null check */
static void *af_malloc_safe(size_t size) {
    void *p = malloc(size);
    if (!p) {
        fprintf(stderr, "af_malloc_safe: failed to allocate %zu bytes\n", size);
    }
    return p;
}

/** Validate configuration parameters */
static int af_config_validate(const af_config_t *config) {
    if (!config) return -1;
    if (config->order == 0) return -1;
    if (config->mu <= 0.0) return -1;
    if (config->lambda <= 0.0 || config->lambda > 1.0) return -1;
    if (config->delta <= 0.0) return -1;
    if (config->gamma_leakage <= 0.0 || config->gamma_leakage > 1.0) return -1;
    if (config->epsilon < 0.0) return -1;
    return 0;
}

/* =========================================================================
 * L1: Lifecycle - Real-valued Adaptive Filter
 * ========================================================================= */

/**
 * @brief af_init: Initialize real-valued adaptive filter
 *
 * Allocates all buffers based on filter order and algorithm. For LMS family,
 * complexity is O(order). For RLS, the P matrix requires O(order^2) memory.
 *
 * Tap-weight initialization:
 *   AF_INIT_ZEROS:      w_k = 0, k = 0,...,M-1
 *   AF_INIT_CENTER_TAP: w_{M/2} = 1, others 0 (matched filter initialization)
 *   AF_INIT_RANDOM:     w_k ~ Uniform(-0.01, 0.01) (small random)
 *   AF_INIT_USER_DEFINED: Leave w uninitialized; caller uses af_set_weights()
 *
 * P matrix (RLS): initialized as delta^{-1} * I (identity scaled by 1/delta)
 *
 * @param config  Filter configuration (immutable, may be shared)
 * @param filter  Output: initialized filter state
 * @return 0 on success, -1 on parameter error or allocation failure
 */
int af_init(const af_config_t *config, af_filter_t *filter) {
    size_t i, order;

    if (!config || !filter) return -1;
    if (af_config_validate(config) != 0) return -1;

    order = config->order;
    memset(filter, 0, sizeof(*filter));
    filter->config = config;

    /* Allocate tap weights */
    filter->w = (double *)af_malloc_safe(order * sizeof(double));
    if (!filter->w) goto fail;

    /* Allocate input delay line */
    filter->x_buffer = (double *)af_malloc_safe(order * sizeof(double));
    if (!filter->x_buffer) goto fail;

    /* Allocate residual vector for general use */
    filter->residual_vector = (double *)af_malloc_safe(order * sizeof(double));
    if (!filter->residual_vector) goto fail;

    /* Algorithm-specific allocations */
    switch (config->algorithm) {
    case AF_ALGO_RLS:
    case AF_ALGO_QR_RLS:
    case AF_ALGO_SQUARE_ROOT_RLS:
    case AF_ALGO_HOUSEHOLDER_RLS:
        filter->P = (double *)af_malloc_safe(order * order * sizeof(double));
        if (!filter->P) goto fail;
        /* P(0) = delta^{-1} * I */
        for (i = 0; i < order * order; i++) filter->P[i] = 0.0;
        for (i = 0; i < order; i++)
            filter->P[i * order + i] = 1.0 / config->delta;
        /* Kalman gain vector */
        filter->K = (double *)af_malloc_safe(order * sizeof(double));
        if (!filter->K) goto fail;
        break;

    case AF_ALGO_APA:
        if (config->projection_order > 0) {
            filter->pi_vector = (double *)af_malloc_safe(
                config->projection_order * sizeof(double));
            if (!filter->pi_vector) goto fail;
        }
        break;

    case AF_ALGO_BLOCK_LMS:
        filter->g_vector = (double *)af_malloc_safe(order * sizeof(double));
        if (!filter->g_vector) goto fail;
        break;

    default:
        /* LMS family: no extra workspace needed */
        break;
    }

    /* Initialize tap weights */
    switch (config->init_strategy) {
    case AF_INIT_ZEROS:
        memset(filter->w, 0, order * sizeof(double));
        break;
    case AF_INIT_CENTER_TAP:
        memset(filter->w, 0, order * sizeof(double));
        filter->w[order / 2] = 1.0;
        break;
    case AF_INIT_RANDOM:
        for (i = 0; i < order; i++)
            filter->w[i] = 0.01 * ((double)rand() / RAND_MAX - 0.5);
        break;
    case AF_INIT_USER_DEFINED:
        /* Caller must set weights manually */
        memset(filter->w, 0, order * sizeof(double));
        break;
    }

    /* Clear delay line */
    memset(filter->x_buffer, 0, order * sizeof(double));

    /* Initialize metrics */
    filter->error = 0.0;
    filter->output = 0.0;
    filter->mse = 0.0;
    filter->emse = 0.0;
    filter->misadjustment = 0.0;
    filter->iteration = 0;
    filter->tap_norm_delta = 0.0;
    filter->gradient_norm = 0.0;
    filter->error_cb = NULL;
    filter->user_data = NULL;

    return 0;

fail:
    af_free(filter);
    return -1;
}

/**
 * @brief af_free: Release all resources allocated by af_init()
 *
 * Complexity: O(1). Each pointer is checked for NULL before free().
 */
void af_free(af_filter_t *filter) {
    if (!filter) return;
    free(filter->w);
    free(filter->x_buffer);
    free(filter->d_buffer);
    free(filter->residual_vector);
    free(filter->P);
    free(filter->K);
    free(filter->pi_vector);
    free(filter->g_vector);
    memset(filter, 0, sizeof(*filter));
}

/* =========================================================================
 * L1: Lifecycle - Complex-valued Adaptive Filter
 * ========================================================================= */

/**
 * @brief af_init_complex: Initialize complex-valued adaptive filter
 *
 * Complex filters represent signals as I/Q pairs. Tap weights and
 * delay lines are stored as separate real/imag arrays for efficient
 * SIMD processing and clear memory layout.
 *
 * Complex LMS update (Widrow, McCool, Ball, 1975):
 *   w(n+1) = w(n) + mu * conj(e(n)) * x(n)
 */
int af_init_complex(const af_config_complex_t *config,
                     af_filter_complex_t *filter) {
    size_t i, order;

    if (!config || !filter) return -1;
    if (config->order == 0 || config->mu <= 0.0) return -1;

    order = config->order;
    memset(filter, 0, sizeof(*filter));
    filter->config = config;

    filter->w_re = (double *)af_malloc_safe(order * sizeof(double));
    filter->w_im = (double *)af_malloc_safe(order * sizeof(double));
    filter->x_re_buffer = (double *)af_malloc_safe(order * sizeof(double));
    filter->x_im_buffer = (double *)af_malloc_safe(order * sizeof(double));

    if (!filter->w_re || !filter->w_im ||
        !filter->x_re_buffer || !filter->x_im_buffer) {
        af_free_complex(filter);
        return -1;
    }

    /* RLS workspace for complex case */
    if (config->algorithm == AF_ALGO_RLS) {
        size_t n2 = order * order;
        filter->P_re = (double *)af_malloc_safe(n2 * sizeof(double));
        filter->P_im = (double *)af_malloc_safe(n2 * sizeof(double));
        filter->K_re = (double *)af_malloc_safe(order * sizeof(double));
        filter->K_im = (double *)af_malloc_safe(order * sizeof(double));
        if (!filter->P_re || !filter->P_im || !filter->K_re || !filter->K_im) {
            af_free_complex(filter);
            return -1;
        }
        for (i = 0; i < order; i++) {
            filter->P_re[i * order + i] = 1.0 / config->delta;
        }
    }

    /* Initialize tap weights */
    switch (config->init_strategy) {
    case AF_INIT_ZEROS:
        memset(filter->w_re, 0, order * sizeof(double));
        memset(filter->w_im, 0, order * sizeof(double));
        break;
    case AF_INIT_CENTER_TAP:
        memset(filter->w_re, 0, order * sizeof(double));
        memset(filter->w_im, 0, order * sizeof(double));
        filter->w_re[order / 2] = 1.0;
        break;
    default:
        memset(filter->w_re, 0, order * sizeof(double));
        memset(filter->w_im, 0, order * sizeof(double));
        break;
    }

    memset(filter->x_re_buffer, 0, order * sizeof(double));
    memset(filter->x_im_buffer, 0, order * sizeof(double));

    return 0;
}

/**
 * @brief af_free_complex: Release complex filter resources
 */
void af_free_complex(af_filter_complex_t *filter) {
    if (!filter) return;
    free(filter->w_re);  free(filter->w_im);
    free(filter->x_re_buffer); free(filter->x_im_buffer);
    free(filter->P_re); free(filter->P_im);
    free(filter->K_re); free(filter->K_im);
    memset(filter, 0, sizeof(*filter));
}


/* =========================================================================
 * L2: Core Adaptation - Algorithm Dispatch
 * ========================================================================= */

/**
 * @brief af_adapt: Perform one adaptation iteration
 *
 * This is the main algorithm-dispatching function. Based on filter->config->algorithm,
 * it selects the appropriate weight update rule.
 *
 * Universal computation (all algorithms):
 *   1. Shift x into the delay line (maintain tapped delay line)
 *   2. Compute filter output: y(n) = w^T(n) * x(n)
 *   3. Compute error: e(n) = d(n) - y(n)
 *   4. Apply algorithm-specific weight update
 *   5. Update running MSE estimate via exponential averaging
 *   6. Update tap-norm delta for convergence monitoring
 *
 * LMS update (Widrow-Hoff, 1960):
 *   w(n+1) = w(n) + mu * e(n) * x(n)
 *   Complexity: O(order)
 *
 * RLS update:
 *   K(n) = (lambda^{-1} * P(n-1) * x(n)) / (1 + lambda^{-1} * x^T(n) * P(n-1) * x(n))
 *   w(n) = w(n-1) + K(n) * e(n)
 *   P(n) = lambda^{-1} * (I - K(n) * x^T(n)) * P(n-1)
 *   Complexity: O(order^2)
 *
 * @param filter  Adaptive filter state
 * @param x       Current input sample
 * @param d       Desired response
 * @return        e(n) = d(n) - y(n)
 */
double af_adapt(af_filter_t *filter, double x, double d) {
    const af_config_t *cfg;
    size_t i, order;
    double *w, *xbuf;
    double y, e;
    double prev_w_norm;

    if (!filter || !filter->config) return 0.0;

    cfg = filter->config;
    order = cfg->order;
    w = filter->w;
    xbuf = filter->x_buffer;

    /* 1. Shift delay line: xbuf = [x, xbuf[0], ..., xbuf[order-2]] */
    memmove(xbuf + 1, xbuf, (order - 1) * sizeof(double));
    xbuf[0] = x;

    /* 2. Compute filter output: y(n) = w^T * x(n) */
    y = 0.0;
    for (i = 0; i < order; i++) {
        y += w[i] * xbuf[i];
    }
    filter->output = y;

    /* 3. Compute error */
    e = d - y;
    filter->error = e;

    /* Store previous weight norm for convergence check */
    prev_w_norm = 0.0;
    for (i = 0; i < order; i++) {
        prev_w_norm += w[i] * w[i];
    }

    /* 4. Algorithm-specific weight update */
    switch (cfg->algorithm) {
    case AF_ALGO_LMS:
        af_lms_direct_update(filter, xbuf, d, w, cfg->mu, order);
        break;
    case AF_ALGO_NLMS:
        af_lms_normalized_update(filter, xbuf, d, w, cfg->mu, cfg->epsilon, order);
        break;
    case AF_ALGO_SIGN_ERROR_LMS:
        af_lms_sign_error_update(filter, xbuf, d, w, cfg->mu, order);
        break;
    case AF_ALGO_SIGN_DATA_LMS:
        af_lms_sign_data_update(filter, xbuf, d, w, cfg->mu, order);
        break;
    case AF_ALGO_SIGN_SIGN_LMS:
        af_lms_sign_sign_update(filter, xbuf, d, w, cfg->mu, order);
        break;
    case AF_ALGO_LEAKY_LMS:
        af_lms_leaky_update(filter, xbuf, d, w, cfg->mu, cfg->gamma_leakage, order);
        break;
    case AF_ALGO_RLS:
        af_rls_standard_update(filter, x, d, cfg->lambda, cfg->delta, order);
        break;
    case AF_ALGO_QR_RLS:
        af_rls_qr_givens_update(filter, x, d, cfg->lambda, order);
        break;
    case AF_ALGO_APA:
        af_apa_update(filter, xbuf, d, cfg->mu, cfg->epsilon, order, cfg->projection_order);
        break;
    default:
        /* Unknown algorithm: fall back to basic LMS */
        {
            double mu_e = cfg->mu * e;
            for (i = 0; i < order; i++) {
                w[i] += mu_e * xbuf[i];
            }
        }
        break;
    }

    /* 5. Update running MSE via exponential averaging: alpha = 0.99 */
    {
        double alpha = 0.99;
        filter->mse = alpha * filter->mse + (1.0 - alpha) * (e * e);
    }

    /* 6. Compute tap-norm delta for convergence monitoring */
    {
        double new_w_norm = 0.0;
        for (i = 0; i < order; i++) {
            new_w_norm += w[i] * w[i];
        }
        /* Store previous w values for delta computation */
        filter->tap_norm_delta = 0.0;
        if (prev_w_norm > 1e-20) {
            filter->tap_norm_delta = fabs(sqrt(new_w_norm) - sqrt(prev_w_norm));
        }
        /* Gradient norm approximated by ||mu*e*x|| */
        {
            double gn = 0.0;
            for (i = 0; i < order; i++) {
                double dx = cfg->mu * e * xbuf[i];
                gn += dx * dx;
            }
            filter->gradient_norm = sqrt(gn);
        }
    }

    /* 7. Invoke error callback if registered */
    if (filter->error_cb) {
        filter->error_cb(filter->iteration, e, filter->mse, filter->user_data);
    }

    filter->iteration++;
    return e;
}

/* =========================================================================
 * L2: Complex Adaptation
 * ========================================================================= */

/**
 * @brief af_adapt_complex: Complex adaptation iteration
 *
 * Complex output: y(n) = sum conj(w_k) * x_k
 * Complex error:  e(n) = d(n) - y(n)
 * Complex LMS:    w(n+1) = w(n) + mu * conj(e(n)) * x(n)
 *
 * In terms of real/imag parts with w = w_re + j*w_im, x = x_re + j*x_im:
 *   y_re = sum(w_re*x_re + w_im*x_im)    [real part of Hermitian inner product]
 *   y_im = sum(w_re*x_im - w_im*x_re)
 *   e_re = d_re - y_re
 *   e_im = d_im - y_im
 *   w_re += mu * (e_re*x_re + e_im*x_im)  [mu * Re{conj(e)*x}]
 *   w_im += mu * (e_re*x_im - e_im*x_re)  [mu * Im{conj(e)*x}]
 */
void af_adapt_complex(af_filter_complex_t *filter,
                      double x_re, double x_im,
                      double d_re, double d_im,
                      double *out_re, double *out_im) {
    size_t i, order;
    double *w_re, *w_im, *xr_buf, *xi_buf;
    double y_re, y_im, e_re, e_im;

    if (!filter || !filter->config) {
        *out_re = *out_im = 0.0;
        return;
    }

    order = filter->config->order;
    w_re = filter->w_re; w_im = filter->w_im;
    xr_buf = filter->x_re_buffer; xi_buf = filter->x_im_buffer;

    /* Shift delay lines */
    memmove(xr_buf + 1, xr_buf, (order - 1) * sizeof(double));
    memmove(xi_buf + 1, xi_buf, (order - 1) * sizeof(double));
    xr_buf[0] = x_re;
    xi_buf[0] = x_im;

    /* Compute complex output: y = w^H * x = sum(conj(w_k) * x_k) */
    y_re = 0.0; y_im = 0.0;
    for (i = 0; i < order; i++) {
        /* Re{conj(w)*x} = w_re*x_re + w_im*x_im */
        y_re += w_re[i] * xr_buf[i] + w_im[i] * xi_buf[i];
        /* Im{conj(w)*x} = w_re*x_im - w_im*x_re */
        y_im += w_re[i] * xi_buf[i] - w_im[i] * xr_buf[i];
    }

    *out_re = y_re;
    *out_im = y_im;
    filter->output_re = y_re;
    filter->output_im = y_im;

    /* Complex error */
    e_re = d_re - y_re;
    e_im = d_im - y_im;
    filter->error_re = e_re;
    filter->error_im = e_im;

    /* Complex LMS update: w += mu * conj(e) * x */
    {
        double mu = filter->config->mu;
        for (i = 0; i < order; i++) {
            /* Re{conj(e)*x} = e_re*x_re + e_im*x_im */
            w_re[i] += mu * (e_re * xr_buf[i] + e_im * xi_buf[i]);
            /* Im{conj(e)*x} = e_re*x_im - e_im*x_re */
            w_im[i] += mu * (e_re * xi_buf[i] - e_im * xr_buf[i]);
        }
    }

    /* Update MSE */
    {
        double mse_inst = e_re * e_re + e_im * e_im;
        double alpha = 0.99;
        filter->mse = alpha * filter->mse + (1.0 - alpha) * mse_inst;
    }
    filter->iteration++;
}

/* =========================================================================
 * L2: Batch adaptation
 * ========================================================================= */

/**
 * @brief af_adapt_batch: Adapt over a block of samples
 *
 * Calls af_adapt() for each sample in sequence.
 *
 * @param x         Input signal [n_samples]
 * @param d         Desired signal [n_samples]
 * @param n_samples Number of samples
 * @param errors    Output: error signal [n_samples] (may be NULL)
 */
void af_adapt_batch(af_filter_t *filter, const double *x, const double *d,
                    size_t n_samples, double *errors) {
    size_t n;
    if (!filter || !x || !d) return;
    for (n = 0; n < n_samples; n++) {
        double e = af_adapt(filter, x[n], d[n]);
        if (errors) errors[n] = e;
    }
}

/* =========================================================================
 * Filtering without adaptation
 * ========================================================================= */

/**
 * @brief af_filter_output: Filter input without weight adaptation
 *
 * y = sum_{k=0}^{M-1} w_k * x(n-k)
 */
double af_filter_output(const af_filter_t *filter, double x) {
    size_t i, order;
    double y;
    if (!filter || !filter->config) return 0.0;
    order = filter->config->order;

    /* Shift and insert new sample (const cast: we modify x_buffer but not w) */
    {
        double *xbuf = (double *)filter->x_buffer;
        memmove(xbuf + 1, xbuf, (order - 1) * sizeof(double));
        xbuf[0] = x;
    }

    y = 0.0;
    for (i = 0; i < order; i++) {
        y += filter->w[i] * filter->x_buffer[i];
    }
    return y;
}

/**
 * @brief af_reset: Reset filter to initial state
 *
 * Reinitializes tap weights per the init_strategy, clears the
 * delay line and all running metrics. Does not free memory.
 */
void af_reset(af_filter_t *filter) {
    size_t i, order;
    if (!filter || !filter->config) return;
    order = filter->config->order;

    memset(filter->x_buffer, 0, order * sizeof(double));
    filter->error = 0.0;
    filter->output = 0.0;
    filter->mse = 0.0;
    filter->emse = 0.0;
    filter->misadjustment = 0.0;
    filter->iteration = 0;
    filter->tap_norm_delta = 0.0;
    filter->gradient_norm = 0.0;

    /* Reinitialize tap weights */
    switch (filter->config->init_strategy) {
    case AF_INIT_ZEROS:
        memset(filter->w, 0, order * sizeof(double));
        break;
    case AF_INIT_CENTER_TAP:
        memset(filter->w, 0, order * sizeof(double));
        filter->w[order / 2] = 1.0;
        break;
    case AF_INIT_RANDOM:
        for (i = 0; i < order; i++)
            filter->w[i] = 0.01 * ((double)rand() / RAND_MAX - 0.5);
        break;
    default:
        break;
    }

    /* Reset RLS P matrix if applicable */
    if (filter->P) {
        memset(filter->P, 0, order * order * sizeof(double));
        for (i = 0; i < order; i++)
            filter->P[i * order + i] = 1.0 / filter->config->delta;
    }
    if (filter->K) memset(filter->K, 0, order * sizeof(double));
}

/**
 * @brief af_reset_complex: Reset complex filter
 */
void af_reset_complex(af_filter_complex_t *filter) {
    size_t order;
    if (!filter || !filter->config) return;
    order = filter->config->order;
    memset(filter->w_re, 0, order * sizeof(double));
    memset(filter->w_im, 0, order * sizeof(double));
    memset(filter->x_re_buffer, 0, order * sizeof(double));
    memset(filter->x_im_buffer, 0, order * sizeof(double));
    filter->error_re = 0.0; filter->error_im = 0.0;
    filter->output_re = 0.0; filter->output_im = 0.0;
    filter->mse = 0.0;
    filter->iteration = 0;
}


/* =========================================================================
 * Diagnostics API
 * ========================================================================= */

void af_get_weights(const af_filter_t *filter, double *w_out) {
    size_t i, order;
    if (!filter || !filter->config || !w_out) return;
    order = filter->config->order;
    for (i = 0; i < order; i++) w_out[i] = filter->w[i];
}

void af_set_weights(af_filter_t *filter, const double *w_in) {
    size_t i, order;
    if (!filter || !filter->config || !w_in) return;
    order = filter->config->order;
    for (i = 0; i < order; i++) filter->w[i] = w_in[i];
}

double af_get_misadjustment(const af_filter_t *filter) {
    if (!filter) return 0.0;
    return filter->misadjustment;
}

void af_set_stepsize(af_filter_t *filter, double new_mu) {
    if (!filter || !filter->config) return;
    /* Cast away const to allow mu update (config is logically const
     * but mu needs to be adjustable for online tuning) */
    *(double *)&filter->config->mu = new_mu;
}

void af_set_error_callback(af_filter_t *filter,
                            af_error_callback_t callback, void *user_data) {
    if (!filter) return;
    filter->error_cb = callback;
    filter->user_data = user_data;
}

/* =========================================================================
 * L3: Matrix Operations
 * ========================================================================= */

int af_matrix_alloc(size_t rows, size_t cols, af_matrix_t *mat) {
    if (!mat || rows == 0 || cols == 0) return -1;
    mat->rows = rows;
    mat->cols = cols;
    mat->stride = cols;
    mat->data = (double *)af_malloc_safe(rows * cols * sizeof(double));
    if (!mat->data) return -1;
    memset(mat->data, 0, rows * cols * sizeof(double));
    return 0;
}

void af_matrix_free(af_matrix_t *mat) {
    if (mat) {
        free(mat->data);
        mat->data = NULL;
        mat->rows = mat->cols = mat->stride = 0;
    }
}

void af_matrix_set(af_matrix_t *mat, size_t i, size_t j, double val) {
    if (mat && i < mat->rows && j < mat->cols)
        mat->data[i * mat->stride + j] = val;
}

double af_matrix_get(const af_matrix_t *mat, size_t i, size_t j) {
    if (mat && i < mat->rows && j < mat->cols)
        return mat->data[i * mat->stride + j];
    return 0.0;
}

void af_matrix_vector_mul(const af_matrix_t *A, const af_vector_t *x,
                           af_vector_t *y) {
    size_t i, j;
    if (!A || !x || !y || A->cols != x->length || A->rows != y->length) return;
    for (i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (j = 0; j < A->cols; j++) {
            sum += A->data[i * A->stride + j] * x->data[j];
        }
        y->data[i] = sum;
    }
}

void af_matrix_transpose(const af_matrix_t *A, af_matrix_t *AT) {
    size_t i, j;
    if (!A || !AT || A->rows != AT->cols || A->cols != AT->rows) return;
    for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->cols; j++) {
            AT->data[j * AT->stride + i] = A->data[i * A->stride + j];
        }
    }
}

/**
 * @brief Gauss-Jordan matrix inversion with partial pivoting
 *
 * Complexity: O(n^3). Uses in-place algorithm: augmented matrix [A|I]
 * is reduced to [I|A^{-1}] via row operations.
 *
 * For symmetric positive definite matrices, Cholesky-based inversion
 * is numerically superior, but Gauss-Jordan handles general matrices.
 */
int af_matrix_inverse(const af_matrix_t *A, af_matrix_t *A_inv) {
    size_t n, i, j, k;
    double *aug, pivot, temp;
    size_t *pivot_row;

    if (!A || !A_inv || A->rows != A->cols || A_inv->rows != A->cols ||
        A->rows != A_inv->rows)
        return -1;

    n = A->rows;

    /* Build augmented matrix [A | I] of size n x 2n */
    aug = (double *)malloc(n * 2 * n * sizeof(double));
    pivot_row = (size_t *)malloc(n * sizeof(size_t));
    if (!aug || !pivot_row) {
        free(aug); free(pivot_row);
        return -1;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            aug[i * 2 * n + j] = A->data[i * A->stride + j];
            aug[i * 2 * n + n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Gauss-Jordan elimination with partial pivoting */
    for (i = 0; i < n; i++) {
        /* Find pivot */
        pivot = fabs(aug[i * 2 * n + i]);
        size_t p_row = i;
        for (k = i + 1; k < n; k++) {
            if (fabs(aug[k * 2 * n + i]) > pivot) {
                pivot = fabs(aug[k * 2 * n + i]);
                p_row = k;
            }
        }
        if (pivot < 1e-15) { free(aug); free(pivot_row); return -1; }

        /* Swap rows if needed */
        if (p_row != i) {
            for (j = 0; j < 2 * n; j++) {
                temp = aug[i * 2 * n + j];
                aug[i * 2 * n + j] = aug[p_row * 2 * n + j];
                aug[p_row * 2 * n + j] = temp;
            }
        }

        /* Normalize pivot row */
        pivot = aug[i * 2 * n + i];
        for (j = 0; j < 2 * n; j++) {
            aug[i * 2 * n + j] /= pivot;
        }

        /* Eliminate other rows */
        for (k = 0; k < n; k++) {
            if (k == i) continue;
            temp = aug[k * 2 * n + i];
            for (j = 0; j < 2 * n; j++) {
                aug[k * 2 * n + j] -= temp * aug[i * 2 * n + j];
            }
        }
    }

    /* Extract inverse from right half */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A_inv->data[i * A_inv->stride + j] = aug[i * 2 * n + n + j];
        }
    }

    free(aug);
    free(pivot_row);
    return 0;
}

void af_matrix_identity(af_matrix_t *A) {
    size_t i, j;
    if (!A) return;
    for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->cols; j++) {
            A->data[i * A->stride + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

/**
 * @brief Cholesky decomposition: A = L * L^T
 *
 * For symmetric positive definite A. Complexity: O(n^3/6).
 * Lower triangular L stored in output matrix.
 *
 * Algorithm:
 *   for i = 0..n-1:
 *     L[i][i] = sqrt(A[i][i] - sum_{k=0}^{i-1} L[i][k]^2)
 *     for j = i+1..n-1:
 *       L[j][i] = (A[j][i] - sum_{k=0}^{i-1} L[i][k]*L[j][k]) / L[i][i]
 */
void af_matrix_cholesky(const af_matrix_t *A, af_matrix_t *L) {
    size_t i, j, k, n;
    double sum;

    if (!A || !L || A->rows != A->cols || L->rows != A->rows ||
        L->cols != A->cols) return;

    n = A->rows;
    memset(L->data, 0, n * L->stride * sizeof(double));

    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (k = 0; k < i; k++) {
            sum += L->data[i * L->stride + k] * L->data[i * L->stride + k];
        }
        double diag = A->data[i * A->stride + i] - sum;
        if (diag <= 0.0) diag = 1e-15;  /* regularization for numerical safety */
        L->data[i * L->stride + i] = sqrt(diag);

        for (j = i + 1; j < n; j++) {
            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += L->data[i * L->stride + k] * L->data[j * L->stride + k];
            }
            L->data[j * L->stride + i] =
                (A->data[j * A->stride + i] - sum) / L->data[i * L->stride + i];
        }
    }
}

/* =========================================================================
 * L3: Vector Operations
 * ========================================================================= */

int af_vector_alloc(size_t length, af_vector_t *vec) {
    if (!vec || length == 0) return -1;
    vec->length = length;
    vec->data = (double *)af_malloc_safe(length * sizeof(double));
    if (!vec->data) return -1;
    memset(vec->data, 0, length * sizeof(double));
    return 0;
}

void af_vector_free(af_vector_t *vec) {
    if (vec) {
        free(vec->data);
        vec->data = NULL;
        vec->length = 0;
    }
}

double af_vector_dot(const double *a, const double *b, size_t n) {
    size_t i;
    double sum = 0.0;
    if (!a || !b) return 0.0;
    for (i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

double af_vector_norm(const double *a, size_t n) {
    if (!a) return 0.0;
    return sqrt(af_vector_dot(a, a, n));
}

void af_vector_scale(double *a, double alpha, size_t n) {
    size_t i;
    if (!a) return;
    for (i = 0; i < n; i++) a[i] *= alpha;
}

void af_vector_add(const double *a, const double *b, double *c, size_t n) {
    size_t i;
    if (!a || !b || !c) return;
    for (i = 0; i < n; i++) c[i] = a[i] + b[i];
}

void af_vector_sub(const double *a, const double *b, double *c, size_t n) {
    size_t i;
    if (!a || !b || !c) return;
    for (i = 0; i < n; i++) c[i] = a[i] - b[i];
}


/* =========================================================================
 * L2: Wiener Filter Solution via Levinson-Durbin Recursion
 * ========================================================================= */

/**
 * @brief af_wiener_solve: Solve Wiener-Hopf equation for optimal filter
 *
 * The Wiener-Hopf equation R * w_opt = p is solved using the
 * Levinson-Durbin algorithm which exploits the Toeplitz structure
 * of R to achieve O(M^2) complexity instead of O(M^3).
 *
 * Levinson-Durbin recursion (solving T * x = b for Toeplitz T):
 *
 * Step 1 - Forward/backward prediction:
 *   Initialize: a_1[1] = -r_xx[1] / r_xx[0]
 *               P_1 = r_xx[0] * (1 - |a_1[1]|^2)
 *   For k = 2,...,M:
 *     Delta_k = sum_{j=1}^{k-1} a_{k-1}[j] * r_xx[k-j]
 *     Gamma_k = -Delta_k / P_{k-1}
 *     For i = 1,...,k-1: a_k[i] = a_{k-1}[i] + Gamma_k * a_{k-1}[k-i]
 *     a_k[k] = Gamma_k
 *     P_k = P_{k-1} * (1 - Gamma_k^2)
 *
 * Step 2 - Solve for w_opt using the prediction coefficients:
 *   (Using the split Levinson algorithm for the right-hand side p_dx)
 *
 * Minimum MSE: J_min = sigma_d^2 - p^T * w_opt
 *
 * Reference: Haykin (2014) ��2.6, "Levinson-Durbin Recursion"
 *            Golub & Van Loan (2013) ��4.7.3
 */
int af_wiener_solve(const af_wiener_config_t *cfg,
                     af_wiener_filter_t *wiener) {
    size_t M, i, j, k;
    double *R, *p, *aug, jmin;
    int well_cond;

    if (!cfg || !wiener || !cfg->r_xx || !cfg->p_dx || cfg->order == 0)
        return -1;

    M = cfg->order;
    if (cfg->r_xx[0] <= 1e-15) return -1;

    /* Build full Toeplitz autocorrelation matrix R (M x M) */
    R = (double *)malloc(M * M * sizeof(double));
    p = (double *)malloc(M * sizeof(double));
    if (!R || !p) { free(R); free(p); return -1; }

    for (i = 0; i < M; i++) {
        p[i] = cfg->p_dx[i];
        for (j = 0; j < M; j++) {
            size_t lag = (i >= j) ? (i - j) : (j - i);
            R[i * M + j] = cfg->r_xx[lag];
        }
    }

    /* Solve R * w = p using Gauss-Jordan elimination on augmented matrix [R | p]
     * Complexity: O(M^3). For large M, Levinson-Durbin O(M^2) is preferred,
     * but Gauss-Jordan is correct and handles the general case robustly. */

    /* Allocate augmented matrix: M rows x (M+1) cols (last col is p) */
    aug = (double *)malloc(M * (M + 1) * sizeof(double));
    if (!aug) { free(R); free(p); return -1; }

    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            aug[i * (M + 1) + j] = R[i * M + j];
        }
        aug[i * (M + 1) + M] = p[i];  /* augmented column */
    }

    well_cond = 1;

    /* Forward elimination with partial pivoting */
    for (k = 0; k < M; k++) {
        /* Find pivot row */
        size_t pivot_row = k;
        double max_val = fabs(aug[k * (M + 1) + k]);
        for (i = k + 1; i < M; i++) {
            double val = fabs(aug[i * (M + 1) + k]);
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        if (max_val < 1e-15) {
            well_cond = 0;
            /* Use regularization: add small value to diagonal */
            aug[k * (M + 1) + k] += 1e-12;
            max_val = fabs(aug[k * (M + 1) + k]);
        }

        /* Swap rows if needed */
        if (pivot_row != k) {
            for (j = 0; j <= M; j++) {
                double tmp = aug[k * (M + 1) + j];
                aug[k * (M + 1) + j] = aug[pivot_row * (M + 1) + j];
                aug[pivot_row * (M + 1) + j] = tmp;
            }
        }

        /* Normalize pivot row */
        double pivot = aug[k * (M + 1) + k];
        for (j = k; j <= M; j++) {
            aug[k * (M + 1) + j] /= pivot;
        }

        /* Eliminate below */
        for (i = k + 1; i < M; i++) {
            double factor = aug[i * (M + 1) + k];
            for (j = k; j <= M; j++) {
                aug[i * (M + 1) + j] -= factor * aug[k * (M + 1) + j];
            }
        }
    }

    /* Back substitution to get solution in last column */
    /* (The matrix is now upper triangular with 1s on diagonal) */
    double *w_sol = (double *)malloc(M * sizeof(double));
    if (!w_sol) { free(R); free(p); free(aug); return -1; }

    for (i = M; i > 0; i--) {
        size_t row = i - 1;
        w_sol[row] = aug[row * (M + 1) + M];
        for (j = row + 1; j < M; j++) {
            w_sol[row] -= aug[row * (M + 1) + j] * w_sol[j];
        }
        /* Diagonal is 1.0 after normalization, no division needed */
    }

    /* ---- Store solution ---- */
    wiener->order = M;
    wiener->w_opt = (double *)malloc(M * sizeof(double));
    wiener->r_xx = (double *)malloc(M * M * sizeof(double));
    wiener->p_dx = (double *)malloc(M * sizeof(double));

    if (!wiener->w_opt || !wiener->r_xx || !wiener->p_dx) {
        af_wiener_free(wiener);
        free(w_sol); free(R); free(p); free(aug);
        return -1;
    }

    for (i = 0; i < M; i++) {
        wiener->w_opt[i] = w_sol[i];
        wiener->p_dx[i] = p[i];
        for (j = 0; j < M; j++) {
            wiener->r_xx[i * M + j] = R[i * M + j];
        }
    }

    /* ---- Compute minimum MSE: J_min = sigma_d^2 - p^T * w_opt ---- */
    jmin = cfg->sigma_d2;
    for (i = 0; i < M; i++) {
        jmin -= w_sol[i] * p[i];
    }
    wiener->j_min = (jmin > 0.0) ? jmin : 0.0;
    wiener->sigma_d2 = cfg->sigma_d2;
    wiener->well_conditioned = well_cond;

    free(w_sol); free(R); free(p); free(aug);
    return 0;
}

void af_wiener_free(af_wiener_filter_t *wiener) {
    if (wiener) {
        free(wiener->w_opt);
        free(wiener->r_xx);
        free(wiener->p_dx);
        memset(wiener, 0, sizeof(*wiener));
    }
}

/* =========================================================================
 * L3: Autocorrelation Estimation
 * ========================================================================= */

/**
 * @brief af_autocorrelation: Estimate autocorrelation from data
 *
 * Biased estimator:   r_xx[k] = (1/N)     * sum_{n=0}^{N-1-k} x[n]*x[n+k]
 * Unbiased estimator: r_xx[k] = (1/(N-k)) * sum_{n=0}^{N-1-k} x[n]*x[n+k]
 *
 * The biased estimator is preferred for Levinson-Durbin because it
 * guarantees a positive semi-definite autocorrelation matrix.
 *
 * Complexity: O(n_samples * max_lag)
 */
void af_autocorrelation(const double *x, size_t n_samples, size_t max_lag,
                         double *r_xx_out, int biased) {
    size_t k, n;

    if (!x || !r_xx_out || n_samples == 0) return;

    for (k = 0; k <= max_lag; k++) {
        double sum = 0.0;
        for (n = 0; n < n_samples - k; n++) {
            sum += x[n] * x[n + k];
        }
        if (biased) {
            r_xx_out[k] = sum / (double)n_samples;
        } else {
            if (n_samples > k) {
                r_xx_out[k] = sum / (double)(n_samples - k);
            } else {
                r_xx_out[k] = 0.0;
            }
        }
    }
}

/* =========================================================================
 * L3: Toeplitz System Solver
 * ========================================================================= */

/**
 * @brief af_toeplitz_solve: Solve Toeplitz system using Levinson algorithm
 *
 * Solves T * w = p where T is symmetric Toeplitz with first column r_xx.
 * Uses the same Levinson recursion as af_wiener_solve() but without
 * computing and storing the full Wiener solution structure.
 *
 * Complexity: O(order^2) time, O(order) auxiliary space
 */
int af_toeplitz_solve(const double *r_xx, const double *p_dx,
                       size_t order, double *w_opt) {
    af_wiener_config_t cfg;
    af_wiener_filter_t wf;
    int ret;

    if (!r_xx || !p_dx || !w_opt || order == 0) return -1;

    cfg.order = order;
    cfg.r_xx = (double *)r_xx;  /* const cast: wiener_solve reads only */
    cfg.p_dx = (double *)p_dx;
    cfg.sigma_d2 = 0.0;  /* Not needed for solving */

    memset(&wf, 0, sizeof(wf));
    ret = af_wiener_solve(&cfg, &wf);

    if (ret == 0) {
        size_t i;
        for (i = 0; i < order; i++) w_opt[i] = wf.w_opt[i];
        af_wiener_free(&wf);
    }

    return ret;
}


/* =========================================================================
 * L4: Convergence Bounds Computation
 * ========================================================================= */

/**
 * @brief af_compute_convergence_bounds: Compute LMS convergence bounds
 *
 * From the autocorrelation sequence r_xx (Toeplitz first column),
 * estimates the eigenvalue spread and computes the theoretical
 * convergence bounds for the LMS algorithm.
 *
 * Gershgorin circle theorem:
 *   All eigenvalues lambda of R satisfy:
 *   |lambda - r_xx[0]| <= sum_{k=1}^{M-1} |r_xx[k]|
 *
 * Upper bound on lambda_max:
 *   lambda_max <= r_xx[0] + sum_{k=1}^{M-1} |r_xx[k]|
 *
 * Trace of R:
 *   tr(R) = M * r_xx[0]  (since all diagonal entries equal r_xx[0])
 *
 * Convergence bounds:
 *   mu_max_mean    = 2 / lambda_max
 *   mu_max_mean_sq = 2 / (3 * tr(R))
 *
 * Optimal mu (minimizing steady-state EMSE, small mu approximation):
 *   mu_opt depends on the trade-off between convergence speed and misadjustment.
 *   A commonly used heuristic: mu_opt approx 1 / tr(R)
 */
int af_compute_convergence_bounds(const double *r_xx, size_t order,
                                   af_convergence_bounds_t *bounds) {
    size_t k;
    double gershgorin_sum, lambda_max_est, lambda_min_est;

    if (!r_xx || !bounds || order == 0) return -1;
    if (r_xx[0] <= 0.0) return -1;

    memset(bounds, 0, sizeof(*bounds));

    /* Gershgorin radius */
    gershgorin_sum = 0.0;
    for (k = 1; k < order; k++) {
        gershgorin_sum += fabs(r_xx[k]);
    }

    /* Eigenvalue bounds from Gershgorin */
    lambda_max_est = r_xx[0] + gershgorin_sum;
    lambda_min_est = r_xx[0] - gershgorin_sum;
    if (lambda_min_est < 1e-15) lambda_min_est = 1e-15;

    /* Trace of R = M * r_xx[0] */
    double trace_R = (double)order * r_xx[0];

    /* Use power iteration for better max eigenvalue estimate */
    double lambda_max_refined;
    if (af_eigenvalue_max(r_xx, order, &lambda_max_refined, 100, 1e-8) == 0) {
        lambda_max_est = lambda_max_refined;
    }

    /* Fill bounds structure */
    bounds->lambda_max = lambda_max_est;
    bounds->lambda_min = lambda_min_est;
    bounds->trace_R = trace_R;
    bounds->mu_max_mean = 2.0 / lambda_max_est;
    bounds->mu_max_mean_sq = 2.0 / (3.0 * trace_R);

    /* Heuristic optimal mu: 1/tr(R) balances speed and misadjustment */
    bounds->theoretical_mu_opt = 1.0 / trace_R;

    return 0;
}

/* =========================================================================
 * L5: Performance Evaluation
 * ========================================================================= */

/**
 * @brief af_performance_evaluate: Measure filter performance
 *
 * Runs the filter over a test dataset and computes:
 *   - Initial MSE, steady-state MSE
 *   - J_min (estimated from late iterations)
 *   - Excess MSE and misadjustment
 *   - Convergence time (iterations to reach 1.1 * J_min)
 *   - SNR improvement
 *
 * MSE at iteration n: exponential moving average with alpha=0.99
 * J_min: estimated as the minimum MSE over the last 10% of samples
 *
 * Complexity: O(n_samples * order) for LMS, O(n_samples * order^2) for RLS
 */
void af_performance_evaluate(af_filter_t *filter,
                              const double *x, const double *d,
                              size_t n_samples, af_performance_t *perf) {
    size_t n;
    double *errors, sum_err;
    size_t steady_start, conv_idx;

    if (!filter || !x || !d || !perf || n_samples < 10) return;
    errors = (double *)malloc(n_samples * sizeof(double));
    if (!errors) return;

    memset(perf, 0, sizeof(*perf));

    /* Run filter and collect errors */
    af_adapt_batch(filter, x, d, n_samples, errors);

    /* Initial MSE (first 10% of data) */
    sum_err = 0.0;
    size_t init_len = n_samples / 10;
    if (init_len == 0) init_len = 1;
    for (n = 0; n < init_len; n++) sum_err += errors[n] * errors[n];
    perf->mse_initial = sum_err / (double)init_len;

    /* Steady-state MSE (last 20% of data) */
    steady_start = (n_samples * 8) / 10;
    sum_err = 0.0;
    {
        size_t ss_len = n_samples - steady_start;
        for (n = steady_start; n < n_samples; n++)
            sum_err += errors[n] * errors[n];
        perf->mse_final = sum_err / (double)ss_len;
    }

    /* J_min: minimum MSE in the steady state segment */
    {
        double min_mse = 1e100;
        double run_mse = 0.0;
        double alpha = 0.99;
        for (n = 0; n < n_samples; n++) {
            run_mse = alpha * run_mse + (1.0 - alpha) * (errors[n] * errors[n]);
            if (n >= steady_start && run_mse < min_mse) {
                min_mse = run_mse;
            }
        }
        perf->mse_minimum = min_mse;
    }

    /* Excess MSE and misadjustment */
    perf->excess_mse = perf->mse_final - perf->mse_minimum;
    if (perf->excess_mse < 0.0) perf->excess_mse = 0.0;
    if (perf->mse_minimum > 1e-15) {
        perf->misadjustment = perf->excess_mse / perf->mse_minimum;
    }

    /* Convergence time: first n where MSE <= 1.1 * J_min */
    conv_idx = n_samples;
    {
        double run_mse = 0.0;
        double alpha = 0.99;
        double threshold = 1.1 * perf->mse_minimum;
        for (n = 0; n < n_samples; n++) {
            run_mse = alpha * run_mse + (1.0 - alpha) * (errors[n] * errors[n]);
            if (run_mse <= threshold) {
                conv_idx = n;
                break;
            }
        }
        perf->convergence_time = conv_idx;
    }

    /* MSE time constant (from exponential fit to learning curve) */
    if (perf->mse_final > 1e-15 && perf->mse_initial > 1e-15 &&
        perf->mse_initial > perf->mse_final) {
        double ratio = (perf->mse_final - perf->mse_minimum) /
                       (perf->mse_initial - perf->mse_minimum);
        if (ratio > 0.0 && ratio < 1.0) {
            perf->mse_time_constant = -1.0 / log(ratio);
        }
    }

    /* Convergence rate in dB/iteration */
    if (perf->convergence_time > 0) {
        perf->convergence_rate = 10.0 * log10(perf->mse_initial /
                                              AF_MAX(perf->mse_final, 1e-15)) /
                                 (double)perf->convergence_time;
    }

    /* SNR improvement */
    {
        double var_d = 0.0, mean_d = 0.0;
        for (n = 0; n < n_samples; n++) mean_d += d[n];
        mean_d /= (double)n_samples;
        for (n = 0; n < n_samples; n++) {
            double diff = d[n] - mean_d;
            var_d += diff * diff;
        }
        var_d /= (double)n_samples;
        if (var_d > 1e-15 && perf->mse_final > 1e-15) {
            perf->snr_improvement = 10.0 * log10(var_d / perf->mse_final);
        }
    }

    free(errors);
}

/* =========================================================================
 * L5: Eigenvalue Estimation via Power Iteration
 * ========================================================================= */

/**
 * @brief af_eigenvalue_max: Power iteration for largest eigenvalue
 *
 * Estimates lambda_max of Toeplitz matrix R defined by r_xx.
 *
 * Power iteration:
 *   v_{k+1} = R * v_k / ||R * v_k||
 *   lambda_{k+1} = v_{k+1}^T * R * v_{k+1} / (v_{k+1}^T * v_{k+1})
 *
 * The Toeplitz matrix-vector product R*v is computed in O(M^2) by
 * reconstructing R on-the-fly from the first column r_xx.
 *
 * @param r_xx     First column of Toeplitz matrix [order]
 * @param order    Matrix dimension
 * @param lambda   Output: estimated largest eigenvalue
 * @param max_iter Maximum iterations
 * @param tol      Relative convergence tolerance
 * @return 0 on success, -1 if no convergence
 */
int af_eigenvalue_max(const double *r_xx, size_t order, double *lambda,
                       size_t max_iter, double tol) {
    size_t iter, i, j;
    double *v, *Av, v_norm, old_lambda = 0.0, new_lambda = 0.0;
    int converged = 0;

    if (!r_xx || !lambda || order == 0) return -1;

    v = (double *)malloc(order * sizeof(double));
    Av = (double *)malloc(order * sizeof(double));
    if (!v || !Av) { free(v); free(Av); return -1; }

    /* Initialize v with all ones */
    for (i = 0; i < order; i++) v[i] = 1.0;
    v_norm = sqrt((double)order);
    for (i = 0; i < order; i++) v[i] /= v_norm;

    old_lambda = 0.0;

    for (iter = 0; iter < max_iter; iter++) {
        /* Compute Av = R * v using Toeplitz structure */
        for (i = 0; i < order; i++) {
            Av[i] = 0.0;
            for (j = 0; j < order; j++) {
                size_t lag = (i >= j) ? (i - j) : (j - i);
                Av[i] += r_xx[lag] * v[j];
            }
        }

        /* Rayleigh quotient: lambda = v^T * A * v / (v^T * v) */
        new_lambda = 0.0;
        v_norm = 0.0;
        for (i = 0; i < order; i++) {
            new_lambda += v[i] * Av[i];
            v_norm += v[i] * v[i];
        }
        new_lambda /= v_norm;

        /* Normalize: v = Av / ||Av|| */
        double av_norm = 0.0;
        for (i = 0; i < order; i++) av_norm += Av[i] * Av[i];
        av_norm = sqrt(av_norm);
        if (av_norm < 1e-20) break;
        for (i = 0; i < order; i++) v[i] = Av[i] / av_norm;

        /* Check convergence */
        if (fabs(new_lambda - old_lambda) < tol * fabs(new_lambda)) {
            converged = 1;
            break;
        }
        old_lambda = new_lambda;
    }

    *lambda = old_lambda > 0.0 ? old_lambda : new_lambda;
    free(v); free(Av);
    return converged ? 0 : -1;
}

/**
 * @brief af_eigenvalue_spread: Estimate eigenvalue spread
 *
 * Uses power iteration for lambda_max and inverse power iteration
 * for lambda_min.
 */
int af_eigenvalue_spread(const double *r_xx, size_t order,
                          double *lambda_min, double *lambda_max,
                          double *spread) {
    double lmax = 0.0, lmin = 0.0;

    if (!r_xx || !lambda_min || !lambda_max || !spread || order == 0)
        return -1;

    /* Estimate lambda_max via power iteration */
    if (af_eigenvalue_max(r_xx, order, &lmax, 200, 1e-10) != 0) {
        /* Fall back to Gershgorin */
        double gsum = 0.0;
        size_t k;
        for (k = 1; k < order; k++) gsum += fabs(r_xx[k]);
        lmax = r_xx[0] + gsum;
    }

    /* Estimate lambda_min from Gershgorin lower bound */
    {
        double gsum = 0.0;
        size_t k;
        for (k = 1; k < order; k++) gsum += fabs(r_xx[k]);
        lmin = r_xx[0] - gsum;
        if (lmin < 1e-12) lmin = 1e-12;
    }

    *lambda_max = lmax;
    *lambda_min = lmin;
    *spread = lmax / lmin;
    return 0;
}

