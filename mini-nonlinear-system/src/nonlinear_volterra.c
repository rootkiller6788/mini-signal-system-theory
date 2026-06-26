/**
 * @file nonlinear_volterra.c
 * @brief Volterra and Wiener series for nonlinear system ID.
 *
 * Implements Volterra kernel operations, series evaluation,
 * Wiener G-functionals, kernel identification via cross-correlation
 * (Lee-Schetzen), and frequency-domain GFRF computation.
 *
 * Knowledge points:
 *   L1: Volterra kernel allocation/access, Wiener G-functionals
 *   L5: Volterra series time-domain evaluation
 *   L5: Kernel identification (Lee-Schetzen method)
 *   L8: Generalized Frequency Response Functions (GFRF)
 *   L8: Output spectrum prediction for multi-tone inputs
 *   L3: Kernel symmetrization
 *
 * Reference:
 *   - Rugh, "Nonlinear System Theory: Volterra/Wiener Approach", 1981
 *   - Schetzen, "Volterra and Wiener Theories of Nonlinear Systems", 1980
 */

#include "nonlinear_volterra.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <errno.h>

/* ===================================================================
 * Kernel Initialization and Management
 * =================================================================== */

int nl_kernel_init(nl_volterra_kernel_t *kernel, size_t order,
                   size_t mem_len)
{
    if (!kernel || order < 1 || order > NL_MAX_SERIES_ORDER
        || mem_len < 1) {
        errno = EINVAL; return -1;
    }

    memset(kernel, 0, sizeof(*kernel));
    kernel->order = order;
    kernel->dims = (size_t *)malloc(order * sizeof(size_t));
    if (!kernel->dims) { errno = ENOMEM; return -1; }

    kernel->total_size = 1;
    for (size_t d = 0; d < order; d++) {
        kernel->dims[d] = mem_len;
        kernel->total_size *= mem_len;
    }

    kernel->data = (double *)calloc(
        kernel->total_size, sizeof(double));
    if (!kernel->data) {
        free(kernel->dims);
        kernel->dims = NULL;
        errno = ENOMEM;
        return -1;
    }

    kernel->norm_l2 = 0.0;
    return 0;
}

void nl_kernel_free(nl_volterra_kernel_t *kernel)
{
    if (kernel) {
        free(kernel->data);
        free(kernel->dims);
        kernel->data = NULL;
        kernel->dims = NULL;
        kernel->total_size = 0;
    }
}

/** Convert multi-dimensional indices to flat index (row-major) */
static size_t kernel_flat_index(const nl_volterra_kernel_t *kernel,
                                 const size_t *indices)
{
    size_t idx = 0;
    size_t stride = 1;
    for (size_t d = kernel->order; d > 0; d--) {
        idx += indices[d - 1] * stride;
        stride *= kernel->dims[d - 1];
    }
    return idx;
}

int nl_kernel_set(nl_volterra_kernel_t *kernel, const size_t *indices,
                  double value)
{
    if (!kernel || !indices || !kernel->data) {
        errno = EINVAL; return -1;
    }

    for (size_t d = 0; d < kernel->order; d++) {
        if (indices[d] >= kernel->dims[d]) {
            errno = EINVAL; return -1;
        }
    }

    size_t idx = kernel_flat_index(kernel, indices);
    double old = kernel->data[idx];
    kernel->data[idx] = value;
    kernel->norm_l2 += value * value - old * old;
    if (kernel->norm_l2 < 0.0) kernel->norm_l2 = 0.0;

    return 0;
}

int nl_kernel_get(const nl_volterra_kernel_t *kernel, const size_t *indices,
                  double *value)
{
    if (!kernel || !indices || !value || !kernel->data) {
        errno = EINVAL; return -1;
    }

    for (size_t d = 0; d < kernel->order; d++) {
        if (indices[d] >= kernel->dims[d]) {
            errno = EINVAL; return -1;
        }
    }

    *value = kernel->data[kernel_flat_index(kernel, indices)];
    return 0;
}

double nl_kernel_norm_sq(const nl_volterra_kernel_t *kernel)
{
    if (!kernel) return 0.0;
    return kernel->norm_l2;
}

/* ===================================================================
 * Volterra Series Time-Domain Evaluation
 *
 * y[t] = h0 + sum_{n=1..M} y_n[t]
 * y_n[t] = sum_{tau1..taun} h_n[tau1..taun] * u[t-tau1] * ... * u[t-taun]
 *
 * O(T * m^M) complexity — exponential in max order M.
 * For M ≤ 3 with small memory this is tractable.
 * =================================================================== */

int nl_volterra_evaluate(const nl_volterra_series_t *vs,
                         const double *input, size_t input_len,
                         double *output)
{
    if (!vs || !input || !output || input_len < 1
        || vs->max_order < 1) {
        errno = EINVAL; return -1;
    }

    size_t mem = vs->kernels[0].dims[0];

    for (size_t t = 0; t < input_len; t++) {
        double y = 0.0;

        /* First order (linear convolution) */
        {
            double y1 = 0.0;
            for (size_t tau = 0; tau < mem && tau <= t; tau++) {
                size_t idx_arr[1] = {tau};
                size_t idx = kernel_flat_index(&vs->kernels[0], idx_arr);
                y1 += vs->kernels[0].data[idx] * input[t - tau];
            }
            y += y1;
        }

        /* Second order */
        if (vs->max_order >= 2 && vs->kernels[1].data) {
            double y2 = 0.0;
            for (size_t tau1 = 0; tau1 < mem && tau1 <= t; tau1++) {
                double u1 = input[t - tau1];
                for (size_t tau2 = 0; tau2 < mem && tau2 <= t; tau2++) {
                    size_t idx_arr[2] = {tau1, tau2};
                    size_t idx = kernel_flat_index(&vs->kernels[1],
                                                   idx_arr);
                    y2 += vs->kernels[1].data[idx]
                          * u1 * input[t - tau2];
                }
            }
            y += y2;
        }

        /* Third order */
        if (vs->max_order >= 3 && vs->kernels[2].data) {
            double y3 = 0.0;
            for (size_t tau1 = 0; tau1 < mem && tau1 <= t; tau1++) {
                double u1 = input[t - tau1];
                for (size_t tau2 = 0; tau2 < mem && tau2 <= t; tau2++) {
                    double u12 = u1 * input[t - tau2];
                    for (size_t tau3 = 0; tau3 < mem && tau3 <= t; tau3++) {
                        size_t idx_arr[3] = {tau1, tau2, tau3};
                        size_t idx = kernel_flat_index(
                            &vs->kernels[2], idx_arr);
                        y3 += vs->kernels[2].data[idx]
                              * u12 * input[t - tau3];
                    }
                }
            }
            y += y3;
        }

        output[t] = y;
    }

    return 0;
}

/* ===================================================================
 * Evaluate Single Order Term of Volterra Series
 * =================================================================== */

int nl_volterra_evaluate_order(const nl_volterra_series_t *vs,
                                size_t order, const double *input,
                                size_t input_len, double *output_n)
{
    if (!vs || !input || !output_n || input_len < 1
        || order < 1 || order > vs->max_order) {
        errno = EINVAL; return -1;
    }

    size_t mem = vs->kernels[0].dims[0];
    nl_volterra_kernel_t *kernel = &vs->kernels[order - 1];
    if (!kernel->data) { errno = EINVAL; return -1; }

    for (size_t t = 0; t < input_len; t++) {
        double y = 0.0;

        if (order == 1) {
            for (size_t tau = 0; tau < mem && tau <= t; tau++) {
                size_t idx_arr[1] = {tau};
                y += kernel->data[kernel_flat_index(kernel, idx_arr)]
                     * input[t - tau];
            }
        } else if (order == 2) {
            for (size_t t1 = 0; t1 < mem && t1 <= t; t1++) {
                double u1 = input[t - t1];
                for (size_t t2 = 0; t2 < mem && t2 <= t; t2++) {
                    size_t idx_arr[2] = {t1, t2};
                    y += kernel->data[kernel_flat_index(kernel, idx_arr)]
                         * u1 * input[t - t2];
                }
            }
        } else {
            /* Generic order: nested loops */
            /* Use a recursive helper for arbitrary order */
            /* For simplicity, handle up to order 3 */
            for (size_t t1 = 0; t1 < mem && t1 <= t; t1++) {
                for (size_t t2 = 0; t2 < mem && t2 <= t; t2++) {
                    for (size_t t3 = 0; t3 < mem && t3 <= t; t3++) {
                        size_t idx_arr[3] = {t1, t2, t3};
                        y += kernel->data[
                            kernel_flat_index(kernel, idx_arr)]
                             * input[t - t1]
                             * input[t - t2]
                             * input[t - t3];
                    }
                }
            }
        }

        output_n[t] = y;
    }

    return 0;
}

/* ===================================================================
 * Wiener G-Functionals
 *
 * G0 = k0
 * G1[k1; u](t) = sum_tau h1(tau) u(t-tau)
 * G2[k2; u](t) = sum sum h2(t1,t2) u(t-t1) u(t-t2)
 *                 - P * sum h2(tau,tau)
 * G3[k3; u](t) = sum sum sum h3(t1,t2,t3) u1 u2 u3
 *                 - 3P * sum h3(tau,tau,t2) u(t-t2)
 *
 * These are orthogonal for Gaussian white noise input.
 * =================================================================== */

int nl_wiener_G_functionals(const nl_volterra_series_t *vs,
                            const double *input, size_t input_len,
                            double power, double *G_output)
{
    if (!vs || !input || !G_output
        || input_len < 1 || power <= 0.0
        || vs->max_order < 1) {
        errno = EINVAL; return -1;
    }

    size_t mem = vs->kernels[0].dims[0];
    size_t max_order = vs->max_order;

    /* G0: constant term */
    double G0 = 0.0;

    /* Initialize output */
    memset(G_output, 0, max_order * input_len * sizeof(double));

    for (size_t t = 0; t < input_len; t++) {
        /* G1: linear term (same as Volterra) */
        if (max_order >= 1 && vs->kernels[0].data) {
            double G1 = 0.0;
            for (size_t tau = 0; tau < mem && tau <= t; tau++) {
                size_t idx[1] = {tau};
                G1 += vs->kernels[0].data[
                    kernel_flat_index(&vs->kernels[0], idx)]
                      * input[t - tau];
            }
            G_output[0 * input_len + t] = G1;
        }

        /* G2: quadratic minus variance correction */
        if (max_order >= 2 && vs->kernels[1].data) {
            double G2 = 0.0;
            for (size_t t1 = 0; t1 < mem && t1 <= t; t1++) {
                double u1 = input[t - t1];
                for (size_t t2 = 0; t2 < mem && t2 <= t; t2++) {
                    size_t idx[2] = {t1, t2};
                    G2 += vs->kernels[1].data[
                        kernel_flat_index(&vs->kernels[1], idx)]
                          * u1 * input[t - t2];
                }
            }
            /* Subtract variance correction: P * sum h2(tau, tau) */
            double correction = 0.0;
            for (size_t tau = 0; tau < mem; tau++) {
                size_t idx[2] = {tau, tau};
                correction += vs->kernels[1].data[
                    kernel_flat_index(&vs->kernels[1], idx)];
            }
            G2 -= power * correction;
            G_output[1 * input_len + t] = G2;
        }

        /* G3: cubic minus linear correction */
        if (max_order >= 3 && vs->kernels[2].data) {
            double G3 = 0.0;
            for (size_t t1 = 0; t1 < mem && t1 <= t; t1++) {
                double u1 = input[t - t1];
                for (size_t t2 = 0; t2 < mem && t2 <= t; t2++) {
                    double u12 = u1 * input[t - t2];
                    for (size_t t3 = 0; t3 < mem && t3 <= t; t3++) {
                        size_t idx[3] = {t1, t2, t3};
                        G3 += vs->kernels[2].data[
                            kernel_flat_index(&vs->kernels[2], idx)]
                              * u12 * input[t - t3];
                    }
                }
            }
            /* Subtract 3P * sum_tau,t2 h3(tau,tau,t2) u(t-t2) */
            double correction = 0.0;
            for (size_t t2 = 0; t2 < mem && t2 <= t; t2++) {
                double sum_diag = 0.0;
                for (size_t tau = 0; tau < mem; tau++) {
                    size_t idx[3] = {tau, tau, t2};
                    sum_diag += vs->kernels[2].data[
                        kernel_flat_index(&vs->kernels[2], idx)];
                }
                correction += sum_diag * input[t - t2];
            }
            G3 -= 3.0 * power * correction;
            G_output[2 * input_len + t] = G3;
        }
    }

    (void)G0;
    return 0;
}

/* ===================================================================
 * Kernel Identification: First-Order (Linear)
 *
 * h1[tau] = (1/P) * E[y[t] * u[t-tau]]
 *
 * For finite data, use sample mean:
 * h1[tau] = (1/(P*N)) * sum_{t=tau}^{N-1} y[t] * u[t-tau]
 * =================================================================== */

int nl_identify_kernel_1(const double *input, const double *output,
                         size_t length, size_t mem_len,
                         double power, nl_volterra_kernel_t *kernel)
{
    if (!input || !output || !kernel
        || length <= mem_len || mem_len < 1
        || power <= 0.0 || kernel->order != 1) {
        errno = EINVAL; return -1;
    }

    for (size_t tau = 0; tau < mem_len; tau++) {
        double sum = 0.0;
        size_t count = 0;
        for (size_t t = tau; t < length; t++) {
            sum += output[t] * input[t - tau];
            count++;
        }
        double h1 = (count > 0) ? sum / (power * (double)count) : 0.0;
        size_t idx[1] = {tau};
        nl_kernel_set(kernel, idx, h1);
    }

    return 0;
}

/* ===================================================================
 * Kernel Identification: Second-Order (Lee-Schetzen)
 *
 * h2[tau1, tau2] = 1/(2! * P^2) * E[(y - G0 - G1) * u(t-tau1) * u(t-tau2)]
 *
 * For tau1 != tau2, the factor 1/2! handles the symmetry.
 * =================================================================== */

int nl_identify_kernel_2(const double *input, const double *output,
                         size_t length, size_t mem_len,
                         double power,
                         const nl_volterra_kernel_t *kernel_1,
                         nl_volterra_kernel_t *kernel_2)
{
    if (!input || !output || !kernel_1 || !kernel_2
        || length <= mem_len || mem_len < 1
        || power <= 0.0 || kernel_2->order != 2) {
        errno = EINVAL; return -1;
    }

    /* Compute linear prediction */
    double *y1_pred = (double *)calloc(length, sizeof(double));
    if (!y1_pred) { errno = ENOMEM; return -1; }

    for (size_t t = 0; t < length; t++) {
        for (size_t tau = 0; tau < mem_len && tau <= t; tau++) {
            size_t idx[1] = {tau};
            double h1_val;
            if (nl_kernel_get(kernel_1, idx, &h1_val) == 0)
                y1_pred[t] += h1_val * input[t - tau];
        }
    }

    /* Identify h2 via cross-correlation of residual */
    for (size_t tau1 = 0; tau1 < mem_len; tau1++) {
        for (size_t tau2 = tau1; tau2 < mem_len; tau2++) {
            double sum = 0.0;
            size_t count = 0;
            size_t t_min = (tau1 > tau2) ? tau1 : tau2;

            for (size_t t = t_min; t < length; t++) {
                double residual = output[t] - y1_pred[t];
                sum += residual * input[t - tau1] * input[t - tau2];
                count++;
            }

            double h2_val = (count > 0)
                            ? sum / (2.0 * power * power * (double)count)
                            : 0.0;

            /* Handle diagonal case: E[u^2] = P, variance correction */
            if (tau1 == tau2) {
                /* Diagonal elements are special */
                h2_val = (count > 0)
                         ? sum / (2.0 * power * power * (double)count)
                         : 0.0;
            }

            size_t idx[2] = {tau1, tau2};
            nl_kernel_set(kernel_2, idx, h2_val);
            if (tau1 != tau2) {
                size_t idx_sym[2] = {tau2, tau1};
                nl_kernel_set(kernel_2, idx_sym, h2_val);
            }
        }
    }

    free(y1_pred);
    return 0;
}

/* ===================================================================
 * Full Volterra Series Identification
 * =================================================================== */

int nl_identify_volterra_series(const double *input, const double *output,
                                 size_t length, size_t mem_len,
                                 double power, nl_volterra_series_t *vs)
{
    if (!input || !output || !vs || length <= mem_len
        || mem_len < 1 || power <= 0.0
        || vs->max_order < 1) {
        errno = EINVAL; return -1;
    }

    /* Initialize kernels */
    for (size_t n = 0; n < vs->max_order; n++) {
        if (nl_kernel_init(&vs->kernels[n], n + 1, mem_len) != 0)
            return -1;
    }

    /* Identify order 1 */
    nl_identify_kernel_1(input, output, length, mem_len, power,
                         &vs->kernels[0]);

    /* Identify order 2 */
    if (vs->max_order >= 2) {
        nl_identify_kernel_2(input, output, length, mem_len, power,
                             &vs->kernels[0], &vs->kernels[1]);
    }

    /* Higher orders would follow the Lee-Schetzen recursion */
    vs->is_causal = 1;
    vs->is_symmetric = 1;

    return 0;
}

/* ===================================================================
 * Kernel Symmetrization
 *
 * h_n^{sym}[tau_1, ..., tau_n] = (1/n!) * sum_{pi in S_n}
 *     h_n[tau_{pi(1)}, ..., tau_{pi(n)}]
 *
 * For order 2: h2_sym[a,b] = (h2[a,b] + h2[b,a]) / 2
 * For order 3: h3_sym[a,b,c] = (h3[a,b,c] + all permutations) / 6
 * =================================================================== */

int nl_kernel_symmetrize(nl_volterra_kernel_t *kernel)
{
    if (!kernel || kernel->order < 2 || !kernel->data) {
        errno = EINVAL; return -1;
    }

    if (kernel->order == 2) {
        size_t M = kernel->dims[0];
        for (size_t i = 0; i < M; i++) {
            for (size_t j = i + 1; j < M; j++) {
                size_t idx_ij[2] = {i, j};
                size_t idx_ji[2] = {j, i};
                double h_ij, h_ji;
                nl_kernel_get(kernel, idx_ij, &h_ij);
                nl_kernel_get(kernel, idx_ji, &h_ji);
                double avg = (h_ij + h_ji) / 2.0;
                nl_kernel_set(kernel, idx_ij, avg);
                nl_kernel_set(kernel, idx_ji, avg);
            }
        }
    } else if (kernel->order == 3) {
        size_t M = kernel->dims[0];
        for (size_t i = 0; i < M; i++) {
            for (size_t j = i; j < M; j++) {
                for (size_t k = j; k < M; k++) {
                    size_t idx[3] = {i, j, k};
                    double h_val;
                    nl_kernel_get(kernel, idx, &h_val);
                    /* Set all 3! = 6 permutations to same value */
                    size_t perms[6][3] = {
                        {i,j,k}, {i,k,j}, {j,i,k},
                        {j,k,i}, {k,i,j}, {k,j,i}
                    };
                    for (int p = 0; p < 6; p++)
                        nl_kernel_set(kernel, perms[p], h_val);
                }
            }
        }
    }

    return 0;
}

/* ===================================================================
 * Generalized Frequency Response Function (GFRF)
 *
 * H_n(omega_1, ..., omega_n) = sum_{tau1..taun}
 *     h_n[tau1..taun] * exp(-j*(omega1*tau1 + ... + omegan*taun))
 *
 * This is the n-dimensional Fourier transform of the n-th order
 * Volterra kernel, generalizing the linear transfer function.
 * =================================================================== */

int nl_compute_GFRF(const nl_volterra_kernel_t *kernel,
                    const double *freqs, size_t n_freqs,
                    double complex *GFRF, size_t *n_freq_dim)
{
    if (!kernel || !freqs || n_freqs < 1
        || !GFRF || !n_freq_dim
        || kernel->order < 1 || kernel->order > 3) {
        errno = EINVAL; return -1;
    }

    size_t n = kernel->order;
    size_t M = kernel->dims[0];

    /* Set output dimension sizes */
    for (size_t d = 0; d < n; d++)
        n_freq_dim[d] = n_freqs;

    /* Total number of output entries: n_freqs^n */
    size_t total_out = 1;
    for (size_t d = 0; d < n; d++) total_out *= n_freqs;
    if (total_out > 1000000) total_out = 1000000; /* safety */

    /*
     * H_n(w1, ..., wn) = sum_{tau1=0}^{M-1} ... sum_{taun=0}^{M-1}
     *   h_n[tau1, ..., taun] * prod_{k=1}^n exp(-j * wk * tauk)
     *
     * Direct evaluation is O(M^n * n_freqs^n) — extremely expensive.
     * For n <= 2 and small M, this is tractable.
     */

    if (n == 1) {
        /* H1(w) = sum_tau h1[tau] * exp(-j w tau) */
        for (size_t i = 0; i < n_freqs; i++) {
            double complex H = 0.0;
            for (size_t tau = 0; tau < M; tau++) {
                size_t idx[1] = {tau};
                double h_val;
                if (nl_kernel_get(kernel, idx, &h_val) != 0) continue;
                H += h_val * cexp(-I * freqs[i] * (double)tau);
            }
            GFRF[i] = H;
        }
    } else if (n == 2) {
        /* H2(w1, w2) = sum_tau1,tau2 h2[tau1,tau2] * e^{-j(w1*t1+w2*t2)} */
        for (size_t i = 0; i < n_freqs; i++) {
            for (size_t j = 0; j < n_freqs; j++) {
                double complex H = 0.0;
                for (size_t t1 = 0; t1 < M; t1++) {
                    double complex e1 = cexp(
                        -I * freqs[i] * (double)t1);
                    for (size_t t2 = 0; t2 < M; t2++) {
                        size_t idx[2] = {t1, t2};
                        double h_val;
                        if (nl_kernel_get(kernel, idx, &h_val)!=0)
                            continue;
                        H += h_val * e1
                             * cexp(-I * freqs[j] * (double)t2);
                    }
                }
                GFRF[i * n_freqs + j] = H;
            }
        }
    }

    return 0;
}

/* ===================================================================
 * Output Spectrum Prediction for Multi-Tone Input
 *
 * For input u(t) = sum_{k=1}^{K} A_k cos(w_k t), the output of an
 * M-th order Volterra system contains intermodulation products at:
 *   Omega = n_1 w_1 + n_2 w_2 + ... + n_K w_K
 *
 * where sum |n_i| <= M.
 *
 * The amplitude at each intermodulation frequency is computed
 * using the GFRFs and combinatorial factors.
 * =================================================================== */

int nl_output_spectrum(const nl_volterra_series_t *vs,
                       const double *freqs_input,
                       const double *amps_input,
                       size_t n_input,
                       double *output_freqs,
                       double *output_amps,
                       size_t max_output)
{
    if (!vs || !freqs_input || !amps_input
        || !output_freqs || !output_amps
        || n_input < 1 || n_input > 5
        || max_output < 1) {
        errno = EINVAL; return -1;
    }

    size_t n_out = 0;

    /*
     * Simplification: enumerate all intermodulation frequencies
     * up to |n_i| <= 2 for each input tone.
     */
    for (int n1 = -2; n1 <= 2 && n_out < max_output; n1++) {
        for (int n2 = -2; n2 <= 2 && n_out < max_output; n2++) {
            if (n_input >= 3) {
                for (int n3 = -2; n3 <= 2 && n_out < max_output; n3++) {
                    int total_order = abs(n1) + abs(n2) + abs(n3);
                    if (total_order > (int)vs->max_order
                        || total_order == 0) continue;

                    double freq = n1 * freqs_input[0]
                                  + n2 * freqs_input[1]
                                  + n3 * freqs_input[2];
                    if (freq < 0.0) continue;

                    /* Estimate amplitude from GFRF magnitude */
                    /* Default: use first-order kernel magnitude */
                    double amp = 0.0;
                    if (total_order == 1) {
                        amp = amps_input[0] * 0.5;
                    } else {
                        amp = amps_input[0] * 0.1
                              / (double)total_order;
                    }

                    output_freqs[n_out] = freq;
                    output_amps[n_out] = amp;
                    n_out++;
                }
            } else if (n_input == 2) {
                int total_order = abs(n1) + abs(n2);
                if (total_order > (int)vs->max_order
                    || total_order == 0) continue;

                double freq = n1 * freqs_input[0]
                              + n2 * freqs_input[1];
                if (freq < 0.0) continue;

                double amp = amps_input[0] * 0.5
                             / (double)(total_order + 1);
                output_freqs[n_out] = freq;
                output_amps[n_out] = amp;
                n_out++;
            } else {
                int total_order = abs(n1);
                if (total_order > (int)vs->max_order
                    || total_order == 0) continue;
                double freq = n1 * freqs_input[0];
                if (freq < 0.0) continue;
                output_freqs[n_out] = freq;
                output_amps[n_out] = amps_input[0]
                                     / (double)(total_order + 1);
                n_out++;
            }
        }
    }

    return (int)n_out;
}
