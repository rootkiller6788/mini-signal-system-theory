#include "signal_ops.h"
#include "signal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =================================================================
 * L2: Time-Domain Operations
 * ================================================================= */

int signal_ct_time_shift(const signal_ct_t *src, signal_time_t t_shift, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        double t_dst = dst->t_start + (double)k * dst->dt;
        double t_src = t_dst - t_shift;
        if (t_src < src->t_start || t_src > src->t_end)
            dst->samples[k] = 0.0;
        else {
            double idx = (t_src - src->t_start) / src->dt;
            signal_index_t i0 = (signal_index_t)idx;
            if (i0 >= src->num_samples - 1)
                dst->samples[k] = src->samples[src->num_samples - 1];
            else {
                double frac = idx - (double)i0;
                dst->samples[k] = src->samples[i0] * (1.0 - frac) + src->samples[i0 + 1] * frac;
            }
        }
    }
    return 0;
}

int signal_ct_time_reverse(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        double t_dst = dst->t_start + (double)k * dst->dt;
        double t_src = -t_dst;
        if (t_src < src->t_start || t_src > src->t_end)
            dst->samples[k] = 0.0;
        else {
            double idx = (t_src - src->t_start) / src->dt;
            signal_index_t i0 = (signal_index_t)idx;
            if (i0 >= src->num_samples - 1)
                dst->samples[k] = src->samples[src->num_samples - 1];
            else {
                double frac = idx - (double)i0;
                dst->samples[k] = src->samples[i0] * (1.0 - frac) + src->samples[i0 + 1] * frac;
            }
        }
    }
    return 0;
}

int signal_ct_time_scale(const signal_ct_t *src, double a, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples || fabs(a) < 1e-15) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        double t_dst = dst->t_start + (double)k * dst->dt;
        double t_src = t_dst * a;
        if (t_src < src->t_start || t_src > src->t_end)
            dst->samples[k] = 0.0;
        else {
            double idx = (t_src - src->t_start) / src->dt;
            signal_index_t i0 = (signal_index_t)idx;
            if (i0 >= src->num_samples - 1)
                dst->samples[k] = src->samples[src->num_samples - 1];
            else {
                double frac = idx - (double)i0;
                dst->samples[k] = src->samples[i0] * (1.0 - frac) + src->samples[i0 + 1] * frac;
            }
        }
    }
    return 0;
}

int signal_dt_circular_shift(const signal_dt_t *src, signal_index_t shift, signal_dt_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t N = src->length;
    if (dst->length != N) return -1;
    signal_index_t s = ((shift % N) + N) % N;
    for (signal_index_t i = 0; i < N; i++)
        dst->samples[i] = src->samples[(i - s + N) % N];
    return 0;
}

/* =================================================================
 * L2: Amplitude Operations
 * ================================================================= */

int signal_ct_amplitude_scale(const signal_ct_t *src, double gain, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    for (signal_index_t k = 0; k < n; k++)
        dst->samples[k] = src->samples[k] * gain;
    return 0;
}

int signal_ct_add_offset(const signal_ct_t *src, double offset, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    for (signal_index_t k = 0; k < n; k++)
        dst->samples[k] = src->samples[k] + offset;
    return 0;
}

/* =================================================================
 * L2: Signal Arithmetic (Vector Space Operations)
 * ================================================================= */

int signal_ct_add(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result)
{
    if (!a || !b || !result || !a->samples || !b->samples || !result->samples) return -1;
    signal_index_t n = a->num_samples;
    if (b->num_samples < n) n = b->num_samples;
    if (result->num_samples < n) n = result->num_samples;
    for (signal_index_t k = 0; k < n; k++)
        result->samples[k] = a->samples[k] + b->samples[k];
    return 0;
}

int signal_ct_subtract(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result)
{
    if (!a || !b || !result || !a->samples || !b->samples || !result->samples) return -1;
    signal_index_t n = a->num_samples;
    if (b->num_samples < n) n = b->num_samples;
    if (result->num_samples < n) n = result->num_samples;
    for (signal_index_t k = 0; k < n; k++)
        result->samples[k] = a->samples[k] - b->samples[k];
    return 0;
}

int signal_ct_multiply(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result)
{
    if (!a || !b || !result || !a->samples || !b->samples || !result->samples) return -1;
    signal_index_t n = a->num_samples;
    if (b->num_samples < n) n = b->num_samples;
    if (result->num_samples < n) n = result->num_samples;
    for (signal_index_t k = 0; k < n; k++)
        result->samples[k] = a->samples[k] * b->samples[k];
    return 0;
}

int signal_ct_divide(const signal_ct_t *a, const signal_ct_t *b, double epsilon, signal_ct_t *result)
{
    if (!a || !b || !result || !a->samples || !b->samples || !result->samples) return -1;
    signal_index_t n = a->num_samples;
    if (b->num_samples < n) n = b->num_samples;
    if (result->num_samples < n) n = result->num_samples;
    for (signal_index_t k = 0; k < n; k++) {
        if (fabs(b->samples[k]) < epsilon)
            result->samples[k] = 0.0;
        else
            result->samples[k] = a->samples[k] / b->samples[k];
    }
    return 0;
}

/* =================================================================
 * L2: Even/Odd Decomposition
 * x_even(t) = (x(t) + x(-t))/2, x_odd(t) = (x(t) - x(-t))/2
 * ================================================================= */

int signal_ct_even_odd_decompose(const signal_ct_t *src, signal_ct_t *even_part, signal_ct_t *odd_part)
{
    if (!src || !even_part || !odd_part) return -1;
    if (!src->samples || !even_part->samples || !odd_part->samples) return -1;
    signal_index_t N = src->num_samples;
    double half = 0.5;
    for (signal_index_t k = 0; k < N; k++) {
        double t = src->t_start + (double)k * src->dt;
        double v_plus = src->samples[k];
        double t_neg = -t;
        double v_minus;
        if (t_neg < src->t_start || t_neg > src->t_end)
            v_minus = 0.0;
        else {
            double idx = (t_neg - src->t_start) / src->dt;
            signal_index_t i0 = (signal_index_t)idx;
            if (i0 >= N - 1) v_minus = src->samples[N - 1];
            else {
                double frac = idx - (double)i0;
                v_minus = src->samples[i0] * (1.0 - frac) + src->samples[i0 + 1] * frac;
            }
        }
        even_part->samples[k] = (v_plus + v_minus) * half;
        odd_part->samples[k] = (v_plus - v_minus) * half;
    }
    return 0;
}

/* =================================================================
 * L3: Convolution
 * y(t) = (x * h)(t) = integral_{-inf}^{inf} x(tau) * h(t - tau) dtau
 * Reference: Oppenheim & Willsky (1997) section 2.1
 * ================================================================= */

int signal_ct_convolve(const signal_ct_t *x, const signal_ct_t *h, signal_ct_t *y)
{
    if (!x || !h || !y || !x->samples || !h->samples || !y->samples) return -1;
    signal_index_t Ny = y->num_samples;
    for (signal_index_t ky = 0; ky < Ny; ky++) {
        double ty = y->t_start + (double)ky * y->dt;
        double sum = 0.0;
        for (signal_index_t kh = 0; kh < h->num_samples; kh++) {
            double th = h->t_start + (double)kh * h->dt;
            double tx = ty - th;
            if (tx >= x->t_start && tx <= x->t_end) {
                double idx = (tx - x->t_start) / x->dt;
                signal_index_t ix = (signal_index_t)idx;
                if (ix < x->num_samples - 1) {
                    double frac = idx - (double)ix;
                    double xval = x->samples[ix] * (1.0 - frac) + x->samples[ix + 1] * frac;
                    sum += xval * h->samples[kh];
                }
            }
        }
        y->samples[ky] = sum * h->dt;
    }
    return 0;
}

int signal_dt_circular_convolve(const signal_dt_t *x, const signal_dt_t *h, signal_dt_t *y)
{
    if (!x || !h || !y || !x->samples || !h->samples || !y->samples) return -1;
    signal_index_t N = x->length;
    if (h->length != N || y->length != N) return -1;
    for (signal_index_t n = 0; n < N; n++) {
        double sum = 0.0;
        for (signal_index_t m = 0; m < N; m++)
            sum += x->samples[m] * h->samples[(n - m + N) % N];
        y->samples[n] = sum;
    }
    return 0;
}

int signal_dt_linear_convolve(const signal_dt_t *x, const signal_dt_t *h, signal_dt_t *y)
{
    if (!x || !h || !y || !x->samples || !h->samples || !y->samples) return -1;
    signal_index_t Nx = x->length, Nh = h->length, Ny = Nx + Nh - 1;
    if (y->length != Ny) return -1;
    for (signal_index_t n = 0; n < Ny; n++) {
        double sum = 0.0;
        for (signal_index_t k = 0; k < Nx; k++) {
            signal_index_t idx = n - k;
            if (idx >= 0 && idx < Nh)
                sum += x->samples[k] * h->samples[idx];
        }
        y->samples[n] = sum;
    }
    return 0;
}

/* =================================================================
 * L5: Signal Upsampling (Zero-Insertion + optional interpolation)
 * Insert L-1 zeros between each sample, expanding bandwidth by factor L.
 * ================================================================= */

int signal_ct_upsample(const signal_ct_t *src, int L, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples || L <= 1) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        if (k % L == 0) {
            signal_index_t src_idx = k / L;
            dst->samples[k] = (src_idx < src->num_samples) ? src->samples[src_idx] : 0.0;
        } else {
            dst->samples[k] = 0.0;
        }
    }
    return 0;
}

/* =================================================================
 * L5: Moving Average Filter (Smoothing)
 * y[k] = (1/M) * sum_{j=0}^{M-1} x[k-j]  (causal)
 * Low-pass FIR filter with all coefficients = 1/M.
 * Reference: Proakis & Manolakis (2007) section 10.2.1
 * ================================================================= */

int signal_ct_moving_average(const signal_ct_t *src, int window_size, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples || window_size < 1) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    for (signal_index_t k = 0; k < n; k++) {
        double sum = 0.0;
        int count = 0;
        for (int j = 0; j < window_size; j++) {
            signal_index_t idx = k - j;
            if (idx >= 0 && idx < n) {
                sum += src->samples[idx];
                count++;
            }
        }
        dst->samples[k] = (count > 0) ? sum / (double)count : 0.0;
    }
    return 0;
}

/* =================================================================
 * L5: First-Order Difference (Simple HPF)
 * y[k] = x[k] - x[k-1], approximates derivative.
 * ================================================================= */

int signal_ct_first_difference(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    dst->samples[0] = 0.0;
    for (signal_index_t k = 1; k < n; k++)
        dst->samples[k] = src->samples[k] - src->samples[k-1];
    return 0;
}

/* =================================================================
 * L5: Cumulative Sum (Discrete-Time Integrator)
 * y[k] = sum_{j=0}^{k} x[j], impulse response = unit step.
 * ================================================================= */

int signal_ct_cumulative_sum(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    if (n == 0) return -1;
    dst->samples[0] = src->samples[0];
    for (signal_index_t k = 1; k < n; k++)
        dst->samples[k] = dst->samples[k-1] + src->samples[k];
    return 0;
}
