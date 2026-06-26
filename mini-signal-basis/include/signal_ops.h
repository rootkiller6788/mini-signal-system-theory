/*
 * signal_ops.h - Signal Operations (L2: Core Concepts / L3: Math Structures)
 *
 * Time-domain operations: shift, scale, reverse
 * Amplitude operations: gain, offset
 * Signal arithmetic: add, subtract, multiply, divide
 * Even/odd decomposition
 * Convolution: linear, circular, via DFT
 *
 * Course: MIT 6.003 Ch.1-2 / Stanford EE102A Ch.1 / Berkeley EE16B
 * Textbook: Oppenheim & Willsky (1997) section 1.2, 2.1
 */
#ifndef SIGNAL_OPS_H
#define SIGNAL_OPS_H
#include "signal_basis.h"

/* L2: Time-Domain Operations */
int signal_ct_time_shift(const signal_ct_t *src, signal_time_t t_shift, signal_ct_t *dst);
int signal_ct_time_reverse(const signal_ct_t *src, signal_ct_t *dst);
int signal_ct_time_scale(const signal_ct_t *src, double a, signal_ct_t *dst);
int signal_dt_circular_shift(const signal_dt_t *src, signal_index_t shift, signal_dt_t *dst);

/* L2: Amplitude Operations */
int signal_ct_amplitude_scale(const signal_ct_t *src, double gain, signal_ct_t *dst);
int signal_ct_add_offset(const signal_ct_t *src, double offset, signal_ct_t *dst);

/* L2: Signal Arithmetic (Vector Space Operations) */
int signal_ct_add(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result);
int signal_ct_subtract(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result);
int signal_ct_multiply(const signal_ct_t *a, const signal_ct_t *b, signal_ct_t *result);
int signal_ct_divide(const signal_ct_t *a, const signal_ct_t *b, double epsilon, signal_ct_t *result);

/* L2: Even/Odd Decomposition - Every signal x(t) = x_even(t) + x_odd(t)
 * x_even(t) = (x(t) + x(-t))/2, x_odd(t) = (x(t) - x(-t))/2
 * Reference: O&W section 1.2.3 */
int signal_ct_even_odd_decompose(const signal_ct_t *src, signal_ct_t *even_part, signal_ct_t *odd_part);

/* L3: Convolution - y(t) = (x * h)(t) = integral x(tau)*h(t-tau) dtau
 * Reference: O&W section 2.1 */
int signal_ct_convolve(const signal_ct_t *x, const signal_ct_t *h, signal_ct_t *y);

/* L3: Circular convolution: y[n] = sum x[k]*h[(n-k) mod N]
 * Computed via DFT: y = IFFT(FFT(x) * FFT(h)) */
int signal_dt_circular_convolve(const signal_dt_t *x, const signal_dt_t *h, signal_dt_t *y);

/* L3: Linear convolution via DFT (zero-padding) */
int signal_dt_linear_convolve(const signal_dt_t *x, const signal_dt_t *h, signal_dt_t *y);


/* L5: Additional signal processing operations */
int signal_ct_upsample(const signal_ct_t *src, int L, signal_ct_t *dst);
int signal_ct_moving_average(const signal_ct_t *src, int window_size, signal_ct_t *dst);
int signal_ct_first_difference(const signal_ct_t *src, signal_ct_t *dst);
int signal_ct_cumulative_sum(const signal_ct_t *src, signal_ct_t *dst);

#endif /* SIGNAL_OPS_H */
