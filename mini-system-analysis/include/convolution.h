/**
 * @file convolution.h
 * @brief L2 Core Concepts - Convolution Operations for LTI System Analysis
 *
 * CT: y(t) = x(t) * h(t) = integral x(tau) h(t-tau) dtau
 * DT: y[n] = x[n] * h[n] = sum x[k] h[n-k]
 *
 * Knowledge Coverage:
 *   L2-C01: CT convolution (trapezoidal integration)
 *   L2-C02: DT linear convolution (direct O(N^2))
 *   L2-C03: DT circular convolution
 *   L2-C04: Overlap-add fast convolution
 *   L2-C05: Overlap-save fast convolution
 *   L2-C06: Deconvolution (frequency-domain)
 *   L2-C07: Cross-correlation and auto-correlation
 *   L2-C08: Convolution properties verification
 *   L2-C09: Impulse/step response relationship
 *   L2-C10: Normalized cross-correlation
 *
 * Course Mapping:
 *   MIT 6.003 Ch.2; Stanford EE102A Ch.4;
 *   Berkeley EE123 Ch.3; Tsinghua Signals Ch.2
 */

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "system_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- CT Convolution ---- */

int ct_convolution_direct(const ct_signal_t *x, const ct_signal_t *h,
                          ct_signal_t *result);

int ct_convolution_fast(const ct_signal_t *x, const ct_signal_t *h,
                        ct_signal_t *result, int decimation_factor);

/* ---- DT Linear Convolution ---- */

int dt_convolution_direct(const dt_signal_t *x, const dt_signal_t *h,
                          dt_signal_t *y);

int dt_convolution_fft(const dt_signal_t *x, const dt_signal_t *h,
                       dt_signal_t *y);

/* ---- Circular Convolution ---- */

int dt_circular_convolution_direct(const dt_signal_t *x, const dt_signal_t *h,
                                    dt_signal_t *y);

int dt_circular_convolution_fft(const dt_signal_t *x, const dt_signal_t *h,
                                 dt_signal_t *y);

/* ---- Overlap-Add / Overlap-Save ---- */

int dt_overlap_add(const dt_signal_t *x, const dt_signal_t *h,
                   int block_size, dt_signal_t *y);

int dt_overlap_save(const dt_signal_t *x, const dt_signal_t *h,
                    int block_size, dt_signal_t *y);

/* ---- Deconvolution ---- */

int dt_deconvolution_freq(const dt_signal_t *y, const dt_signal_t *h,
                          double epsilon, dt_signal_t *x_est);

/* ---- Correlation ---- */

int dt_cross_correlation(const dt_signal_t *x, const dt_signal_t *y,
                         dt_signal_t *r_xy);

int dt_auto_correlation(const dt_signal_t *x, dt_signal_t *r_xx);

int dt_cross_correlation_normalized(const dt_signal_t *x, const dt_signal_t *y,
                                     dt_signal_t *rho_xy);

/* ---- Impulse/Step Relationship ---- */

int dt_step_from_impulse(const dt_signal_t *h, dt_signal_t *s);

int dt_impulse_from_step(const dt_signal_t *s, dt_signal_t *h);

/* ---- Utility ---- */

double dt_convolution_commutativity_check(const dt_signal_t *x,
                                           const dt_signal_t *h);

int dt_signal_zero_pad(const dt_signal_t *src, size_t new_len,
                       dt_signal_t *dst);

#ifdef __cplusplus
}
#endif

#endif /* CONVOLUTION_H */
