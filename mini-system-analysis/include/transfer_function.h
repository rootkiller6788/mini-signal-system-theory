/**
 * @file transfer_function.h
 * @brief L3 Math Structures - Transfer Functions, Pole-Zero, Partial Fractions
 *
 * H(s) = N(s)/D(s) in CT Laplace domain, H(z) = N(z)/D(z) in DT Z domain.
 *
 * Knowledge Coverage:
 *   L3-M01: H(s) evaluation at complex s
 *   L3-M02: H(z) evaluation at complex z
 *   L3-M03: Pole computation via companion matrix QR
 *   L3-M04: Zero computation
 *   L3-M05: Partial fraction expansion (residues)
 *   L3-M06: Bilinear transform CT->DT
 *   L3-M07: Impulse invariance CT->DT
 *   L3-M08: ZOH discretization
 *   L3-M09: Matched Z-transform
 *   L3-M10: ZPK (Zero-Pole-Gain) decomposition
 *   L3-M11: Bode plot generation
 *   L3-M12: Nyquist plot generation
 *   L3-M13: Group delay computation
 *   L3-M14: System order analysis
 *
 * Course Mapping:
 *   MIT 6.003 Ch.9-10; Stanford EE102A Ch.7-8;
 *   Berkeley EE123 Ch.6; Tsinghua Signals Ch.4-5
 */

#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include "system_defs.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Transfer Function Evaluation ---- */

double complex ct_tf_evaluate(const ct_transfer_function_t *tf,
                               double complex s);

double complex dt_tf_evaluate(const dt_transfer_function_t *tf,
                               double complex z);

int ct_tf_evaluate_vector(const ct_transfer_function_t *tf,
                           const double complex *s_values, int num_points,
                           double complex *result);

int dt_tf_evaluate_vector(const dt_transfer_function_t *tf,
                           const double complex *z_values, int num_points,
                           double complex *result);

/* ---- Frequency Response ---- */

int ct_frequency_response(const ct_transfer_function_t *tf,
                          frequency_response_t *fr);

int dt_frequency_response(const dt_transfer_function_t *tf,
                          double sample_rate,
                          frequency_response_t *fr);

/* ---- Pole-Zero Analysis ---- */

int ct_tf_find_poles(const ct_transfer_function_t *tf,
                     pole_zero_collection_t *pz);

int ct_tf_find_zeros(const ct_transfer_function_t *tf,
                     pole_zero_collection_t *pz);

int dt_tf_find_poles(const dt_transfer_function_t *tf,
                     pole_zero_collection_t *pz);

int dt_tf_find_zeros(const dt_transfer_function_t *tf,
                     pole_zero_collection_t *pz);

int ct_tf_pole_zero_analysis(const ct_transfer_function_t *tf,
                              pole_zero_collection_t *pz);

int dt_tf_pole_zero_analysis(const dt_transfer_function_t *tf,
                              pole_zero_collection_t *pz);

/* ---- Partial Fraction Expansion ---- */

int ct_tf_partial_fraction(const ct_transfer_function_t *tf,
                           double complex *residues,
                           double complex *poles,
                           double *direct_term, int *direct_len);

int dt_tf_partial_fraction(const dt_transfer_function_t *tf,
                           double complex *residues,
                           double complex *poles,
                           double *direct_term, int *direct_len);

/* ---- CT to DT Transformations ---- */

int ct_to_dt_bilinear(const ct_transfer_function_t *ct_tf,
                      double sample_time,
                      dt_transfer_function_t *dt_tf);

int ct_to_dt_impulse_invariance(const ct_transfer_function_t *ct_tf,
                                 double sample_time,
                                 dt_transfer_function_t *dt_tf);

int ct_to_dt_zoh(const ct_transfer_function_t *ct_tf,
                 double sample_time,
                 dt_transfer_function_t *dt_tf);

int ct_to_dt_matched_z(const ct_transfer_function_t *ct_tf,
                        double sample_time,
                        dt_transfer_function_t *dt_tf);

/* ---- System Order Analysis ---- */

int ct_tf_order_analysis(const ct_transfer_function_t *ct_tf,
                          system_order_info_t *info);

int dt_tf_order_analysis(const dt_transfer_function_t *dt_tf,
                          system_order_info_t *info);

/* ---- ZPK Decomposition ---- */

double ct_tf_zpk_gain(const ct_transfer_function_t *tf);

double dt_tf_zpk_gain(const dt_transfer_function_t *tf);

/* ---- Bode / Nyquist / Group Delay ---- */

int ct_tf_bode_plot(const ct_transfer_function_t *tf,
                    frequency_response_t *fr);

int ct_tf_nyquist_plot(const ct_transfer_function_t *tf,
                       frequency_response_t *fr);

int ct_tf_group_delay(const ct_transfer_function_t *tf,
                      const double *frequencies, int num_freq,
                      double *group_delay);

int dt_tf_group_delay(const dt_transfer_function_t *tf,
                      const double *frequencies, int num_freq,
                      double *group_delay);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFER_FUNCTION_H */
