/**
 * @file stability.h
 * @brief L4 Fundamental Laws - System Stability Analysis
 *
 * Stability is the most critical property of any dynamic system.
 *
 * Knowledge Coverage:
 *   L4-S01: BIBO stability test (CT: poles in LHP; DT: poles in unit circle)
 *   L4-S02: Routh-Hurwitz criterion (polynomial root location without solving)
 *   L4-S03: Jury stability test (DT equivalent of Routh-Hurwitz)
 *   L4-S04: Nyquist stability criterion (encirclements of -1 point)
 *   L4-S05: Gain margin and phase margin from Bode/Nyquist
 *   L4-S06: Root locus analysis (Evans root locus)
 *   L4-S07: Lyapunov stability (direct method)
 *   L4-S08: Marginal stability detection (poles on jw-axis / unit circle)
 *   L4-S09: Relative stability measures
 *   L4-S10: Stability of feedback interconnection
 *
 * Course Mapping:
 *   MIT 6.003 Ch.11; Stanford EE102A Ch.9;
 *   Berkeley EE123 Ch.7; Tsinghua Signals Ch.6
 */

#ifndef STABILITY_H
#define STABILITY_H

#include "system_defs.h"
#include "state_space.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- BIBO Stability ---- */

/** Check BIBO stability for CT system: all poles must have Re(s) < 0.
 *  Returns 1 if BIBO stable, 0 otherwise. */
int ct_is_bibo_stable(const ct_transfer_function_t *tf);

/** Check BIBO stability for DT system: all poles must have |z| < 1.
 *  Returns 1 if BIBO stable, 0 otherwise. */
int dt_is_bibo_stable(const dt_transfer_function_t *tf);

/** BIBO stability from pole list. */
int ct_is_bibo_stable_from_poles(const pole_zero_collection_t *pz);

int dt_is_bibo_stable_from_poles(const pole_zero_collection_t *pz);

/* ---- Routh-Hurwitz Criterion ---- */

/** Routh-Hurwitz table for polynomial a_n s^n + ... + a_0.
 *  Returns number of RHP poles (0 = stable, >0 = unstable).
 *  Also reports any special cases (entire row of zeros, first element zero).
 *  Complexity: O(n^2).
 *
 *  @param coeffs   Coefficients a_0 ... a_n (length n+1).
 *  @param order    Polynomial degree n.
 *  @param rh_table Output Routh array (caller-allocated, (n+1) x ceil((n+1)/2)).
 *  @param special_case  Output: 0=none, 1=zero first element, 2=entire row zero.
 *  @return Number of sign changes in first column (= RHP poles). */
int routh_hurwitz(const double *coeffs, int order,
                  double *rh_table, int *special_case);

/** Convenience: apply Routh-Hurwitz to CT transfer function denominator. */
int ct_routh_hurwitz(const ct_transfer_function_t *tf,
                     int *rhp_poles, int *special_case);

/* ---- Jury Stability Test (DT) ---- */

/** Jury stability test for DT polynomial a_n z^n + ... + a_0.
 *  All roots inside unit circle iff Jury conditions are met.
 *  Complexity: O(n^2).
 *
 *  @param coeffs  Coefficients a_0 ... a_n (length n+1).
 *  @param order   Polynomial degree n.
 *  @return 1 if stable (all roots |z|<1), 0 otherwise. */
int jury_stability_test(const double *coeffs, int order);

/** Convenience: apply Jury test to DT transfer function denominator. */
int dt_jury_stability_test(const dt_transfer_function_t *tf);

/* ---- Nyquist Stability Criterion ---- */

/** Nyquist stability analysis for open-loop transfer function L(s).
 *  N = Z - P, where N = encirclements of -1, Z = closed-loop RHP poles,
 *  P = open-loop RHP poles.
 *
 *  @param loop_tf   Open-loop transfer function L(s) = G(s)H(s).
 *  @param fr        Frequency response of L(jw) (pre-computed).
 *  @param open_loop_rhp_poles  Number of open-loop RHP poles P.
 *  @param encirclements  Output: net CCW encirclements of -1.
 *  @param closed_loop_rhp_poles Output: Z = P + N.
 *  @return 0 on success, -1 if insufficient frequency data. */
int nyquist_stability(const ct_transfer_function_t *loop_tf,
                      const frequency_response_t *fr,
                      int open_loop_rhp_poles,
                      int *encirclements,
                      int *closed_loop_rhp_poles);

/* ---- Gain and Phase Margin ---- */

/** Compute gain margin (dB) and phase margin (degrees) from Bode data.
 *  Gain margin: gain at -180 deg phase crossover.
 *  Phase margin: 180 + phase at 0 dB gain crossover.
 *
 *  @param fr             Frequency response data.
 *  @param gain_margin_db Output gain margin in dB.
 *  @param phase_margin_deg Output phase margin in degrees.
 *  @param gm_freq        Frequency at which gain margin is measured.
 *  @param pm_freq        Frequency at which phase margin is measured.
 *  @return 0 on success, -1 if margins cannot be determined. */
int compute_stability_margins(const frequency_response_t *fr,
                               double *gain_margin_db,
                               double *phase_margin_deg,
                               double *gm_freq,
                               double *pm_freq);

/* ---- Root Locus ---- */

/** Compute root locus: closed-loop pole trajectories as gain K varies.
 *  H_cl(s) = K*N(s) / (D(s) + K*N(s)).
 *
 *  @param tf       Open-loop transfer function.
 *  @param K_values Array of gain values.
 *  @param num_K    Number of gain values.
 *  @param poles    Output: (num_K x den_order) complex pole locations.
 *  @return 0 on success. */
int root_locus(const ct_transfer_function_t *tf,
               const double *K_values, int num_K,
               double complex *poles);

/* ---- Lyapunov Direct Method ---- */

/** Check CT stability via Lyapunov equation: A'*P + P*A = -Q.
 *  System is stable iff P > 0 (positive definite) for any Q > 0.
 *
 *  @param A      System matrix n x n (row-major).
 *  @param n      State dimension.
 *  @param Q      Positive definite matrix n x n (typically I).
 *  @param P      Output solution n x n.
 *  @return 1 if P is positive definite (system stable), 0 otherwise. */
int ct_lyapunov_stability(const double *A, int n,
                           const double *Q, double *P);

/* ---- Marginal Stability ---- */

/** Check for marginal stability: poles on jw-axis (CT) or unit circle (DT).
 *  Returns number of poles exactly on the stability boundary. */
int ct_count_marginal_poles(const pole_zero_collection_t *pz);

int dt_count_marginal_poles(const pole_zero_collection_t *pz);

/* ---- Feedback Stability ---- */

/** Check stability of negative feedback interconnection:
 *  H_cl = G / (1 + G*H). Stable iff 1+G(s)H(s) has no RHP zeros.
 *
 *  @param G  Forward path transfer function.
 *  @param H  Feedback path transfer function.
 *  @return 1 if closed-loop stable, 0 otherwise. */
int feedback_stability_check(const ct_transfer_function_t *G,
                              const ct_transfer_function_t *H);

/* ---- Relative Stability ---- */

/** Compute damping ratio and natural frequency for dominant pole pair. */
int ct_dominant_pole_analysis(const pole_zero_collection_t *pz,
                               double *damping_ratio,
                               double *natural_freq);

/** Settling time estimate: t_s ≈ 4/(ζ·ω_n) for 2% criterion. */
double ct_settling_time_estimate(double damping_ratio, double natural_freq);

/** Peak overshoot estimate: M_p = exp(-π·ζ/√(1-ζ²)) for 0 ≤ ζ < 1. */
double ct_overshoot_estimate(double damping_ratio);

#ifdef __cplusplus
}
#endif

#endif /* STABILITY_H */
