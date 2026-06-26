/**
 * @file nonlinear_describing.h
 * @brief Describing function (DF) analysis for nonlinear systems.
 *
 * The describing function method is a quasi-linearization technique
 * that approximates a nonlinear element by its fundamental harmonic
 * response when driven by a sinusoidal input.
 *
 * L2: Core Concepts — describing function approximation, limit cycle prediction
 * L5: Algorithms — analytical DF formulas, numerical DF via FFT,
 *     harmonic balance, dual-input describing functions (DIDF)
 * L6: Canonical Problems — limit cycle existence and stability,
 *     oscillation amplitude/frequency prediction
 *
 * Course Mapping:
 *   - MIT 6.302: Feedback Systems (DF chapter)
 *   - Stanford EE 364A: Convex Optimization (LMIs for DF)
 *   - Tsinghua: 自动控制原理 (describing function method)
 *
 * Reference:
 *   - Gelb & Vander Velde, "Multiple-Input Describing Functions and
 *     Nonlinear System Design", McGraw-Hill, 1968
 *   - Atherton, "Nonlinear Control Engineering", Van Nostrand, 1975
 */

#ifndef NONLINEAR_DESCRIBING_H
#define NONLINEAR_DESCRIBING_H

#include "nonlinear_core.h"
#include "nonlinear_stability.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L2, L5: Describing Function Computation
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute the describing function N(A) for a memoryless nonlinearity.
 *
 * If the nonlinearity is odd and memoryless, the describing function is:
 *
 *   N(A) = (2/(πA)) ∫₀^π f(A sin θ) sin θ dθ + j·0
 *
 * For general memoryless nonlinearities:
 *   N(A) = (1/(πA)) ∫₀^{2π} f(A sin θ) e^{-jθ} dθ
 *
 * This function uses numerical quadrature (composite Simpson's rule)
 * to evaluate the integral.
 *
 * Complexity: O(n_quad) where n_quad is the number of quadrature points.
 *
 * @param nl       Memoryless nonlinearity descriptor.
 * @param A        Input amplitude (A > 0).
 * @param n_quad   Number of quadrature points (must be even, ≥ 100 recommended).
 * @param df       Output: describing function result.
 * @return 0 on success, -1 on error (A ≤ 0, NULL pointers, etc.).
 */
int nl_describing_function(const nl_nonlinearity_t *nl, double A,
                           int n_quad, nl_describing_function_t *df);

/**
 * @brief Compute the analytical describing function for ideal relay (L5).
 *
 * For ideal relay: f(x) = M·sgn(x), the describing function is:
 *   N(A) = 4M / (π A)
 *
 * This is a purely real, amplitude-dependent gain.
 *
 * Reference: Slotine & Li, §5.2.
 *
 * @param M   Relay output magnitude (M > 0).
 * @param A   Input amplitude (A > 0).
 * @return N(A) value (real).
 */
double nl_df_ideal_relay(double M, double A);

/**
 * @brief Compute the analytical describing function for dead-zone (L5).
 *
 * Dead-zone: f(x) = 0 for |x| ≤ δ, f(x) = x - δ·sgn(x) for |x| > δ.
 *
 * For A > δ:
 *   N(A) = 1 - (2/π)[arcsin(δ/A) + (δ/A)√(1 - (δ/A)²)]
 *
 * For A ≤ δ: N(A) = 0.
 *
 * @param delta  Dead-zone half-width (δ > 0).
 * @param A      Input amplitude.
 * @return N(A) value (real).
 */
double nl_df_dead_zone(double delta, double A);

/**
 * @brief Compute the analytical describing function for saturation (L5).
 *
 * Saturation: f(x) = kx for |x| ≤ a/k, f(x) = a·sgn(x) for |x| > a/k.
 *
 * For A > a/k:
 *   N(A) = (2k/π)[arcsin(a/(kA)) + (a/(kA))√(1 - (a/(kA))²)]
 *
 * For A ≤ a/k: N(A) = k (linear region).
 *
 * @param k   Linear gain (k > 0).
 * @param a   Saturation level (a > 0).
 * @param A   Input amplitude.
 * @return N(A) value (real).
 */
double nl_df_saturation(double k, double a, double A);

/**
 * @brief Compute the analytical describing function for backlash (L5).
 *
 * Backlash: hysteresis loop with gap width 2δ and slope k.
 *
 * For A > δ:
 *   N(A) = k/2 - (k/π)[arcsin(1 - 2δ/A) + 2(1 - 2δ/A)√((δ/A)(1 - δ/A))]
 *          + j·(4kδ/πA)(δ/A - 1)
 *
 * Note: Backlash has a phase lag (non-zero imaginary part).
 *
 * @param k       Slope of the hysteresis loop.
 * @param delta   Half gap width.
 * @param A       Input amplitude.
 * @param real_part Output: real part of N(A).
 * @param imag_part Output: imaginary part of N(A).
 * @return 0 on success.
 */
int nl_df_backlash(double k, double delta, double A,
                   double *real_part, double *imag_part);

/**
 * @brief Compute the describing function for relay with hysteresis (L5).
 *
 * Relay with hysteresis:
 *
 *   N(A) = (4M/(πA))·√(1 - (h/A)²) - j·(4Mh/(πA²))   for A ≥ h
 *   N(A) = 0                                            for A < h
 *
 * where M is the relay output magnitude and h is half the hysteresis width.
 *
 * @param M    Relay output magnitude.
 * @param h    Half hysteresis width.
 * @param A    Input amplitude.
 * @param real_part Output: real part.
 * @param imag_part Output: imaginary part.
 * @return 0 on success.
 */
int nl_df_relay_hysteresis(double M, double h, double A,
                           double *real_part, double *imag_part);

/**
 * @brief Predict limit cycle using describing function method (L6).
 *
 * The limit cycle condition is: G(jω)·N(A) = -1, or equivalently:
 *   G(jω) = -1/N(A)
 *
 * This function solves for the intersection of G(jω) and -1/N(A) loci.
 * It uses a grid search followed by Newton refinement.
 *
 * @param lin     Linear part of the system (G(s) = C(sI-A)⁻¹B + D).
 * @param nl      Nonlinearity descriptor.
 * @param A_range Amplitude search range [A_min, A_max].
 * @param nA      Number of amplitude points for grid search.
 * @param w_range Frequency search range [ω_min, ω_max] (rad/s).
 * @param nw      Number of frequency points for grid search.
 * @param A_pred  Output: predicted limit cycle amplitude.
 * @param w_pred  Output: predicted limit cycle frequency.
 * @param stable  Output: 1 if stable limit cycle, 0 if unstable.
 * @return Number of limit cycle solutions found (0, 1, or -1 on error).
 */
int nl_df_limit_cycle(const nl_linear_system_t *lin,
                      const nl_nonlinearity_t *nl,
                      const double *A_range, size_t nA,
                      const double *w_range, size_t nw,
                      double *A_pred, double *w_pred,
                      int *stable);

/**
 * @brief Compute -1/N(A) locus for stability analysis (L6).
 *
 * The -1/N(A) locus is the critical locus for the describing function
 * method. Intersection with the Nyquist plot of G(jω) indicates
 * a potential limit cycle.
 *
 * @param nl   Nonlinearity descriptor.
 * @param A    Array of amplitudes.
 * @param nA   Number of amplitude points.
 * @param locus_real Output: real parts of -1/N(A).
 * @param locus_imag Output: imaginary parts of -1/N(A).
 * @return 0 on success, -1 on error.
 */
int nl_df_critical_locus(const nl_nonlinearity_t *nl,
                         const double *A, size_t nA,
                         double *locus_real, double *locus_imag);

/**
 * @brief Harmonic balance: solve for steady-state periodic solution (L5, L8).
 *
 * Extends the DF method to include multiple harmonics for higher accuracy.
 * Approximates the periodic solution as:
 *   x(t) ≈ a₀ + Σₖ [aₖ cos(kωt) + bₖ sin(kωt)]
 *
 * and solves for the Fourier coefficients using Galerkin projection.
 *
 * Complexity: O(N_harmonics³ · N_iterations).
 *
 * @param sys          Nonlinear system descriptor.
 * @param N_harmonics  Number of harmonics to include.
 * @param w_guess      Initial frequency guess (rad/s).
 * @param coeffs       Output: Fourier coefficients (length 2*N_harmonics+1).
 * @param w_sol        Output: converged frequency.
 * @return 0 on success, -1 on non-convergence.
 */
int nl_harmonic_balance(nl_nonlinear_system_t *sys,
                        int N_harmonics, double w_guess,
                        double *coeffs, double *w_sol);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_DESCRIBING_H */
