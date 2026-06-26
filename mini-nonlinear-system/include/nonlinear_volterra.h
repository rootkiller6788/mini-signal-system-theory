/**
 * @file nonlinear_volterra.h
 * @brief Volterra and Wiener series for nonlinear system identification.
 *
 * Volterra series extend the concept of linear convolution to
 * nonlinear systems, providing a polynomial functional expansion
 * that is universal for time-invariant nonlinear systems with
 * fading memory.
 *
 * L1: Definitions — Volterra kernel, Wiener kernel, G-functional
 * L3: Math Structures — functional power series, multilinear operators,
 *     symmetric kernels, frequency-domain kernels
 * L5: Algorithms — kernel identification via cross-correlation
 *     (Lee-Schetzen), orthogonal Wiener series via Gram-Schmidt
 *     with Gaussian white noise input
 * L8: Advanced Topics — sparse Volterra, frequency-domain Volterra,
 *     nonlinear frequency response functions (GFRF)
 *
 * Course Mapping:
 *   - MIT 6.435: System Identification (nonlinear methods)
 *   - Tsinghua: 系统辨识 (Volterra/Wiener methods)
 *
 * Reference:
 *   - Rugh, "Nonlinear System Theory: The Volterra/Wiener Approach",
 *     Johns Hopkins, 1981
 *   - Schetzen, "The Volterra and Wiener Theories of Nonlinear Systems",
 *     Wiley, 1980 (reprint Krieger, 2006)
 */

#ifndef NONLINEAR_VOLTERRA_H
#define NONLINEAR_VOLTERRA_H

#include "nonlinear_core.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L1, L5: Volterra Kernel Operations
 * -------------------------------------------------------------------------*/

/**
 * @brief Allocate and initialize a single Volterra kernel (L1).
 *
 * For a kernel of order n with memory length M in each dimension:
 *   h_n[τ₁, ..., τ_n],   τ_i ∈ {0, 1, ..., M-1}
 *
 * The kernel is stored as a flattened array with dims = (M, M, ..., M)
 * (n copies). For causal and symmetric kernels, the effective number
 * of parameters is reduced.
 *
 * @param kernel  Kernel to initialize.
 * @param order   Kernel order n ≥ 1.
 * @param mem_len Memory length M for each τ dimension.
 * @return 0 on success, -1 on allocation failure.
 */
int nl_kernel_init(nl_volterra_kernel_t *kernel, size_t order,
                   size_t mem_len);

/**
 * @brief Free a single Volterra kernel.
 */
void nl_kernel_free(nl_volterra_kernel_t *kernel);

/**
 * @brief Set a specific kernel coefficient (L1).
 *
 * Sets h_n[τ₁, ..., τ_n] = value.
 * Automatic symmetrization can be requested for symmetric kernels.
 *
 * @param kernel  Kernel.
 * @param indices Array of n indices (τ₁, ..., τ_n), each in [0, M-1].
 * @param value   Coefficient value.
 * @return 0 on success, -1 if indices out of bounds.
 */
int nl_kernel_set(nl_volterra_kernel_t *kernel, const size_t *indices,
                  double value);

/**
 * @brief Get a specific kernel coefficient (L1).
 *
 * @param kernel  Kernel.
 * @param indices Array of n indices.
 * @param value   Output: coefficient value.
 * @return 0 on success.
 */
int nl_kernel_get(const nl_volterra_kernel_t *kernel, const size_t *indices,
                  double *value);

/**
 * @brief Compute the L² norm of a kernel (L3).
 *
 * ||h_n||² = Σ_{τ₁} ... Σ_{τ_n} h_n[τ₁,...,τ_n]²
 *
 * For symmetric kernels, the sum includes multiplicity.
 *
 * @param kernel Kernel.
 * @return L² norm squared.
 */
double nl_kernel_norm_sq(const nl_volterra_kernel_t *kernel);

/* ---------------------------------------------------------------------------
 * L5: Volterra Series Computation
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute the Volterra series output for a given input sequence (L5).
 *
 * y[t] = h₀ + Σ_{n=1}^{M} Σ_{τ₁=0}^{m-1} ... Σ_{τ_n=0}^{m-1}
 *        h_n[τ₁,...,τ_n] · u[t-τ₁] · ... · u[t-τ_n]
 *
 * Complexity: O(M · m^M) — exponential in max order M,
 * making this feasible only for low-order (M ≤ 3) or small memory.
 *
 * For higher orders, use frequency-domain methods or sparse Volterra.
 *
 * @param vs       Volterra series model.
 * @param input    Input sequence.
 * @param input_len Length of input.
 * @param output   Output array (pre-allocated, length input_len).
 * @return 0 on success.
 */
int nl_volterra_evaluate(const nl_volterra_series_t *vs,
                         const double *input, size_t input_len,
                         double *output);

/**
 * @brief Evaluate a single term (order n) of the Volterra series (L5).
 *
 * y_n[t] = Σ_{τ₁} ... Σ_{τ_n} h_n[τ₁,...,τ_n] · u[t-τ₁] · ... · u[t-τ_n]
 *
 * This allows analyzing the contribution of each order separately.
 *
 * @param vs       Volterra series.
 * @param order    Kernel order to evaluate (1 ≤ order ≤ max_order).
 * @param input    Input sequence.
 * @param input_len Length.
 * @param output_n Output: n-th order contribution.
 * @return 0 on success.
 */
int nl_volterra_evaluate_order(const nl_volterra_series_t *vs,
                                size_t order, const double *input,
                                size_t input_len, double *output_n);

/* ---------------------------------------------------------------------------
 * L5, L8: Wiener Series (Orthogonal Expansion)
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute Wiener G-functionals via Gram-Schmidt orthogonalization (L5, L8).
 *
 * Wiener G-functionals are orthogonal with respect to a Gaussian white
 * noise input with power spectral density P. This orthogonality
 * simplifies system identification.
 *
 * G₀ = k₀
 * G₁[k₁; u] = Σ h₁(τ) u(t-τ)
 * G₂[k₂; u] = ΣΣ h₂(τ₁,τ₂) u(t-τ₁)u(t-τ₂) - P Σ h₂(τ₁,τ₁)
 * G₃[k₃; u] = ΣΣΣ h₃(τ₁,τ₂,τ₃) u(t-τ₁)u(t-τ₂)u(t-τ₃)
 *              - 3P Σ h₃(τ₁,τ₁,τ₂) u(t-τ₂)
 * etc.
 *
 * Reference: Wiener, "Nonlinear Problems in Random Theory", MIT Press, 1958.
 *
 * @param vs         Volterra/Wiener series (kernels in vs).
 * @param input      Input sequence (should be Gaussian white noise
 *                   for proper orthogonality).
 * @param input_len  Length.
 * @param power      Input power spectral density P.
 * @param G_output   Output: array of G-functional outputs
 *                   (max_order × input_len, row-major).
 * @return 0 on success.
 */
int nl_wiener_G_functionals(const nl_volterra_series_t *vs,
                            const double *input, size_t input_len,
                            double power, double *G_output);

/* ---------------------------------------------------------------------------
 * L5, L8: System Identification via Cross-Correlation
 * -------------------------------------------------------------------------*/

/**
 * @brief Identify first-order (linear) kernel via cross-correlation (L5).
 *
 * h₁[τ] = (1/P) · E[y[t] · u[t-τ]]
 *
 * For Gaussian white noise input with variance P.
 *
 * @param input      Gaussian white noise input.
 * @param output     System output.
 * @param length     Data length.
 * @param mem_len    Memory length of kernel.
 * @param power      Input power P.
 * @param kernel     Output: identified linear kernel.
 * @return 0 on success.
 */
int nl_identify_kernel_1(const double *input, const double *output,
                         size_t length, size_t mem_len,
                         double power, nl_volterra_kernel_t *kernel);

/**
 * @brief Identify second-order kernel via Lee-Schetzen method (L5, L8).
 *
 * For Wiener G-functionals with Gaussian white noise input:
 *   h₂[τ₁, τ₂] = (1/(2! P²)) · E[ (y[t] - G₀ - G₁[t]) · u[t-τ₁] u[t-τ₂] ]
 *
 * where G₀ and G₁ are the identified lower-order contributions.
 *
 * This function implements the cross-correlation estimation using
 * the diagonal and off-diagonal terms separately to handle symmetry.
 *
 * Reference: Lee & Schetzen, "Measurement of the Wiener kernels of a
 * nonlinear system by cross-correlation", Int. J. Control, 1965.
 *
 * @param input     Gaussian white noise input.
 * @param output    System output.
 * @param length    Data length.
 * @param mem_len   Memory length.
 * @param power     Input power P.
 * @param kernel_1  Previously identified first-order kernel.
 * @param kernel_2  Output: identified second-order kernel.
 * @return 0 on success.
 */
int nl_identify_kernel_2(const double *input, const double *output,
                         size_t length, size_t mem_len,
                         double power,
                         const nl_volterra_kernel_t *kernel_1,
                         nl_volterra_kernel_t *kernel_2);

/**
 * @brief Identify Volterra series up to given order (L8).
 *
 * Sequentially identifies kernels of orders 1 through max_order
 * using the cross-correlation method with Gaussian white noise input.
 *
 * @param input     Gaussian white noise.
 * @param output    System output.
 * @param length    Data length.
 * @param mem_len   Memory length per dimension.
 * @param power     Input power.
 * @param vs        Output: identified Volterra series (pre-initialized).
 * @return 0 on success.
 */
int nl_identify_volterra_series(const double *input, const double *output,
                                 size_t length, size_t mem_len,
                                 double power, nl_volterra_series_t *vs);

/* ---------------------------------------------------------------------------
 * L3, L8: Frequency-Domain Volterra (Generalized FRF)
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute the Generalized Frequency Response Function (GFRF) (L8).
 *
 * The n-th order GFRF H_n(ω₁, ..., ω_n) is the n-dimensional Fourier
 * transform of the n-th order Volterra kernel:
 *
 *   H_n(ω₁,...,ω_n) = Σ_{τ₁}...Σ_{τ_n}
 *       h_n[τ₁,...,τ_n] · e^{-j(ω₁τ₁ + ... + ω_nτ_n)}
 *
 * This is the frequency-domain representation that generalizes
 * the linear transfer function to nonlinear systems.
 *
 * @param kernel    Volterra kernel.
 * @param freqs     Frequency vector (length n_freqs).
 * @param n_freqs   Number of frequency points.
 * @param GFRF      Output: H_n(ω₁, ..., ω_n) as complex numbers
 *                  (n_freqs^n elements, multidimensional).
 * @param n_freq_dim Array of n copies of n_freqs for dimension sizes.
 * @return 0 on success.
 */
int nl_compute_GFRF(const nl_volterra_kernel_t *kernel,
                    const double *freqs, size_t n_freqs,
                    double complex *GFRF, size_t *n_freq_dim);

/**
 * @brief Compute output spectrum using the GFRFs (nonlinear frequency response) (L8).
 *
 * For an M-th order Volterra system with sinusoidal inputs at
 * frequencies ω_1, ..., ω_K, the output spectrum contains intermodulation
 * products at frequencies:
 *
 *   Ω = n₁ω₁ + n₂ω₂ + ... + n_K ω_K
 *
 * where Σ |n_i| ≤ M.
 *
 * This function predicts the output spectrum.
 *
 * @param vs          Volterra series.
 * @param freqs_input Input frequencies.
 * @param amps_input  Input amplitudes.
 * @param n_input     Number of input tones.
 * @param output_freqs Output frequencies.
 * @param output_amps  Output amplitudes.
 * @param max_output   Max output tones.
 * @return Number of output tones found.
 */
int nl_output_spectrum(const nl_volterra_series_t *vs,
                       const double *freqs_input,
                       const double *amps_input,
                       size_t n_input,
                       double *output_freqs,
                       double *output_amps,
                       size_t max_output);

/**
 * @brief Symmetrize a Volterra kernel (L3).
 *
 * A Volterra kernel can always be symmetrized without changing the
 * input-output behavior:
 *
 *   h_n^{sym}[τ₁,...,τ_n] = (1/n!) Σ_{π∈S_n} h_n[τ_{π(1)},...,τ_{π(n)}]
 *
 * where S_n is the symmetric group of permutations.
 *
 * Symmetric kernels have n! times fewer independent parameters.
 *
 * @param kernel Kernel to symmetrize (modified in place).
 * @return 0 on success.
 */
int nl_kernel_symmetrize(nl_volterra_kernel_t *kernel);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_VOLTERRA_H */
