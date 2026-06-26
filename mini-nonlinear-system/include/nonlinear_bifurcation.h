/**
 * @file nonlinear_bifurcation.h
 * @brief Bifurcation analysis for nonlinear dynamical systems.
 *
 * Bifurcation theory studies qualitative changes in system behavior
 * as parameters vary. This module implements numerical methods for
 * detecting and classifying bifurcations in ODE systems.
 *
 * L2: Core Concepts — bifurcation types, stability exchange,
 *     normal forms, codimension
 * L5: Algorithms — numerical continuation, branch switching,
 *     bifurcation detection (test functions), normal form computation
 * L6: Canonical Problems — saddle-node, pitchfork, Hopf bifurcations
 *
 * Course Mapping:
 *   - MIT 18.385: Nonlinear Dynamics and Chaos (bifurcation theory)
 *   - Cornell MAE 6780: Nonlinear Dynamics
 *   - Tsinghua: 非线性动力学 (bifurcation and chaos)
 *
 * Reference:
 *   - Kuznetsov, "Elements of Applied Bifurcation Theory", 3rd ed.,
 *     Springer, 2004
 *   - Seydel, "Practical Bifurcation and Stability Analysis", 3rd ed.,
 *     Springer, 2010
 *   - Guckenheimer & Holmes, "Nonlinear Oscillations, Dynamical Systems,
 *     and Bifurcations of Vector Fields", Springer, 1983
 */

#ifndef NONLINEAR_BIFURCATION_H
#define NONLINEAR_BIFURCATION_H

#include "nonlinear_core.h"
#include "nonlinear_stability.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L5: Numerical Continuation and Bifurcation Detection
 * -------------------------------------------------------------------------*/

/**
 * @brief Find an equilibrium point using Newton-Raphson iteration (L5).
 *
 * Solves f(x, μ) = 0 for x given parameter μ.
 * The Newton iteration is:
 *   x_{k+1} = x_k - J(x_k)⁻¹ f(x_k)
 *
 * where J is the Jacobian of f with respect to x.
 *
 * Complexity: O(dim³ · iterations) due to linear solve per step.
 *
 * @param sys    Nonlinear system.
 * @param x0     Initial guess (length = dim).
 * @param params Parameter values (copied to sys->params).
 * @param nparam Number of parameters.
 * @param x_eq   Output: equilibrium point.
 * @param tol    Convergence tolerance for ||f(x)||.
 * @return 0 on success, -1 on non-convergence, -2 on singular Jacobian.
 */
int nl_find_equilibrium(nl_nonlinear_system_t *sys,
                        const double *x0, const double *params,
                        size_t nparam, nl_state_t *x_eq, double tol);

/**
 * @brief Pseudo-arclength continuation: trace equilibrium curve (L5, L8).
 *
 * Follows the equilibrium manifold x*(μ) in (x, μ) space using
 * pseudo-arclength continuation. This allows passing through folds
 * (saddle-node bifurcations) where the Jacobian becomes singular
 * and naive parameter continuation would fail.
 *
 * The extended system is:
 *   f(x, μ) = 0
 *   (x - x_prev)ᵀ (dx/ds)_prev + (μ - μ_prev) (dμ/ds)_prev - Δs = 0
 *
 * Reference: Keller, "Lectures on Numerical Methods in Bifurcation
 * Problems", Springer, 1987.
 *
 * @param sys        System (uses param[0] as the active parameter).
 * @param x_start    Initial equilibrium state.
 * @param mu_start   Initial parameter value.
 * @param ds         Step size in arclength.
 * @param n_steps    Number of continuation steps.
 * @param branch_x   Output: array of equilibrium states (length n_steps).
 * @param branch_mu  Output: array of parameter values (length n_steps).
 * @return Number of successfully computed steps, -1 on failure.
 */
int nl_continuation(nl_nonlinear_system_t *sys,
                    const nl_state_t *x_start, double mu_start,
                    double ds, size_t n_steps,
                    nl_state_t *branch_x, double *branch_mu);

/**
 * @brief Detect saddle-node (fold) bifurcation using test function (L5, L6).
 *
 * At a saddle-node bifurcation, the Jacobian has a zero eigenvalue:
 *   det J(x*, μ*) = 0
 *
 * This function computes the determinant of the Jacobian as a
 * test function ψ(x, μ) = det J(x, μ). A sign change in ψ along the
 * equilibrium branch indicates a saddle-node bifurcation.
 *
 * Uses LU decomposition for determinant computation.
 *
 * @param sys          System (uses param[0] as active parameter).
 * @param branch_x     Array of equilibrium states from continuation.
 * @param branch_mu    Array of parameter values.
 * @param n_pts        Number of branch points.
 * @param bif_pts      Output: detected bifurcation points.
 * @param max_bif      Maximum number of bifurcations to detect.
 * @return Number of bifurcations detected (≥0), -1 on error.
 */
int nl_detect_saddle_node(nl_nonlinear_system_t *sys,
                          const nl_state_t *branch_x,
                          const double *branch_mu, size_t n_pts,
                          nl_bifurcation_point_t *bif_pts,
                          size_t max_bif);

/**
 * @brief Detect Andronov-Hopf bifurcation using test function (L5, L6).
 *
 * At a Hopf bifurcation, the Jacobian has a pair of purely imaginary
 * eigenvalues ±jω with nonzero speed (transversality) and no other
 * eigenvalues on the imaginary axis (non-resonance).
 *
 * The test function uses the product of eigenvalue pairs:
 * for each pair (λ, λ̄), compute λ + λ̄ = 2Re(λ). A sign change
 * indicates crossing the imaginary axis.
 *
 * Reference: Marsden & McCracken, "The Hopf Bifurcation and Its
 * Applications", Springer, 1976.
 *
 * @param sys       System.
 * @param branch_x  Equilibrium states.
 * @param branch_mu Parameter values.
 * @param n_pts     Number of points.
 * @param bif_pts   Output: detected bifurcations.
 * @param max_bif   Max bifurcations.
 * @return Number detected (≥0), -1 on error.
 */
int nl_detect_hopf(nl_nonlinear_system_t *sys,
                   const nl_state_t *branch_x,
                   const double *branch_mu, size_t n_pts,
                   nl_bifurcation_point_t *bif_pts,
                   size_t max_bif);

/**
 * @brief Compute the normal form coefficients for a Hopf bifurcation (L5, L8).
 *
 * At a Hopf bifurcation, the dynamics on the center manifold can be
 * reduced to the normal form in polar coordinates:
 *
 *   ṙ = α (μ - μ_c) r + a r³ + O(r⁵)
 *   θ̇ = ω_c + β (μ - μ_c) + b r² + O(r⁴)
 *
 * where a is the first Lyapunov coefficient:
 *   a < 0 → supercritical (stable limit cycle)
 *   a > 0 → subcritical (unstable limit cycle)
 *
 * This function computes a using center manifold reduction and
 * normalization. The computation requires the Jacobian and its
 * eigenvectors, plus second and third order derivatives of f.
 *
 * @param sys       System at the bifurcation point.
 * @param x_bif     Equilibrium at bifurcation.
 * @param mu_bif    Parameter value at bifurcation.
 * @param omega_c   Output: critical frequency ω_c.
 * @param a_coeff   Output: first Lyapunov coefficient a.
 * @return 0 on success, -1 on error.
 */
int nl_hopf_normal_form(nl_nonlinear_system_t *sys,
                        const nl_state_t *x_bif, double mu_bif,
                        double *omega_c, double *a_coeff);

/**
 * @brief Compute the determinant of a square matrix via LU decomposition (L3).
 *
 * det(A) = ∏_{i=0}^{n-1} L_{ii} · U_{ii}  (product of diagonal entries
 * of L and U, where U_{ii} has the pivots and L_{ii} = 1 for unit
 * lower triangular).
 *
 * Complexity: O(n³).
 *
 * @param A    Input matrix (row-major, n x n).
 * @param n    Dimension.
 * @param det  Output: determinant value.
 * @return 0 on success, -1 if matrix is singular (det ≈ 0).
 */
int nl_matrix_determinant(const double *A, size_t n, double *det);

/**
 * @brief Check transversality condition for bifurcation (L4, L6).
 *
 * For a bifurcation at (x*, μ*), the transversality condition requires
 * that the critical eigenvalue λ(μ) crosses the imaginary axis with
 * nonzero speed: d(Re λ)/dμ ≠ 0 at μ = μ*.
 *
 * This is computed via finite difference approximation using the
 * secular equation for the eigenvalue.
 *
 * @param sys     System.
 * @param x_bif   Equilibrium at bifurcation.
 * @param mu_bif  Parameter at bifurcation.
 * @param dmu     Perturbation for finite difference.
 * @param speed   Output: d(Re λ)/dμ at bifurcation.
 * @return 0 on success, nonzero if transversality fails.
 */
int nl_check_transversality(nl_nonlinear_system_t *sys,
                             const nl_state_t *x_bif,
                             double mu_bif, double dmu,
                             double *speed);

/**
 * @brief Compute the bifurcation diagram data for a 1D parameter sweep (L6).
 *
 * Generates a bifurcation diagram by sweeping a parameter and
 * computing the asymptotic behavior (fixed points, periodic orbits,
 * or chaos) for each parameter value.
 *
 * For each parameter value μ_k, the system is integrated from an
 * initial condition for T_transient time (to eliminate transients),
 * then the local extrema of the state are recorded for T_record time.
 *
 * @param sys          System (param[0] swept).
 * @param x0           Initial condition.
 * @param mu_range     [μ_min, μ_max].
 * @param n_mu         Number of parameter values.
 * @param T_transient  Transient integration time.
 * @param T_record     Recording time.
 * @param dt           Integration step.
 * @param diagram_x    Output: local extrema (2D array, n_mu × max_pts).
 * @param n_pts_per_mu Output: number of extrema per μ.
 * @param max_pts      Maximum extrema per parameter value.
 * @return 0 on success.
 */
int nl_bifurcation_diagram(nl_nonlinear_system_t *sys,
                           const double *x0,
                           const double *mu_range, size_t n_mu,
                           double T_transient, double T_record,
                           double dt, double *diagram_x,
                           size_t *n_pts_per_mu, size_t max_pts);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_BIFURCATION_H */
