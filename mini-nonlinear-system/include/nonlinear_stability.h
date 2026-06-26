/**
 * @file nonlinear_stability.h
 * @brief Lyapunov stability analysis for nonlinear systems.
 *
 * Implements Lyapunov's direct method (second method), the indirect
 * method (first method/linearization), and related stability criteria
 * including the circle criterion and Popov criterion.
 *
 * L2: Core Concepts — stability of equilibria, Lyapunov functions
 * L4: Fundamental Laws — Lyapunov stability theorems, Circle criterion,
 *     Popov criterion, Invariance principle (LaSalle)
 * L5: Algorithms — Lyapunov function construction, linear matrix inequality
 *     (LMI) approach for quadratic Lyapunov functions.
 *
 * Course Mapping:
 *   - MIT 6.241: Dynamic Systems and Control (Lyapunov theory)
 *   - Stanford ENGR 209A: Analysis and Control of Nonlinear Systems
 *   - ETH 227-0216: Nonlinear Systems and Control
 *   - Berkeley EE 222: Nonlinear Systems
 *
 * Reference:
 *   - Khalil, "Nonlinear Systems", Ch.3-7, Prentice Hall, 2002
 *   - Slotine & Li, "Applied Nonlinear Control", Prentice Hall, 1991
 *   - Vidyasagar, "Nonlinear Systems Analysis", SIAM, 2002
 */

#ifndef NONLINEAR_STABILITY_H
#define NONLINEAR_STABILITY_H

#include "nonlinear_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L2, L4: Lyapunov Stability Analysis
 * -------------------------------------------------------------------------*/

/**
 * @brief Evaluate the state derivative for an autonomous system.
 *
 * Computes dx/dt = f(x) for the given nonlinear system at state x.
 * This is the fundamental operation for all stability analysis.
 *
 * Complexity: O(dim * cost(f)) where cost(f) depends on the system.
 *
 * @param sys  Nonlinear system descriptor.
 * @param x    Current state (length = sys->dim).
 * @param dxdt Output: state derivative (length = sys->dim, pre-allocated).
 * @return 0 on success, -1 on null pointer or dimension mismatch.
 */
int nl_evaluate_f(nl_nonlinear_system_t *sys,
                  const double *x, double *dxdt);

/**
 * @brief Evaluate a quadratic Lyapunov function V(x) = xᵀ P x (L4).
 *
 * The matrix P must be symmetric positive definite. This is the
 * most common form of Lyapunov function for linear and linearized systems.
 *
 * Reference: Lyapunov's direct method (1892).
 *
 * @param P    Symmetric matrix (row-major, dim x dim).
 * @param x    State vector (length = dim).
 * @param dim  State dimension.
 * @return V(x) = xᵀ P x. Returns 0 if P or x is NULL.
 */
double nl_quadratic_lyapunov(const double *P, const double *x, size_t dim);

/**
 * @brief Compute V_dot = ∇V · f(x) for a quadratic Lyapunov function (L4).
 *
 * For V(x) = xᵀ P x, we have V_dot = 2 xᵀ P f(x) = 2 Σ_i Σ_j x_i P_ij f_j(x).
 *
 * Theorem (Lyapunov): If V(x) > 0 and V_dot(x) ≤ 0 for all x ≠ 0 near
 * the origin, then the origin is stable. If V_dot(x) < 0, it is
 * asymptotically stable.
 *
 * @param P     Symmetric positive definite matrix (row-major).
 * @param x     State vector.
 * @param sys   Nonlinear system (provides f(x)).
 * @param dim   State dimension.
 * @return V_dot value, or 0 if any pointer is NULL.
 */
double nl_quadratic_lyapunov_dot(const double *P, const double *x,
                                  nl_nonlinear_system_t *sys,
                                  size_t dim);

/**
 * @brief Compute the Jacobian matrix of f at a given point (L3, L4).
 *
 * Uses central finite differences:
 *   J_{ij} ≈ (f_i(x + h e_j) - f_i(x - h e_j)) / (2h)
 *
 * This is the core operation for Lyapunov's indirect method:
 * the eigenvalues of J at the equilibrium determine local stability.
 *
 * Theorem (Lyapunov's Indirect Method): Let x* be an equilibrium.
 * If all eigenvalues of J(x*) have negative real parts, x* is
 * asymptotically stable. If any eigenvalue has positive real part,
 * x* is unstable.
 *
 * Complexity: O(dim² * cost(f)).
 *
 * @param sys  Nonlinear system.
 * @param x    Point at which to compute Jacobian.
 * @param J    Output: Jacobian matrix (row-major, dim x dim, pre-allocated).
 * @return 0 on success, -1 on error.
 */
int nl_compute_jacobian(nl_nonlinear_system_t *sys,
                        const double *x, nl_jacobian_t *J);

/**
 * @brief Classify stability of an equilibrium using eigenvalues of J (L4).
 *
 * Implements Lyapunov's first method (indirect method) for hyperbolic
 * equilibria. Examines the real parts of eigenvalues to classify
 * the equilibrium as stable node, saddle, focus, center, etc.
 *
 * Reference: Hirsch & Smale, "Differential Equations, Dynamical Systems,
 * and Linear Algebra", Academic Press, 1974.
 *
 * @param J     Jacobian at the equilibrium (populated by nl_compute_jacobian).
 * @param tol   Tolerance for zero check (e.g., 1e-6).
 * @return Stability classification.
 */
nl_stability_t nl_classify_equilibrium(const nl_jacobian_t *J, double tol);

/**
 * @brief Compute eigenvalues of a real n×n matrix using QR algorithm (L3, L5).
 *
 * Implements the Francis QR double-shift algorithm for real matrices.
 * Handles both real and complex conjugate eigenvalue pairs.
 *
 * Complexity: O(n³).
 *
 * @param A         Input matrix (row-major, n x n).
 * @param n         Matrix dimension.
 * @param real_part Output: real parts of eigenvalues.
 * @param imag_part Output: imaginary parts of eigenvalues.
 * @param max_iter  Maximum iterations for QR convergence.
 * @param tol       Convergence tolerance.
 * @return 0 on success, -1 on non-convergence.
 */
int nl_eigenvalues(double *A, size_t n, double *real_part, double *imag_part,
                   int max_iter, double tol);

/**
 * @brief Solve the Lyapunov equation Aᵀ P + P A = -Q for P (L4, L5).
 *
 * This is the continuous-time algebraic Lyapunov equation.
 * If for any Q > 0 there exists P > 0, then A is Hurwitz (stable).
 *
 * Uses the Bartels-Stewart algorithm via real Schur decomposition.
 * The solution P serves as the quadratic Lyapunov function matrix.
 *
 * Complexity: O(n³).
 *
 * @param A  System matrix (row-major, n x n).
 * @param Q  Right-hand side matrix (symmetric positive definite).
 * @param P  Output: solution matrix (row-major, n x n, pre-allocated).
 * @param n  Matrix dimension.
 * @return 0 on success, -1 if no solution or numerical issues.
 */
int nl_solve_lyapunov_equation(const double *A, const double *Q,
                                double *P, size_t n);

/**
 * @brief Apply Lyapunov's direct method to check stability (L4).
 *
 * Constructs a quadratic Lyapunov function V(x) = xᵀ P x by solving
 * the Lyapunov equation with the linearization at the equilibrium,
 * then verifies V_dot < 0 in a neighborhood.
 *
 * @param sys         Nonlinear system.
 * @param equilibrium Equilibrium point.
 * @param P           Output: Lyapunov matrix (if found).
 * @param radius      Output: estimated region of attraction radius.
 * @return NL_ASYMPTOTICALLY_STABLE if proven stable,
 *         NL_UNSTABLE if proven unstable,
 *         NL_NONHYPERBOLIC if inconclusive.
 */
nl_stability_t nl_direct_method(nl_nonlinear_system_t *sys,
                                 const nl_state_t *equilibrium,
                                 double *P, double *radius);

/**
 * @brief Circle criterion for Lur'e systems: absolute stability test (L4, L6).
 *
 * Theorem (Circle Criterion): For a Lur'e system with nonlinearity φ(y)
 * in sector [k₁, k₂], if the Nyquist plot of the linear part G(s) does
 * not intersect the critical disk D(-1/k₁, -1/k₂) and encircles it
 * the correct number of times, the system is absolutely stable.
 *
 * This function checks the frequency-domain condition:
 *   Re[ (1 + k₂ G(jω)) / (1 + k₁ G(jω)) ] > 0   for all ω.
 *
 * Reference: Zames, "On the input-output stability of time-varying
 * nonlinear feedback systems", IEEE TAC, 1966.
 *
 * @param lure    Lur'e system to analyze.
 * @param freqs   Frequency points for evaluation (rad/s).
 * @param nfreqs  Number of frequency points.
 * @param margin  Output: stability margin (minimum of circle condition).
 * @return 1 if absolutely stable, 0 if test is inconclusive, -1 on error.
 */
int nl_circle_criterion(const nl_lure_system_t *lure,
                        const double *freqs, size_t nfreqs,
                        double *margin);

/**
 * @brief Popov criterion for Lur'e systems (L4, L6).
 *
 * Theorem (Popov): For a Lur'e system with nonlinearity φ(y) in sector
 * [0, k], if there exists η ≥ 0 such that
 *   Re[(1 + j ω η) G(jω)] + 1/k ≥ ε > 0   for all ω,
 * then the origin is globally asymptotically stable.
 *
 * This is less conservative than the circle criterion for time-invariant
 * nonlinearities.
 *
 * Reference: Popov, "Absolute stability of nonlinear systems of automatic
 * control", Automation and Remote Control, 1961.
 *
 * @param lure    Lur'e system.
 * @param freqs   Frequency vector (rad/s).
 * @param nfreqs  Number of frequencies.
 * @param eta     Output: optimal Popov multiplier η*.
 * @param margin  Output: stability margin.
 * @return 1 if Popov-stable, 0 if inconclusive, -1 on error.
 */
int nl_popov_criterion(const nl_lure_system_t *lure,
                       const double *freqs, size_t nfreqs,
                       double *eta, double *margin);

/**
 * @brief Compute transfer function G(jω) = C(jωI - A)⁻¹ B + D (L3).
 *
 * Evaluates the frequency response of the linear part of a Lur'e system.
 *
 * @param lin   Linear system descriptor.
 * @param omega Frequency in rad/s.
 * @param re    Output: real part of G(jω).
 * @param im    Output: imaginary part of G(jω).
 * @return 0 on success, -1 on error.
 */
int nl_transfer_function(const nl_linear_system_t *lin, double omega,
                         double *re, double *im);

/**
 * @brief Check if a matrix P is positive definite (L3).
 *
 * Uses Cholesky decomposition attempt: if Cholesky succeeds,
 * P is positive definite. Returns the minimum eigenvalue estimate
 * via Gershgorin circle theorem for diagnostic purposes.
 *
 * @param P    Matrix (row-major, n x n).
 * @param n    Dimension.
 * @param min_eig Output: estimated minimum eigenvalue.
 * @return 1 if positive definite, 0 otherwise.
 */
int nl_is_positive_definite(const double *P, size_t n, double *min_eig);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_STABILITY_H */
