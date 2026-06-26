/**
 * @file nonlinear_phase.h
 * @brief Phase plane / phase space analysis for 2D nonlinear systems.
 *
 * Provides tools for constructing and analyzing phase portraits of
 * planar (2D) nonlinear autonomous systems:
 *   dx/dt = f₁(x, y)
 *   dy/dt = f₂(x, y)
 *
 * Phase plane analysis is the graphical method for understanding
 * global behavior of second-order nonlinear systems. The Poincaré-Bendixson
 * theorem guarantees that bounded trajectories in the plane can only
 * approach fixed points, limit cycles, or heteroclinic/homoclinic orbits.
 *
 * L2: Core Concepts — phase plane, nullclines, trajectories, separatrices
 * L3: Math Structures — vector fields, flow maps
 * L4: Fundamental Laws — Poincaré-Bendixson theorem, Bendixson's criterion
 *     for non-existence of limit cycles
 * L5: Algorithms — phase portrait generation, nullcline computation,
 *     equilibrium finding, separatrix tracing
 * L6: Canonical Problems — Van der Pol oscillator, Duffing oscillator,
 *     pendulum, predator-prey (Lotka-Volterra), FitzHugh-Nagumo
 *
 * Course Mapping:
 *   - MIT 18.03: Differential Equations (phase plane)
 *   - MIT 6.003: Signals and Systems (nonlinear systems)
 *   - Berkeley EE 16A/B: phase plane for circuits
 *   - Tsinghua: 信号与系统 (phase plane for nonlinear systems)
 *
 * Reference:
 *   - Strogatz, "Nonlinear Dynamics and Chaos", Ch.5-7, 1994
 *   - Guckenheimer & Holmes, "Nonlinear Oscillations...", Ch.1, 1983
 */

#ifndef NONLINEAR_PHASE_H
#define NONLINEAR_PHASE_H

#include "nonlinear_core.h"
#include "nonlinear_stability.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L5: Phase Portrait Generation
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute the phase portrait for a 2D system (L5, L6).
 *
 * Populates the nl_phase_portrait_t structure with:
 * - Vector field on a regular grid (dx/dt, dy/dt at each point)
 * - Nullclines (curves where dx/dt = 0 or dy/dt = 0)
 * - Equilibrium points (intersections of nullclines)
 * - Representative trajectories from selected initial conditions
 *
 * The vector field grid enables visualization via quiver plots.
 * Nullclines partition the phase plane into regions where the
 * vector field has consistent sign patterns.
 *
 * @param sys  2D nonlinear system (dim = 2 required).
 * @param pp   Phase portrait structure (already initialized via
 *             nl_phase_portrait_init).
 * @return 0 on success, -1 if system is not 2D.
 */
int nl_compute_phase_portrait(nl_nonlinear_system_t *sys,
                               nl_phase_portrait_t *pp);

/**
 * @brief Find all equilibrium points in a region using nullcline intersection (L5).
 *
 * For 2D systems, equilibria are intersections of the x-nullcline
 * (f₁ = 0) and y-nullcline (f₂ = 0). This function binary-searches
 * along nullclines to locate intersection points, then refines
 * using Newton's method.
 *
 * @param sys  2D system.
 * @param pp   Phase portrait (nullclines already computed).
 * @param tol  Tolerance for f₁ = f₂ = 0.
 * @return Number of equilibria found, -1 on error.
 */
int nl_find_equilibria_2d(nl_nonlinear_system_t *sys,
                          nl_phase_portrait_t *pp, double tol);

/**
 * @brief Compute a single trajectory via Runge-Kutta 4(5) (Dormand-Prince) (L5).
 *
 * Integrates dx/dt = f(x) from x0 for time T using adaptive step size
 * Dormand-Prince method (RK5(4)7M by default for tolerances).
 *
 * This function uses a fixed-step classical RK4 for simplicity and
 * predictability in phase portrait generation. For high-precision
 * needs, use nl_integrate_adaptive.
 *
 * Complexity: O(T/dt · cost(f)).
 *
 * @param sys       System.
 * @param x0        Initial state (length = dim).
 * @param T         Total integration time.
 * @param dt        Fixed step size.
 * @param traj      Output: trajectory (already initialized).
 * @param n_steps   Output: number of steps actually computed.
 * @return 0 on success.
 */
int nl_integrate_rk4(nl_nonlinear_system_t *sys,
                     const double *x0, double T, double dt,
                     nl_trajectory_t *traj, size_t *n_steps);

/**
 * @brief Adaptive-step Dormand-Prince 5(4) integration (L5).
 *
 * Higher-order adaptive integration with error control.
 *
 * @param sys    System.
 * @param x0     Initial state.
 * @param T      Total time.
 * @param dt0    Initial step size.
 * @param tol    Local error tolerance per step.
 * @param traj   Output: trajectory.
 * @return 0 on success.
 */
int nl_integrate_adaptive(nl_nonlinear_system_t *sys,
                          const double *x0, double T, double dt0,
                          double tol, nl_trajectory_t *traj);

/* ---------------------------------------------------------------------------
 * L4: Poincaré-Bendixson Theorem and Bendixson's Criterion
 * -------------------------------------------------------------------------*/

/**
 * @brief Apply Bendixson's criterion: check for absence of limit cycles (L4).
 *
 * Theorem (Bendixson's Criterion): For a 2D system dx/dt = f₁(x,y),
 * dy/dt = f₂(x,y), if the divergence ∂f₁/∂x + ∂f₂/∂y is not identically
 * zero and does not change sign in a simply connected region D, then
 * there are no closed orbits (limit cycles) entirely contained in D.
 *
 * This function computes the divergence on the grid and checks
 * for sign changes.
 *
 * Complexity: O(nx · ny · cost(f)).
 *
 * @param pp       Phase portrait (vector field already computed).
 * @param has_cycle Output: 0 if Bendixson guarantees NO cycles,
 *                  1 if criterion is inconclusive.
 * @return 0 on success.
 */
int nl_bendixson_criterion(const nl_phase_portrait_t *pp, int *has_cycle);

/**
 * @brief Apply Dulac's criterion (generalized Bendixson) (L4).
 *
 * Theorem (Dulac's Criterion): If there exists a C¹ function B(x,y) > 0
 * (a Dulac function) in a simply connected region D such that
 * ∂(B f₁)/∂x + ∂(B f₂)/∂y does not change sign in D, then there are
 * no closed orbits in D.
 *
 * Common Dulac functions: B = 1, B = 1/(xy), B = x^α y^β.
 *
 * @param sys      2D system.
 * @param pp       Phase portrait (for region definition).
 * @param dulac_func  Dulac function evaluator (can be NULL for B=1).
 * @param has_cycle Output: 0 if criterion guarantees NO cycles.
 * @return 0 on success.
 */
int nl_dulac_criterion(nl_nonlinear_system_t *sys,
                       const nl_phase_portrait_t *pp,
                       double (*dulac_func)(double x, double y),
                       int *has_cycle);

/* ---------------------------------------------------------------------------
 * L5, L6: Specialized Phase Analysis Functions
 * -------------------------------------------------------------------------*/

/**
 * @brief Classify all equilibria in a phase portrait (L5).
 *
 * Computes the Jacobian at each equilibrium and determines its
 * stability type (node, focus, saddle, center, etc.).
 *
 * @param sys         2D system.
 * @param pp          Phase portrait with equilibria.
 * @param stability   Output: stability classification per equilibrium.
 * @return Number classified.
 */
int nl_classify_equilibria_2d(nl_nonlinear_system_t *sys,
                               nl_phase_portrait_t *pp,
                               nl_stability_t *stability);

/**
 * @brief Compute the index of a closed curve for 2D vector field (L3, L5).
 *
 * The Poincaré index of an isolated equilibrium is defined as:
 *   I = (1/2π) ∮_Γ dθ
 *
 * where θ = atan2(f₂, f₁). For 2D systems:
 *   - Saddle: I = -1
 *   - Node, focus, center: I = +1
 *
 * The sum of indices of all equilibria inside a closed curve
 * determines the possible existence of limit cycles inside.
 *
 * @param sys       2D system.
 * @param center_x  x-coordinate of closed curve center.
 * @param center_y  y-coordinate.
 * @param radius    Radius of circular integration path.
 * @param n_pts     Number of integration points.
 * @param index     Output: Poincaré index (should be integer).
 * @return 0 on success.
 */
int nl_poincare_index(nl_nonlinear_system_t *sys,
                      double center_x, double center_y,
                      double radius, size_t n_pts,
                      double *index);

/**
 * @brief Trace the stable and unstable manifolds of a saddle point (L5, L8).
 *
 * For a 2D saddle, the eigenvectors of the Jacobian define the
 * tangent directions of the stable (W^s) and unstable (W^u) manifolds.
 *
 * This function integrates forward from points near the saddle along
 * the unstable eigenvector, and backward along the stable eigenvector,
 * to trace the global manifolds.
 *
 * @param sys      2D system.
 * @param saddle   Equilibrium state (must be a saddle).
 * @param T_fwd    Forward integration time.
 * @param T_bwd    Backward integration time.
 * @param dt       Step size.
 * @param W_u      Output: unstable manifold trajectory.
 * @param W_s      Output: stable manifold trajectory.
 * @return 0 on success.
 */
int nl_trace_manifolds(nl_nonlinear_system_t *sys,
                       const nl_state_t *saddle,
                       double T_fwd, double T_bwd, double dt,
                       nl_trajectory_t *W_u, nl_trajectory_t *W_s);

/**
 * @brief Construct a Poincaré map (first return map) for a 2D system (L6, L8).
 *
 * For a system with a periodic orbit, the Poincaré map P maps a point
 * on a cross-section Σ to its first return to Σ under the flow.
 *
 * A fixed point of the Poincaré map corresponds to a periodic orbit.
 * The eigenvalues (Floquet multipliers) of DP determine the stability
 * of the periodic orbit.
 *
 * This is the 2D analogue of the Poincaré section method.
 *
 * @param sys        2D system.
 * @param section_x  x-coordinate defining the section y = section_y.
 * @param x_range    Range of x values on the section.
 * @param n_pts      Number of initial points on the section.
 * @param T_max      Maximum integration time for each return.
 * @param dt         Integration step.
 * @param return_x   Output: first return x-coordinates.
 * @param return_t   Output: return times.
 * @return Number of successful returns, -1 on error.
 */
int nl_poincare_map(nl_nonlinear_system_t *sys,
                    double section_x, const double *x_range,
                    size_t n_pts, double T_max, double dt,
                    double *return_x, double *return_t);

/**
 * @brief Compute isoclines (curves of constant slope) in phase plane (L5).
 *
 * The dy/dx = f₂/f₁ isoclines are curves of constant slope.
 * Special cases:
 *   - Slope = 0     → y-nullcline (f₂ = 0)
 *   - Slope = ∞  → x-nullcline (f₁ = 0)
 *
 * This function computes isoclines for arbitrary slopes.
 *
 * @param sys    2D system.
 * @param slope  Target slope (can be INFINITY for vertical).
 * @param x_vals Array of x coordinates for evaluation.
 * @param n_pts  Number of points.
 * @param y_vals Output: y coordinate of isocline at each x.
 * @return 0 on success.
 */
int nl_compute_isocline(nl_nonlinear_system_t *sys,
                        double slope, const double *x_vals,
                        size_t n_pts, double *y_vals);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_PHASE_H */
