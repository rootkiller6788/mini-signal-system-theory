/**
 * @file nonlinear_applications.h
 * @brief Applications of nonlinear system theory (L7, L8).
 *
 * L7: GPS PLL, MEMS resonators, PA DPD, sigma-delta ADC, Lotka-Volterra
 * L8: Stochastic resonance, Lyapunov ROA estimation
 *
 * Course Mapping:
 *   - MIT 6.450 / Stanford EE359: GPS PLL
 *   - Berkeley EE 222: MEMS nonlinear dynamics
 *   - Tsinghua: 通信原理 (PA linearization), 系统辨识
 *
 * Reference:
 *   - Kaplan & Hegarty, "Understanding GPS/GNSS", 2017
 *   - Younis, "MEMS Linear and Nonlinear Statics and Dynamics", 2011
 *   - Gammaitoni et al., "Stochastic Resonance", Rev. Mod. Phys., 1998
 */

#ifndef NONLINEAR_APPLICATIONS_H
#define NONLINEAR_APPLICATIONS_H

#include "nonlinear_core.h"
#include "nonlinear_volterra.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L7: GPS Costas PLL Lock Detection
 * ========================================================================= */

/**
 * @brief Detect PLL lock from I/Q correlator outputs (L7).
 *
 * Lock indicator: L = |ΣI| / Σ√(I²+Q²). If L > threshold, locked.
 * Also returns RMS phase error for tracking jitter estimation.
 *
 * Reference: Kaplan & Hegarty, "Understanding GPS/GNSS", Ch.5.
 */
int nl_costas_pll_lock_detect(const double *i_samples,
                               const double *q_samples,
                               size_t N, double lock_threshold,
                               int *is_locked,
                               double *phase_error_rms);

/* =========================================================================
 * L7: MEMS Duffing Resonator
 * ========================================================================= */

/**
 * @brief Compute amplitude-frequency backbone curve for Duffing resonator.
 *
 * Solves: A²[(ω₀²-Ω²+(3α/4)A²)² + (γΩ)²] = F²
 * via fixed-point iteration. Hardening (α>0) bends right → hysteresis.
 *
 * Reference: Younis, "MEMS Linear and Nonlinear Statics and Dynamics", 2011.
 */
int nl_mems_duffing_backbone(double omega0, double alpha,
                              double gamma, double F,
                              const double *omega_range, size_t n_freqs,
                              double *amplitude, double *phase);

/* =========================================================================
 * L7: PA Digital Predistortion (Volterra DPD)
 * ========================================================================= */

/**
 * @brief Estimate Volterra DPD coefficients for PA linearization.
 *
 * Indirect learning architecture: post-distorter on (y/G)→u,
 * then use as predistorter. Memory polynomial model.
 *
 * Reference: Morgan et al., IEEE TSP, 2006.
 */
int nl_pa_volterra_dpd(const double complex *input,
                        const double complex *output,
                        size_t N, size_t mem_len, int order,
                        double gain, double *coeffs);

/* =========================================================================
 * L6, L7: FitzHugh-Nagumo Neuron
 * ========================================================================= */

/** FitzHugh-Nagumo parameters */
typedef struct {
    double a, b, eps, I_ext;
} nl_fhn_params_t;

/**
 * @brief Simulate FHN neuron: dv/dt=v-v³/3-w+I_ext, dw/dt=ε(v+a-bw).
 *
 * Detects spikes (v crosses +1.0 from below).
 * Exhibits: quiescence, tonic spiking, bistability.
 *
 * Reference: FitzHugh, Biophys. J., 1961.
 */
int nl_fitzhugh_nagumo(const nl_fhn_params_t *params,
                        const double *x0, double T, double dt,
                        nl_trajectory_t *traj, int *spike_count);

/* =========================================================================
 * L8: Stochastic Resonance
 * ========================================================================= */

/**
 * @brief Compute SNR vs noise for bistable system (SR phenomenon).
 *
 * Potential: U(x) = -a x²/2 + b x⁴/4.
 * Langevin: dx/dt = a x - b x³ + A cos(2πf₀t) + √(2D) ξ(t).
 *
 * Finds optimal D maximizing output SNR (Monte Carlo simulation).
 *
 * Reference: Gammaitoni et al., Rev. Mod. Phys., 1998.
 */
int nl_stochastic_resonance_snr(double a, double b, double A,
                                 double f0, const double *D_values,
                                 size_t nD, double T_sim, double dt,
                                 double *snr_db, double *D_opt);

/* =========================================================================
 * L7: Sigma-Delta Modulator Limit Cycle
 * ========================================================================= */

/**
 * @brief Detect idle tones / limit cycles in ΣΔ modulator.
 *
 * v[n] = v[n-1] + u - sgn(v[n-1]) for first-order.
 * Periodic output → idle tones → SNR degradation.
 */
int nl_sigma_delta_limit_cycle(int order, double input_dc,
                                size_t n_samples,
                                size_t *lc_period,
                                double *lc_amplitude,
                                double *tone_freq);

/* =========================================================================
 * L6, L7: Lotka-Volterra Predator-Prey
 * ========================================================================= */

typedef struct {
    double alpha, beta, gamma, delta;
} nl_lotka_volterra_params_t;

/**
 * @brief Simulate Lotka-Volterra: dx/dt=αx-βxy, dy/dt=δxy-γy.
 *
 * Equilibrium (γ/δ, α/β) is a center → neutral stability.
 * All trajectories (except equilibrium) are closed orbits.
 *
 * Applications: fishery management, epidemiology, chemical oscillations.
 * Reference: Lotka (1920), Volterra (1926).
 */
int nl_lotka_volterra(const nl_lotka_volterra_params_t *params,
                       const double *x0, double T, double dt,
                       nl_trajectory_t *traj,
                       double *avg_prey, double *avg_pred);

/* =========================================================================
 * L8: Region of Attraction Estimation (Lyapunov sublevel sets)
 * ========================================================================= */

/**
 * @brief Estimate ROA via largest Lyapunov sublevel set with V_dot<0.
 *
 * 1. Sample grid around equilibrium
 * 2. For each point, compute V and V_dot
 * 3. c_max = min{V(x) | V_dot(x) ≥ 0} — ROA inner estimate.
 *
 * Reference: Chesi, "Domain of Attraction", Springer, 2011.
 */
int nl_estimate_roa(nl_nonlinear_system_t *sys,
                     const nl_state_t *equilibrium,
                     const double *P, double grid_radius,
                     size_t n_grid, double *c_max);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_APPLICATIONS_H */
