/**
 * @file nonlinear_applications.c
 * @brief Application-level nonlinear system analysis (L7, L8).
 */
#include "nonlinear_applications.h"
#include "nonlinear_stability.h"
#include "nonlinear_phase.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
/* Note: complex.h is included via nonlinear_volterra.h -> nonlinear_applications.h
   Avoid local variable name I which conflicts with complex.h I macro */

/* ===================================================================
 * GPS Costas PLL Lock Detection
 * L = |mean(I)| / mean(sqrt(I^2+Q^2)) in [0, 1]
 * sigma_phi = atan(RMS(Q) / RMS(I))
 * App: Boeing 787, Tesla GPS, iPhone GNSS.
 * =================================================================== */

int nl_costas_pll_lock_detect(const double *i_samples,
                               const double *q_samples,
                               size_t N, double lock_threshold,
                               int *is_locked,
                               double *phase_error_rms)
{
    if (!i_samples || !q_samples || !is_locked
        || !phase_error_rms || N < 10) {
        errno = EINVAL; return -1;
    }
    double sum_I = 0.0, sum_mag = 0.0;
    double sum_Q2 = 0.0, sum_I2 = 0.0;
    for (size_t k = 0; k < N; k++) {
        double ival = i_samples[k], qval = q_samples[k];
        sum_I += ival;
        sum_mag += sqrt(ival * ival + qval * qval);
        sum_I2 += ival * ival;
        sum_Q2 += qval * qval;
    }
    double lock_indicator = (sum_mag > 1e-12)
                            ? fabs(sum_I) / sum_mag : 0.0;
    *is_locked = (lock_indicator > lock_threshold) ? 1 : 0;
    double rms_I = sqrt(sum_I2 / (double)N);
    double rms_Q = sqrt(sum_Q2 / (double)N);
    *phase_error_rms = (rms_I > 1e-12)
                       ? atan2(rms_Q, rms_I) : M_PI / 2.0;
    return 0;
}

/* ===================================================================
 * MEMS Duffing Resonator Backbone Curve
 * A^2 = F^2 / [(w0^2 - Omg^2 + 0.75*alpha*A^2)^2 + (gamma*Omg)^2]
 * Fixed-point iteration with continuation.
 * App: MEMS gyroscopes, oscillators, inertial sensors.
 * =================================================================== */

int nl_mems_duffing_backbone(double omega0, double alpha,
                              double gamma, double F,
                              const double *omega_range, size_t n_freqs,
                              double *amplitude, double *phase)
{
    if (!omega_range || !amplitude || !phase || n_freqs < 2
        || omega0 <= 0.0 || gamma <= 0.0 || F <= 0.0) {
        errno = EINVAL; return -1;
    }
    double A = 0.0;
    for (size_t k = 0; k < n_freqs; k++) {
        double Omega = omega_range[0]
            + (omega_range[1] - omega_range[0]) * (double)k
              / (double)(n_freqs - 1);
        double A_prev = (A > 1e-6) ? A
                        : (F / (omega0 * omega0 + 1.0));
        int converged = 0;
        for (int iter = 0; iter < 200; iter++) {
            double detuning = omega0 * omega0 - Omega * Omega
                            + 0.75 * alpha * A_prev * A_prev;
            double denom = sqrt(detuning * detuning
                                + (gamma * Omega) * (gamma * Omega));
            if (denom < 1e-15) { A_prev = F / 1e-12; break; }
            double A_new = F / denom;
            if (fabs(A_new - A_prev) < 1e-8 * fmax(A_new, 1.0)) {
                A_prev = A_new; converged = 1; break;
            }
            A_prev = A_new;
        }
        A = (converged || A_prev > 0.0) ? A_prev : 0.0;
        amplitude[k] = A;
        double detuning = omega0 * omega0 - Omega * Omega
                        + 0.75 * alpha * A * A;
        phase[k] = atan2(-gamma * Omega, detuning);
    }
    return 0;
}


/* ===================================================================
 * PA Volterra DPD (Indirect Learning Architecture)
 * Memory polynomial: y[n] = sum c_{k,q} x[n-q] |x[n-q]|^{k-1}
 * DPD via least squares on (y/G) -> u mapping.
 * App: 5G NR base stations, 4G LTE eNodeB, satellite transponders.
 * =================================================================== */

int nl_pa_volterra_dpd(const double complex *input,
                        const double complex *output,
                        size_t N, size_t mem_len, int order,
                        double gain, double *coeffs)
{
    if (!input || !output || !coeffs || N <= mem_len
        || mem_len < 1 || order < 1 || order > 7 || gain <= 0.0) {
        errno = EINVAL; return -1;
    }
    int n_odd = (order + 1) / 2;
    double complex *z = (double complex *)malloc(
        N * sizeof(double complex));
    if (!z) { errno = ENOMEM; return -1; }
    for (size_t n = 0; n < N; n++)
        z[n] = output[n] / gain;
    for (size_t q = 0; q < mem_len; q++) {
        size_t start = (q > mem_len) ? q : mem_len;
        double *A_mat = (double *)calloc(
            (size_t)(n_odd * n_odd), sizeof(double));
        double *b_vec = (double *)calloc((size_t)n_odd, sizeof(double));
        if (!A_mat || !b_vec) {
            free(A_mat); free(b_vec); free(z); return -1;
        }
        for (size_t n = start; n < N; n++) {
            double x_mag = cabs(z[n - q]);
            double y_mag = cabs(input[n]);
            if (x_mag < 1e-12) continue;
            double basis[4] = {0};
            double xk = x_mag;
            for (int k_idx = 0; k_idx < n_odd; k_idx++) {
                basis[k_idx] = xk;
                xk *= x_mag * x_mag;
            }
            for (int ki = 0; ki < n_odd; ki++) {
                for (int kj = 0; kj < n_odd; kj++)
                    A_mat[ki * n_odd + kj] += basis[ki] * basis[kj];
                b_vec[ki] += basis[ki] * y_mag;
            }
        }
        for (int k = 0; k < n_odd; k++) {
            int max_r = k;
            double max_v = fabs(A_mat[k * n_odd + k]);
            for (int r = k + 1; r < n_odd; r++) {
                double av = fabs(A_mat[r * n_odd + k]);
                if (av > max_v) { max_v = av; max_r = r; }
            }
            if (max_v < 1e-15) continue;
            if (max_r != k) {
                for (int c = 0; c < n_odd; c++) {
                    double t = A_mat[k * n_odd + c];
                    A_mat[k * n_odd + c] = A_mat[max_r * n_odd + c];
                    A_mat[max_r * n_odd + c] = t;
                }
                double t = b_vec[k]; b_vec[k] = b_vec[max_r];
                b_vec[max_r] = t;
            }
            for (int r = k + 1; r < n_odd; r++) {
                double f = A_mat[r * n_odd + k]
                          / A_mat[k * n_odd + k];
                for (int c = k; c < n_odd; c++)
                    A_mat[r * n_odd + c] -= f * A_mat[k * n_odd + c];
                b_vec[r] -= f * b_vec[k];
            }
        }
        for (int k = n_odd - 1; k >= 0; k--) {
            double s = b_vec[k];
            for (int c = k + 1; c < n_odd; c++)
                s -= A_mat[k * n_odd + c] * b_vec[c];
            b_vec[k] = (fabs(A_mat[k * n_odd + k]) > 1e-15)
                      ? s / A_mat[k * n_odd + k] : 0.0;
        }
        for (int k_idx = 0; k_idx < n_odd; k_idx++)
            coeffs[k_idx * mem_len + q] = b_vec[k_idx];
        free(A_mat); free(b_vec);
    }
    free(z);
    return 0;
}


/* ===================================================================
 * FitzHugh-Nagumo Neuron Model
 * dv/dt = v - v^3/3 - w + I_ext
 * dw/dt = eps * (v + a - b*w)
 * v: fast membrane potential, w: slow recovery
 * Spike: v crosses +1.0 from below with positive slope.
 * App: neural prosthetics, cardiac pacemaker modeling.
 * =================================================================== */

static int fhn_ode(const double *x, double t, const double *params,
                    double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 2) return -1;
    double a = params[0], b = params[1], eps = params[2];
    double I_ext = params[3];
    double v = x[0], w = x[1];
    dxdt[0] = v - v * v * v / 3.0 - w + I_ext;
    dxdt[1] = eps * (v + a - b * w);
    return 0;
}

int nl_fitzhugh_nagumo(const nl_fhn_params_t *params,
                        const double *x0, double T, double dt,
                        nl_trajectory_t *traj, int *spike_count)
{
    if (!params || !x0 || !traj || !spike_count
        || T <= 0.0 || dt <= 0.0 || params->eps <= 0.0) {
        errno = EINVAL; return -1;
    }
    nl_nonlinear_system_t sys;
    double sys_params[4] = {params->a, params->b,
                            params->eps, params->I_ext};
    sys.dim = 2; sys.f = fhn_ode; sys.jacobian = NULL;
    memcpy(sys.params, sys_params, 4 * sizeof(double));
    sys.num_params = 4; sys.is_autonomous = 1;
    size_t n_steps;
    int ret = nl_integrate_rk4(&sys, x0, T, dt, traj, &n_steps);
    if (ret != 0) return ret;
    *spike_count = 0;
    int was_below = (traj->points[0].x[0] < 1.0) ? 1 : 0;
    for (size_t i = 1; i < traj->num_points; i++) {
        double v_prev = traj->points[i - 1].x[0];
        double v_curr = traj->points[i].x[0];
        if (was_below && v_curr >= 1.0 && v_curr > v_prev) {
            (*spike_count)++;
            was_below = 0;
        } else if (v_curr < 1.0) {
            was_below = 1;
        }
    }
    return 0;
}


/* ===================================================================
 * Stochastic Resonance SNR Estimation
 * Bistable potential: U(x) = -a*x^2/2 + b*x^4/4
 * Barrier height Delta_U = a^2/(4b).
 * Langevin: dx/dt = a*x - b*x^3 + A*cos(2*pi*f0*t) + sqrt(2*D)*xi(t)
 * Euler-Maruyama integration with Box-Muller N(0,1) generation.
 * SNR detection via synchronous demodulation at f0.
 * App: weak signal detection (F-35 radar, MRI, sonar).
 * =================================================================== */

/* Simple LCG: X_{n+1} = (1103515245 * X_n + 12345) mod 2^31 */
static unsigned int lcg_rand(unsigned int *seed)
{
    *seed = (1103515245u * (*seed) + 12345u) & 0x7fffffffu;
    return *seed;
}

static double rand_normal(unsigned int *seed)
{
    double u1 = (double)lcg_rand(seed) / (double)0x7fffffffu;
    double u2 = (double)lcg_rand(seed) / (double)0x7fffffffu;
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int nl_stochastic_resonance_snr(double a, double b, double A,
                                 double f0, const double *D_values,
                                 size_t nD, double T_sim, double dt,
                                 double *snr_db, double *D_opt)
{
    if (!D_values || !snr_db || !D_opt || nD < 2
        || a <= 0.0 || b <= 0.0 || T_sim <= 0.0 || dt <= 0.0) {
        errno = EINVAL; return -1;
    }
    double best_snr = -INFINITY;
    *D_opt = D_values[0];
    unsigned int seed = 12345u;
    for (size_t id = 0; id < nD; id++) {
        double D = D_values[id];
        size_t n_steps = (size_t)(T_sim / dt);
        if (n_steps < 100) n_steps = 100;
        double x = 0.0;
        double sqrt_2Ddt = sqrt(2.0 * D * dt);
        double omega = 2.0 * M_PI * f0;
        double cos_sum = 0.0, sin_sum = 0.0, sum_x2 = 0.0;
        size_t n_steady = n_steps / 4;
        for (size_t n = 0; n < n_steps; n++) {
            double t = (double)n * dt;
            double drift = a * x - b * x * x * x
                          + A * cos(omega * t);
            x += dt * drift + sqrt_2Ddt * rand_normal(&seed);
            if (n >= n_steady) {
                cos_sum += x * cos(omega * t);
                sin_sum += x * sin(omega * t);
                sum_x2 += x * x;
            }
        }
        size_t n_eff = n_steps - n_steady;
        if (n_eff < 10) { snr_db[id] = -100.0; continue; }
        double sig_power = 2.0 * (cos_sum * cos_sum
                                   + sin_sum * sin_sum)
                           / ((double)n_eff * (double)n_eff);
        double noise_power = sum_x2 / (double)n_eff - sig_power;
        if (noise_power < 1e-15) noise_power = 1e-15;
        double snr = sig_power / noise_power;
        snr_db[id] = 10.0 * log10(snr);
        if (snr_db[id] > best_snr) {
            best_snr = snr_db[id];
            *D_opt = D;
        }
    }
    return 0;
}


/* ===================================================================
 * Sigma-Delta Modulator Limit Cycle Detection
 * First-order SDM: v[n] = v[n-1] + u - sgn(v[n-1])
 * Idle tones at rational DC inputs degrade in-band SNR.
 * Detection via autocorrelation of output bit sequence.
 * App: audio SD ADC (iPhone, Tesla), MEMS microphones.
 * =================================================================== */

int nl_sigma_delta_limit_cycle(int order, double input_dc,
                                size_t n_samples,
                                size_t *lc_period,
                                double *lc_amplitude,
                                double *tone_freq)
{
    (void)order; /* higher-order SDM reserved for extension */
    if (!lc_period || !lc_amplitude || !tone_freq
        || n_samples < 100 || fabs(input_dc) >= 1.0) {
        errno = EINVAL; return -1;
    }
    int *output_bits = (int *)malloc(n_samples * sizeof(int));
    if (!output_bits) { errno = ENOMEM; return -1; }
    double v = 0.0;
    for (size_t n = 0; n < n_samples; n++) {
        int q_out = (v >= 0.0) ? 1 : -1;
        output_bits[n] = q_out;
        v = v + input_dc - (double)q_out;
    }
    size_t max_lag = n_samples / 2;
    double *acf = (double *)calloc(max_lag, sizeof(double));
    if (!acf) { free(output_bits); return -1; }
    for (size_t lag = 1; lag < max_lag; lag++) {
        double sum = 0.0;
        for (size_t n = 0; n < n_samples - lag; n++)
            sum += (double)output_bits[n] * (double)output_bits[n + lag];
        acf[lag] = sum / (double)(n_samples - lag);
    }
    *lc_period = 0;
    double best_acf = 0.0;
    for (size_t lag = 2; lag < max_lag; lag++) {
        if (acf[lag] > acf[lag - 1] && acf[lag] > acf[lag + 1]
            && acf[lag] > 0.5) {
            if (acf[lag] > best_acf) {
                best_acf = acf[lag];
                *lc_period = lag;
            }
        }
    }
    if (*lc_period > 0) {
        *tone_freq = 1.0 / (double)(*lc_period);
        double amp_sum = 0.0;
        for (size_t n = 0; n < *lc_period && n < n_samples; n++)
            amp_sum += (double)(output_bits[n] * output_bits[n]);
        *lc_amplitude = sqrt(amp_sum / (double)(*lc_period));
    } else {
        *tone_freq = 0.0;
        *lc_amplitude = 0.0;
    }
    free(acf);
    free(output_bits);
    return 0;
}


/* ===================================================================
 * Lotka-Volterra Predator-Prey Dynamics
 * dx/dt = alpha*x - beta*x*y     (prey)
 * dy/dt = delta*x*y - gamma*y    (predator)
 * Equilibrium (gamma/delta, alpha/beta) is a center (neutral stability).
 * First integral: V = delta*x - gamma*ln(x) + beta*y - alpha*ln(y)
 * Volterra's law: time averages equal equilibrium values.
 * App: fishery management, epidemiology, Easter Island model.
 * =================================================================== */

static int lv_ode(const double *x, double t, const double *params,
                   double *dxdt, size_t dim)
{
    (void)t;
    if (dim != 2) return -1;
    double alpha = params[0], beta = params[1];
    double gamma = params[2], delta = params[3];
    dxdt[0] = alpha * x[0] - beta * x[0] * x[1];
    dxdt[1] = delta * x[0] * x[1] - gamma * x[1];
    return 0;
}

int nl_lotka_volterra(const nl_lotka_volterra_params_t *params,
                       const double *x0, double T, double dt,
                       nl_trajectory_t *traj,
                       double *avg_prey, double *avg_pred)
{
    if (!params || !x0 || !traj || !avg_prey || !avg_pred
        || T <= 0.0 || dt <= 0.0) {
        errno = EINVAL; return -1;
    }
    nl_nonlinear_system_t sys;
    double sys_params[4] = {params->alpha, params->beta,
                            params->gamma, params->delta};
    sys.dim = 2; sys.f = lv_ode; sys.jacobian = NULL;
    memcpy(sys.params, sys_params, 4 * sizeof(double));
    sys.num_params = 4; sys.is_autonomous = 1;
    size_t n_steps;
    int ret = nl_integrate_rk4(&sys, x0, T, dt, traj, &n_steps);
    if (ret != 0) return ret;
    double sum_prey = 0.0, sum_pred = 0.0;
    for (size_t i = 0; i < traj->num_points; i++) {
        sum_prey += traj->points[i].x[0];
        sum_pred += traj->points[i].x[1];
    }
    *avg_prey = sum_prey / (double)traj->num_points;
    *avg_pred = sum_pred / (double)traj->num_points;
    return 0;
}


/* ===================================================================
 * Region of Attraction Estimation via Lyapunov Sublevel Sets
 * V(x) = (x-x*)^T P (x-x*)
 * c_max = min{V(x_i) : V_dot(x_i) >= -eps} over grid sampling.
 * For n <= 3: uniform grid. For higher dim: random shell sampling.
 * App: safe flight envelope (F-35), chemical reactor stability,
 *      power grid transient stability.
 * =================================================================== */

int nl_estimate_roa(nl_nonlinear_system_t *sys,
                     const nl_state_t *equilibrium,
                     const double *P, double grid_radius,
                     size_t n_grid, double *c_max)
{
    if (!sys || !equilibrium || !P || !c_max
        || grid_radius <= 0.0 || n_grid < 3
        || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }
    size_t n = sys->dim;
    *c_max = INFINITY;
    if (n <= 3) {
        size_t total_pts = 1;
        for (size_t d = 0; d < n; d++) total_pts *= n_grid;
        if (total_pts > 50000) total_pts = 50000;
        double dx = 2.0 * grid_radius / (double)(n_grid - 1);
        double *x_pt = (double *)malloc(n * sizeof(double));
        if (!x_pt) { errno = ENOMEM; return -1; }
        for (size_t idx = 0; idx < total_pts; idx++) {
            size_t rem = idx;
            for (size_t d = 0; d < n; d++) {
                size_t id = rem % n_grid; rem /= n_grid;
                x_pt[d] = equilibrium->x[d]
                         - grid_radius + (double)id * dx;
            }
            double V_val = 0.0;
            for (size_t i = 0; i < n; i++) {
                double di = x_pt[i] - equilibrium->x[i];
                double row_sum = 0.0;
                for (size_t j = 0; j < n; j++)
                    row_sum += P[i * n + j]
                              * (x_pt[j] - equilibrium->x[j]);
                V_val += di * row_sum;
            }
            double V_dot = nl_quadratic_lyapunov_dot(P, x_pt, sys, n);
            if (V_dot >= -1e-8 && V_val > 1e-6) {
                if (V_val < *c_max) *c_max = V_val;
            }
        }
        free(x_pt);
    } else {
        unsigned int seed = 54321u;
        size_t n_samples = 10000;
        double *x_pt = (double *)malloc(n * sizeof(double));
        if (!x_pt) { errno = ENOMEM; return -1; }
        for (size_t s = 0; s < n_samples; s++) {
            double r_sq = 0.0;
            for (size_t d = 0; d < n; d++) {
                x_pt[d] = rand_normal(&seed);
                r_sq += x_pt[d] * x_pt[d];
            }
            double r = sqrt(r_sq);
            if (r < 1e-12) continue;
            for (int shell = 1; shell <= 5; shell++) {
                double radius = grid_radius * (double)shell / 5.0;
                for (size_t d = 0; d < n; d++)
                    x_pt[d] = equilibrium->x[d]
                             + radius * x_pt[d] / r;
                double V_val = 0.0;
                for (size_t i = 0; i < n; i++) {
                    double di = x_pt[i] - equilibrium->x[i];
                    double row_sum = 0.0;
                    for (size_t j = 0; j < n; j++)
                        row_sum += P[i * n + j]
                                  * (x_pt[j] - equilibrium->x[j]);
                    V_val += di * row_sum;
                }
                double V_dot = nl_quadratic_lyapunov_dot(
                    P, x_pt, sys, n);
                if (V_dot >= -1e-8 && V_val > 1e-6) {
                    if (V_val < *c_max) *c_max = V_val;
                }
            }
        }
        free(x_pt);
    }
    if (!isfinite(*c_max)) *c_max = 1.0;
    if (*c_max < 1e-8) *c_max = 1e-8;
    return 0;
}
