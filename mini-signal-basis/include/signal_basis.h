/*
 * signal_basis.h - Core Signal Type Definitions (L1: Definitions)
 *
 * Fundamental building blocks of signal theory:
 *   - Continuous-time and Discrete-time signal data structures
 *   - Elementary signals: unit step, Dirac delta, ramp, exponential,
 *     sinusoid, complex exponential, rectangular/triangular/Gaussian pulse, sinc
 *   - Signal parameters: amplitude, frequency, phase, period, energy, power
 *   - Signal classifications: periodic/aperiodic, causal, even/odd, energy/power
 *
 * Course: MIT 6.003 Ch.1 / Stanford EE102A Ch.1-2 / Berkeley EE123 Ch.1
 *         ETH 227-0427 Ch.1 / TU Munich / Tsinghua Ch.1
 * Textbook: Oppenheim & Willsky (1997) Ch.1; Proakis & Manolakis (2007) Ch.1
 */
#ifndef SIGNAL_BASIS_H
#define SIGNAL_BASIS_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* L1: Fundamental type aliases */
typedef double signal_value_t;
typedef double signal_time_t;
typedef double signal_freq_t;
typedef double signal_phase_t;
typedef int64_t signal_index_t;

/* L1: Signal classifications (enums) */
typedef enum { SIGNAL_DOMAIN_CONTINUOUS = 0, SIGNAL_DOMAIN_DISCRETE = 1 } signal_domain_t;
typedef enum { SIGNAL_PERIODICITY_APERIODIC = 0, SIGNAL_PERIODICITY_PERIODIC = 1 } signal_periodicity_t;
typedef enum { SIGNAL_CAUSALITY_NONCAUSAL = 0, SIGNAL_CAUSALITY_CAUSAL = 1, SIGNAL_CAUSALITY_ANTICAUSAL = 2 } signal_causality_t;
typedef enum { SIGNAL_CLASS_UNKNOWN = 0, SIGNAL_CLASS_ENERGY = 1, SIGNAL_CLASS_POWER = 2, SIGNAL_CLASS_NEITHER = 3 } signal_energy_class_t;
typedef enum { SIGNAL_TYPE_DETERMINISTIC = 0, SIGNAL_TYPE_STOCHASTIC = 1 } signal_randomness_t;
typedef enum { SIGNAL_SYMMETRY_NONE = 0, SIGNAL_SYMMETRY_EVEN = 1, SIGNAL_SYMMETRY_ODD = 2 } signal_symmetry_t;

/* L1: Continuous-Time Signal - Primary Data Structure
 * x(t) on [t_start, t_end], uniformly sampled: x[k] = x(t_start + k*dt) */
typedef struct {
    signal_value_t *samples;
    signal_time_t   t_start;
    signal_time_t   t_end;
    signal_time_t   dt;
    signal_index_t  num_samples;
    signal_domain_t domain;
} signal_ct_t;

/* L1: Discrete-Time Signal - Primary Data Structure
 * x[n] on [n_start, n_end] */
typedef struct {
    signal_value_t *samples;
    signal_index_t  n_start;
    signal_index_t  n_end;
    signal_index_t  length;
    signal_domain_t domain;
} signal_dt_t;

/* L1: Elementary Signal Parameter Structs */
typedef struct { signal_value_t amplitude; signal_freq_t frequency; signal_phase_t phase; } signal_sinusoid_params_t;
typedef struct { signal_value_t amplitude; signal_value_t decay_rate; } signal_exponential_params_t;
typedef struct { signal_value_t amplitude; signal_freq_t frequency; signal_phase_t phase; } signal_cisoid_params_t;
typedef struct { signal_value_t amplitude; signal_time_t center; signal_time_t width; } signal_rect_params_t;
typedef struct { signal_value_t amplitude; signal_time_t center; signal_time_t width; } signal_tri_params_t;
typedef struct { signal_value_t amplitude; signal_time_t center; signal_value_t sigma; } signal_gaussian_params_t;
typedef struct { signal_value_t amplitude; signal_time_t center; signal_time_t zero_crossing; } signal_sinc_params_t;

/* L1: Signal Property Functions (O&W section 1.2) */
double signal_ct_energy(const signal_ct_t *s);
double signal_ct_power(const signal_ct_t *s);
double signal_dt_energy(const signal_dt_t *s);
double signal_dt_power(const signal_dt_t *s);
signal_energy_class_t signal_ct_classify(const signal_ct_t *s);
double signal_ct_rms(const signal_ct_t *s);
double signal_dt_rms(const signal_dt_t *s);
double signal_ct_peak(const signal_ct_t *s);
signal_periodicity_t signal_ct_detect_periodicity(const signal_ct_t *s, double threshold);
double signal_ct_fundamental_period(const signal_ct_t *s);
signal_causality_t signal_ct_check_causality(const signal_ct_t *s, double epsilon);
signal_symmetry_t signal_ct_check_symmetry(const signal_ct_t *s, double epsilon);

/* L1: Memory Management */
signal_ct_t *signal_ct_alloc(signal_time_t t_start, signal_time_t t_end, signal_time_t dt);
signal_dt_t *signal_dt_alloc(signal_index_t n_start, signal_index_t n_end);
void signal_ct_free(signal_ct_t *s);
void signal_dt_free(signal_dt_t *s);
signal_ct_t *signal_ct_copy(const signal_ct_t *src);
signal_dt_t *signal_dt_copy(const signal_dt_t *src);

/* L1: Elementary Signal Generators - each implements one independent definition from O&W Ch.1 */
void signal_ct_fill_unit_step(signal_ct_t *s, signal_time_t step_time);
void signal_ct_fill_dirac_delta(signal_ct_t *s, signal_time_t t0);
void signal_ct_fill_unit_ramp(signal_ct_t *s);
void signal_ct_fill_sinusoid(signal_ct_t *s, const signal_sinusoid_params_t *params);
void signal_ct_fill_cisoid(signal_ct_t *s, const signal_cisoid_params_t *params);
void signal_ct_fill_exponential(signal_ct_t *s, const signal_exponential_params_t *params);
void signal_ct_fill_rect_pulse(signal_ct_t *s, const signal_rect_params_t *params);
void signal_ct_fill_tri_pulse(signal_ct_t *s, const signal_tri_params_t *params);
void signal_ct_fill_gaussian(signal_ct_t *s, const signal_gaussian_params_t *params);
void signal_ct_fill_sinc(signal_ct_t *s, const signal_sinc_params_t *params);

/* L2: Inner Product Space Operations */
double sinc_normalized(double x);
double signal_ct_inner_product(const signal_ct_t *a, const signal_ct_t *b);
double signal_ct_norm(const signal_ct_t *a);
double signal_ct_mean(const signal_ct_t *s);
double signal_ct_variance(const signal_ct_t *s);


/* ===================================================================
 * L2: Additional Signal Operations
 * =================================================================== */

/** Numerical derivative: dx/dt via central difference (boundary: forward/backward) */
int signal_ct_derivative(const signal_ct_t *src, signal_ct_t *dst);

/** Numerical integral: cumulative trapezoidal rule */
int signal_ct_integral(const signal_ct_t *src, signal_ct_t *dst);

/** Linear interpolation to a new sampling rate */
int signal_ct_interpolate(const signal_ct_t *src, signal_ct_t *dst);

/** Decimation (downsampling) by an integer factor */
int signal_ct_decimate(const signal_ct_t *src, int factor, signal_ct_t *dst);

/** Zero-crossing rate (sign changes per second) */
double signal_ct_zero_crossing_rate(const signal_ct_t *s);

/** Moving-window peak envelope detection */
int signal_ct_envelope(const signal_ct_t *s, signal_ct_t *envelope, int window_size);

/** Signal-to-Noise Ratio between original and noisy signal (in dB) */
double signal_ct_snr(const signal_ct_t *clean, const signal_ct_t *noisy);
#endif /* SIGNAL_BASIS_H */
