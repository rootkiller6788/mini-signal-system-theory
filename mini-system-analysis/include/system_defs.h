/**
 * @file system_defs.h
 * @brief L1 Definitions — Core System Classification and Signal Types
 *
 * This header defines the fundamental data structures used throughout
 * system analysis. Each struct maps to a core definition from signal
 * and system theory (Oppenheim & Willsky, 1997).
 *
 * Knowledge Coverage:
 *   L1-D01: Continuous-time signal representation
 *   L1-D02: Discrete-time signal representation with sampling info
 *   L1-D03: System type classification enum (LTI, LTV, NL, etc.)
 *   L1-D04: System property flags (causal, stable, memoryless, invertible)
 *   L1-D05: Impulse response — continuous-time h(t)
 *   L1-D06: Impulse response — discrete-time h[n]
 *   L1-D07: Step response — continuous-time s(t)
 *   L1-D08: Step response — discrete-time s[n]
 *   L1-D09: Transfer function representation (rational in s)
 *   L1-D10: Transfer function representation (rational in z)
 *   L1-D11: Pole/Zero list for system characterization
 *   L1-D12: Frequency response magnitude/phase pair
 *   L1-D13: System order and relative degree
 *
 * Course Mapping:
 *   MIT 6.003 Ch.1-2; Stanford EE102A Ch.1-3;
 *   Berkeley EE123 Ch.1-2; 清华 信号与系统 Ch.1-2
 */

#ifndef SYSTEM_DEFS_H
#define SYSTEM_DEFS_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Signal Domain Enumeration ────────────────────────────────────────── */

typedef enum {
    DOMAIN_CONTINUOUS = 0,
    DOMAIN_DISCRETE   = 1,
    DOMAIN_DIGITAL    = 2,
    DOMAIN_MIXED      = 3
} signal_domain_t;

/* ── System Type Classification ───────────────────────────────────────── */

typedef enum {
    SYSTYPE_LTI     = 0,
    SYSTYPE_LTV     = 1,
    SYSTYPE_NL_TI   = 2,
    SYSTYPE_NL_TV   = 3,
    SYSTYPE_HYBRID  = 4,
    SYSTYPE_STOCHASTIC = 5
} system_type_t;

/* ── System Property Bitmask ──────────────────────────────────────────── */

typedef enum {
    PROP_LINEAR         = (1 << 0),
    PROP_TIME_INVARIANT = (1 << 1),
    PROP_CAUSAL         = (1 << 2),
    PROP_MEMORYLESS     = (1 << 3),
    PROP_INVERTIBLE     = (1 << 4),
    PROP_BIBO_STABLE    = (1 << 5),
    PROP_FIR            = (1 << 6),
    PROP_IIR            = (1 << 7),
    PROP_MINIMUM_PHASE  = (1 << 8),
    PROP_PASSIVE        = (1 << 9),
    PROP_LOSSLESS       = (1 << 10),
    PROP_LUMPED         = (1 << 11),
    PROP_PROPER         = (1 << 12),
    PROP_STRICTLY_PROPER= (1 << 13)
} system_property_t;

/* ── Continuous-Time Signal ───────────────────────────────────────────── */

#define SIG_MAX_SAMPLES 65536

typedef struct {
    double  *samples;
    size_t   num_samples;
    double   t_start;
    double   t_end;
    double   sample_rate;
    int      owns_buffer;
} ct_signal_t;

/* ── Discrete-Time Signal ─────────────────────────────────────────────── */

typedef struct {
    double  *samples;
    size_t   num_samples;
    int      owns_buffer;
} dt_signal_t;

/* ── Complex Signal Extension ─────────────────────────────────────────── */

typedef struct {
    double complex *samples;
    size_t          num_samples;
    int             owns_buffer;
} dt_complex_signal_t;

/* ── Impulse Response (CT) ────────────────────────────────────────────── */

typedef struct {
    ct_signal_t  response;
    double       energy;
    double       integral;
    double       peak_time;
    double       peak_value;
    double       settling_time;
} ct_impulse_response_t;

/* ── Impulse Response (DT) ────────────────────────────────────────────── */

typedef struct {
    dt_signal_t  response;
    double       energy;
    double       abs_sum;
    int          length;
} dt_impulse_response_t;

/* ── Step Response ────────────────────────────────────────────────────── */

typedef struct {
    ct_signal_t  response;
    double       steady_state;
    double       rise_time;
    double       overshoot_pct;
    double       delay_time;
    double       peak_time;
    int          is_oscillatory;
} ct_step_response_t;

typedef struct {
    dt_signal_t  response;
    double       steady_state;
    int          deadbeat;
    int          settle_sample;
} dt_step_response_t;

/* ── Transfer Function (CT, Laplace Domain) ───────────────────────────── */

typedef struct {
    double  *num_coeffs;
    double  *den_coeffs;
    int      num_order;
    int      den_order;
    int      owns_coeffs;
} ct_transfer_function_t;

/* ── Transfer Function (DT, Z Domain) ─────────────────────────────────── */

typedef struct {
    double  *num_coeffs;
    double  *den_coeffs;
    int      num_order;
    int      den_order;
    int      owns_coeffs;
} dt_transfer_function_t;

/* ── Pole and Zero Identification ─────────────────────────────────────── */

typedef struct {
    double complex value;
    int            is_real;
    int            is_pole;
    double         damping;
    double         wn;
} pole_zero_t;

typedef struct {
    pole_zero_t *poles;
    int          num_poles;
    pole_zero_t *zeros;
    int          num_zeros;
    double       dc_gain;
    int          is_stable;
} pole_zero_collection_t;

/* ── Frequency Response ───────────────────────────────────────────────── */

typedef struct {
    double  frequency;
    double  magnitude;
    double  magnitude_db;
    double  phase_rad;
    double  phase_deg;
    double  group_delay;
} freq_point_t;

typedef struct {
    freq_point_t *points;
    int            num_points;
    double         f_start;
    double         f_end;
    int            log_scale;
    double         bandwidth;
    double         resonance_freq;
    double         resonance_peak_db;
} frequency_response_t;

/* ── System Order and Classification ──────────────────────────────────── */

typedef struct {
    int     order;
    int     type_number;
    int     relative_degree;
    int     is_minimum_phase;
    int     is_allpass;
    double  dc_gain;
    double  infinite_gain;
} system_order_info_t;

/* ── Lifecycle Functions ──────────────────────────────────────────────── */

ct_signal_t  ct_signal_alloc(size_t num_samples, double t_start,
                             double t_end, double sample_rate);
void         ct_signal_free(ct_signal_t *sig);

dt_signal_t  dt_signal_alloc(size_t num_samples);
void         dt_signal_free(dt_signal_t *sig);

dt_complex_signal_t dt_complex_signal_alloc(size_t num_samples);
void                dt_complex_signal_free(dt_complex_signal_t *sig);

ct_transfer_function_t ct_tf_alloc(int num_order, int den_order);
void                   ct_tf_free(ct_transfer_function_t *tf);

dt_transfer_function_t dt_tf_alloc(int num_order, int den_order);
void                   dt_tf_free(dt_transfer_function_t *tf);

frequency_response_t freq_resp_alloc(int num_points, double f_start,
                                      double f_end, int log_scale);
void                 freq_resp_free(frequency_response_t *fr);

pole_zero_collection_t pz_collection_alloc(int num_poles, int num_zeros);
void                   pz_collection_free(pole_zero_collection_t *pz);

/* ── Signal Operations (extended) ──────────────────────────────────────── */

double ct_signal_energy(const ct_signal_t *sig);
double dt_signal_energy(const dt_signal_t *sig);
double ct_signal_power(const ct_signal_t *sig);
double dt_signal_power(const dt_signal_t *sig);

int ct_signal_is_valid(const ct_signal_t *sig);
int dt_signal_is_valid(const dt_signal_t *sig);

void ct_signal_zero(ct_signal_t *sig);
void dt_signal_zero(dt_signal_t *sig);

int ct_signal_copy(const ct_signal_t *src, ct_signal_t *dst);
int dt_signal_copy(const dt_signal_t *src, dt_signal_t *dst);

int dt_signal_set_unit_impulse(dt_signal_t *sig, size_t impulse_index);
int dt_signal_set_unit_step(dt_signal_t *sig, size_t step_index);
int ct_signal_set_rectangular_pulse(ct_signal_t *sig, double start_time, double width, double amplitude);

int dt_signal_from_ct(const ct_signal_t *ct_sig, dt_signal_t *dt_sig, double sample_interval);
int ct_signal_from_dt(const dt_signal_t *dt_sig, ct_signal_t *ct_sig, double sample_rate);

int ct_signal_append(const ct_signal_t *a, const ct_signal_t *b, ct_signal_t *result);

double dt_signal_max_abs_diff(const dt_signal_t *a, const dt_signal_t *b);
double dt_signal_rms(const dt_signal_t *sig);
int dt_signal_is_approx_equal(const dt_signal_t *a, const dt_signal_t *b, double tolerance);

int gen_log_spaced_frequencies(double f_start, double f_end, int num_points, double *frequencies);
int gen_linear_spaced_frequencies(double f_start, double f_end, int num_points, double *frequencies);

#ifdef __cplusplus
}
#endif

#endif /* SYSTEM_DEFS_H */
