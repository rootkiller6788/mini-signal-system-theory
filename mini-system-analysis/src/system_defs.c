/**
 * @file system_defs.c
 * @brief L1 Definitions Implementation - Signal and System Lifecycle
 *
 * Implements memory management and basic operations for all core
 * data structures defined in system_defs.h.
 *
 * Each function handles allocation, initialization, boundary checks,
 * and cleanup. All functions zero-initialize outputs and validate
 * inputs to prevent undefined behavior.
 *
 * Knowledge Implementation:
 *   L1-I01: ct_signal_alloc / ct_signal_free
 *   L1-I02: dt_signal_alloc / dt_signal_free
 *   L1-I03: dt_complex_signal_alloc / dt_complex_signal_free
 *   L1-I04: ct_tf_alloc / ct_tf_free
 *   L1-I05: dt_tf_alloc / dt_tf_free
 *   L1-I06: freq_resp_alloc / freq_resp_free
 *   L1-I07: pz_collection_alloc / pz_collection_free
 *   L1-I08: ct_signal_copy (deep copy)
 *   L1-I09: dt_signal_copy (deep copy)
 *   L1-I10: ct_signal_zero (fill with zeros)
 *   L1-I11: dt_signal_zero (fill with zeros)
 *   L1-I12: ct_signal_is_valid (boundary validation)
 *   L1-I13: dt_signal_is_valid
 *   L1-I14: Signal energy computation (Parseval-related)
 *   L1-I15: Signal power computation (average power)
 *   L1-I16: linear spaced time vector generation
 *   L1-I17: log spaced frequency vector generation
 *   L1-I18: dt_signal_from_ct (sampling / decimation)
 *   L1-I19: ct_signal_from_dt (reconstruction)
 *   L1-I20: ct_signal_append (concatenate two CT signals)
 */

#include "system_defs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ────────────────────────────────────────────────────────────────────────
 *  CT Signal Lifecycle (L1-I01)
 * ──────────────────────────────────────────────────────────────────────── */

ct_signal_t ct_signal_alloc(size_t num_samples, double t_start,
                            double t_end, double sample_rate)
{
    ct_signal_t sig;
    memset(&sig, 0, sizeof(sig));

    if (num_samples == 0 || num_samples > SIG_MAX_SAMPLES) {
        return sig; /* invalid - zeroed struct */
    }

    sig.samples = (double *)calloc(num_samples, sizeof(double));
    if (!sig.samples) {
        return sig;
    }

    sig.num_samples = num_samples;
    sig.t_start = t_start;
    sig.t_end = t_end;
    sig.sample_rate = sample_rate;
    sig.owns_buffer = 1;
    return sig;
}

void ct_signal_free(ct_signal_t *sig)
{
    if (!sig) return;
    if (sig->owns_buffer && sig->samples) {
        free(sig->samples);
        sig->samples = NULL;
    }
    sig->num_samples = 0;
    sig->owns_buffer = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  DT Signal Lifecycle (L1-I02)
 * ──────────────────────────────────────────────────────────────────────── */

dt_signal_t dt_signal_alloc(size_t num_samples)
{
    dt_signal_t sig;
    memset(&sig, 0, sizeof(sig));

    if (num_samples == 0 || num_samples > SIG_MAX_SAMPLES) {
        return sig;
    }

    sig.samples = (double *)calloc(num_samples, sizeof(double));
    if (!sig.samples) {
        return sig;
    }

    sig.num_samples = num_samples;
    sig.owns_buffer = 1;
    return sig;
}

void dt_signal_free(dt_signal_t *sig)
{
    if (!sig) return;
    if (sig->owns_buffer && sig->samples) {
        free(sig->samples);
        sig->samples = NULL;
    }
    sig->num_samples = 0;
    sig->owns_buffer = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Complex DT Signal Lifecycle (L1-I03)
 * ──────────────────────────────────────────────────────────────────────── */

dt_complex_signal_t dt_complex_signal_alloc(size_t num_samples)
{
    dt_complex_signal_t sig;
    memset(&sig, 0, sizeof(sig));

    if (num_samples == 0 || num_samples > SIG_MAX_SAMPLES) {
        return sig;
    }

    sig.samples = (double complex *)calloc(num_samples, sizeof(double complex));
    if (!sig.samples) {
        return sig;
    }

    sig.num_samples = num_samples;
    sig.owns_buffer = 1;
    return sig;
}

void dt_complex_signal_free(dt_complex_signal_t *sig)
{
    if (!sig) return;
    if (sig->owns_buffer && sig->samples) {
        free(sig->samples);
        sig->samples = NULL;
    }
    sig->num_samples = 0;
    sig->owns_buffer = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Transfer Function Lifecycle (L1-I04, L1-I05)
 * ──────────────────────────────────────────────────────────────────────── */

ct_transfer_function_t ct_tf_alloc(int num_order, int den_order)
{
    ct_transfer_function_t tf;
    memset(&tf, 0, sizeof(tf));

    if (num_order < 0 || den_order < 0 || den_order > 100) {
        return tf;
    }

    tf.num_coeffs = (double *)calloc((size_t)(num_order + 1), sizeof(double));
    tf.den_coeffs = (double *)calloc((size_t)(den_order + 1), sizeof(double));

    if (!tf.num_coeffs || !tf.den_coeffs) {
        free(tf.num_coeffs);
        free(tf.den_coeffs);
        memset(&tf, 0, sizeof(tf));
        return tf;
    }

    tf.num_order = num_order;
    tf.den_order = den_order;
    tf.owns_coeffs = 1;
    return tf;
}

void ct_tf_free(ct_transfer_function_t *tf)
{
    if (!tf) return;
    if (tf->owns_coeffs) {
        free(tf->num_coeffs);
        free(tf->den_coeffs);
    }
    tf->num_coeffs = NULL;
    tf->den_coeffs = NULL;
    tf->num_order = 0;
    tf->den_order = 0;
    tf->owns_coeffs = 0;
}

dt_transfer_function_t dt_tf_alloc(int num_order, int den_order)
{
    dt_transfer_function_t tf;
    memset(&tf, 0, sizeof(tf));

    if (num_order < 0 || den_order < 0 || den_order > 100) {
        return tf;
    }

    tf.num_coeffs = (double *)calloc((size_t)(num_order + 1), sizeof(double));
    tf.den_coeffs = (double *)calloc((size_t)(den_order + 1), sizeof(double));

    if (!tf.num_coeffs || !tf.den_coeffs) {
        free(tf.num_coeffs);
        free(tf.den_coeffs);
        memset(&tf, 0, sizeof(tf));
        return tf;
    }

    tf.num_order = num_order;
    tf.den_order = den_order;
    tf.owns_coeffs = 1;
    return tf;
}

void dt_tf_free(dt_transfer_function_t *tf)
{
    if (!tf) return;
    if (tf->owns_coeffs) {
        free(tf->num_coeffs);
        free(tf->den_coeffs);
    }
    tf->num_coeffs = NULL;
    tf->den_coeffs = NULL;
    tf->num_order = 0;
    tf->den_order = 0;
    tf->owns_coeffs = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Frequency Response Lifecycle (L1-I06)
 * ──────────────────────────────────────────────────────────────────────── */

frequency_response_t freq_resp_alloc(int num_points, double f_start,
                                      double f_end, int log_scale)
{
    frequency_response_t fr;
    memset(&fr, 0, sizeof(fr));

    if (num_points <= 0 || num_points > 100000) {
        return fr;
    }

    fr.points = (freq_point_t *)calloc((size_t)num_points, sizeof(freq_point_t));
    if (!fr.points) {
        return fr;
    }

    fr.num_points = num_points;
    fr.f_start = f_start;
    fr.f_end = f_end;
    fr.log_scale = log_scale;
    return fr;
}

void freq_resp_free(frequency_response_t *fr)
{
    if (!fr) return;
    free(fr->points);
    fr->points = NULL;
    fr->num_points = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Pole-Zero Collection Lifecycle (L1-I07)
 * ──────────────────────────────────────────────────────────────────────── */

pole_zero_collection_t pz_collection_alloc(int num_poles, int num_zeros)
{
    pole_zero_collection_t pz;
    memset(&pz, 0, sizeof(pz));

    if (num_poles < 0 || num_zeros < 0) return pz;

    if (num_poles > 0) {
        pz.poles = (pole_zero_t *)calloc((size_t)num_poles, sizeof(pole_zero_t));
        if (!pz.poles) { memset(&pz, 0, sizeof(pz)); return pz; }
    }
    if (num_zeros > 0) {
        pz.zeros = (pole_zero_t *)calloc((size_t)num_zeros, sizeof(pole_zero_t));
        if (!pz.zeros) {
            free(pz.poles);
            memset(&pz, 0, sizeof(pz));
            return pz;
        }
    }

    pz.num_poles = num_poles;
    pz.num_zeros = num_zeros;
    return pz;
}

void pz_collection_free(pole_zero_collection_t *pz)
{
    if (!pz) return;
    free(pz->poles);
    free(pz->zeros);
    pz->poles = NULL;
    pz->zeros = NULL;
    pz->num_poles = 0;
    pz->num_zeros = 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Copy Operations (L1-I08, L1-I09)
 * ──────────────────────────────────────────────────────────────────────── */

int ct_signal_copy(const ct_signal_t *src, ct_signal_t *dst)
{
    if (!src || !dst || !src->samples || src->num_samples == 0) {
        return -1;
    }

    ct_signal_t tmp = ct_signal_alloc(src->num_samples,
                                       src->t_start, src->t_end,
                                       src->sample_rate);
    if (!tmp.samples) return -1;

    memcpy(tmp.samples, src->samples, src->num_samples * sizeof(double));

    if (dst->owns_buffer && dst->samples) {
        free(dst->samples);
    }
    *dst = tmp;
    return 0;
}

int dt_signal_copy(const dt_signal_t *src, dt_signal_t *dst)
{
    if (!src || !dst || !src->samples || src->num_samples == 0) {
        return -1;
    }

    dt_signal_t tmp = dt_signal_alloc(src->num_samples);
    if (!tmp.samples) return -1;

    memcpy(tmp.samples, src->samples, src->num_samples * sizeof(double));

    if (dst->owns_buffer && dst->samples) {
        free(dst->samples);
    }
    *dst = tmp;
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Initialization (L1-I10, L1-I11)
 * ──────────────────────────────────────────────────────────────────────── */

void ct_signal_zero(ct_signal_t *sig)
{
    if (!sig || !sig->samples) return;
    memset(sig->samples, 0, sig->num_samples * sizeof(double));
}

void dt_signal_zero(dt_signal_t *sig)
{
    if (!sig || !sig->samples) return;
    memset(sig->samples, 0, sig->num_samples * sizeof(double));
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Validation (L1-I12, L1-I13)
 * ──────────────────────────────────────────────────────────────────────── */

int ct_signal_is_valid(const ct_signal_t *sig)
{
    if (!sig) return 0;
    if (!sig->samples) return 0;
    if (sig->num_samples == 0) return 0;
    if (sig->t_end <= sig->t_start) return 0;
    if (sig->sample_rate <= 0.0) return 0;

    /* Check for NaN/Inf in first few samples */
    for (size_t i = 0; i < sig->num_samples && i < 10; i++) {
        if (!isfinite(sig->samples[i])) return 0;
    }
    return 1;
}

int dt_signal_is_valid(const dt_signal_t *sig)
{
    if (!sig) return 0;
    if (!sig->samples) return 0;
    if (sig->num_samples == 0) return 0;

    for (size_t i = 0; i < sig->num_samples && i < 10; i++) {
        if (!isfinite(sig->samples[i])) return 0;
    }
    return 1;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Energy (L1-I14): E = Σ|x[n]|^2 (DT) or ∫|x(t)|^2 dt (CT)
 * ──────────────────────────────────────────────────────────────────────── */

double ct_signal_energy(const ct_signal_t *sig)
{
    if (!sig || !sig->samples || sig->num_samples < 2) return 0.0;

    double dt = (sig->t_end - sig->t_start) / (double)(sig->num_samples - 1);
    double energy = 0.0;

    /* Trapezoidal integration of |x(t)|^2 */
    for (size_t i = 0; i < sig->num_samples - 1; i++) {
        double avg_sq = 0.5 * (sig->samples[i] * sig->samples[i] +
                                sig->samples[i + 1] * sig->samples[i + 1]);
        energy += avg_sq * dt;
    }
    return energy;
}

double dt_signal_energy(const dt_signal_t *sig)
{
    if (!sig || !sig->samples) return 0.0;

    double energy = 0.0;
    for (size_t i = 0; i < sig->num_samples; i++) {
        energy += sig->samples[i] * sig->samples[i];
    }
    return energy;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Power (L1-I15): P = lim_{T->inf} (1/2T) ∫_{-T}^{T} |x(t)|^2 dt
 *  For finite signals: P = E / duration
 * ──────────────────────────────────────────────────────────────────────── */

double ct_signal_power(const ct_signal_t *sig)
{
    if (!sig || !sig->samples || sig->num_samples < 2) return 0.0;

    double energy = ct_signal_energy(sig);
    double duration = sig->t_end - sig->t_start;
    if (duration <= 0.0) return 0.0;
    return energy / duration;
}

double dt_signal_power(const dt_signal_t *sig)
{
    if (!sig || !sig->samples || sig->num_samples == 0) return 0.0;

    double energy = dt_signal_energy(sig);
    return energy / (double)sig->num_samples;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Time/Frequency Vector Generation (L1-I16, L1-I17)
 * ──────────────────────────────────────────────────────────────────────── */

int ct_signal_linear_time_vector(ct_signal_t *sig)
{
    if (!sig || !sig->samples || sig->num_samples < 2) return -1;

    double dt = (sig->t_end - sig->t_start) / (double)(sig->num_samples - 1);
    /* Note: this fills the signal samples with time values, which is
     * for the time axis, not the signal itself. Used internally. */
    for (size_t i = 0; i < sig->num_samples; i++) {
        sig->samples[i] = sig->t_start + (double)i * dt;
    }
    return 0;
}

int gen_log_spaced_frequencies(double f_start, double f_end, int num_points,
                                double *frequencies)
{
    if (!frequencies || num_points < 2 || f_start <= 0.0 || f_end <= f_start) {
        return -1;
    }

    double log_start = log10(f_start);
    double log_end = log10(f_end);
    double step = (log_end - log_start) / (double)(num_points - 1);

    for (int i = 0; i < num_points; i++) {
        frequencies[i] = pow(10.0, log_start + (double)i * step);
    }
    return 0;
}

int gen_linear_spaced_frequencies(double f_start, double f_end, int num_points,
                                   double *frequencies)
{
    if (!frequencies || num_points < 2 || f_end <= f_start) {
        return -1;
    }

    double step = (f_end - f_start) / (double)(num_points - 1);
    for (int i = 0; i < num_points; i++) {
        frequencies[i] = f_start + (double)i * step;
    }
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Domain Conversion (L1-I18, L1-I19)
 * ──────────────────────────────────────────────────────────────────────── */

int dt_signal_from_ct(const ct_signal_t *ct_sig, dt_signal_t *dt_sig,
                       double sample_interval)
{
    if (!ct_sig || !dt_sig || sample_interval <= 0.0) return -1;

    double ct_dt = (ct_sig->t_end - ct_sig->t_start) /
                   (double)(ct_sig->num_samples - 1);
    double T = sample_interval;

    /* Number of DT samples: floor((t_end - t_start) / T) + 1 */
    size_t dt_len = (size_t)((ct_sig->t_end - ct_sig->t_start) / T) + 1;
    if (dt_len < 1 || dt_len > SIG_MAX_SAMPLES) return -1;

    dt_signal_t tmp = dt_signal_alloc(dt_len);
    if (!tmp.samples) return -1;

    for (size_t n = 0; n < dt_len; n++) {
        double t = ct_sig->t_start + (double)n * T;
        /* Linear interpolation: find nearest CT sample */
        double idx_f = (t - ct_sig->t_start) / ct_dt;
        size_t idx = (size_t)idx_f;
        if (idx >= ct_sig->num_samples - 1) {
            tmp.samples[n] = ct_sig->samples[ct_sig->num_samples - 1];
        } else {
            double frac = idx_f - (double)idx;
            tmp.samples[n] = ct_sig->samples[idx] * (1.0 - frac) +
                             ct_sig->samples[idx + 1] * frac;
        }
    }

    if (dt_sig->owns_buffer && dt_sig->samples) free(dt_sig->samples);
    *dt_sig = tmp;
    return 0;
}

int ct_signal_from_dt(const dt_signal_t *dt_sig, ct_signal_t *ct_sig,
                       double sample_rate)
{
    if (!dt_sig || !ct_sig || sample_rate <= 0.0) return -1;

    double T = 1.0 / sample_rate;
    double duration = (double)(dt_sig->num_samples - 1) * T;

    ct_signal_t tmp = ct_signal_alloc(dt_sig->num_samples, 0.0,
                                       duration, sample_rate);
    if (!tmp.samples) return -1;

    /* Zero-order hold reconstruction (staircase) */
    memcpy(tmp.samples, dt_sig->samples,
           dt_sig->num_samples * sizeof(double));

    if (ct_sig->owns_buffer && ct_sig->samples) free(ct_sig->samples);
    *ct_sig = tmp;
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Signal Concatenation (L1-I20)
 * ──────────────────────────────────────────────────────────────────────── */

int ct_signal_append(const ct_signal_t *a, const ct_signal_t *b,
                      ct_signal_t *result)
{
    if (!a || !b || !result) return -1;
    if (!a->samples || !b->samples) return -1;
    if (a->sample_rate != b->sample_rate) return -1;

    size_t total = a->num_samples + b->num_samples;
    /* double gap = b->t_start - a->t_end; -- unused, kept for documentation */
    double total_duration = b->t_end - a->t_start;

    ct_signal_t tmp = ct_signal_alloc(total, a->t_start,
                                       a->t_start + total_duration,
                                       a->sample_rate);
    if (!tmp.samples) return -1;

    memcpy(tmp.samples, a->samples, a->num_samples * sizeof(double));
    memcpy(tmp.samples + a->num_samples, b->samples,
           b->num_samples * sizeof(double));

    if (result->owns_buffer && result->samples) free(result->samples);
    *result = tmp;
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Unit Step / Unit Impulse Generators
 * ──────────────────────────────────────────────────────────────────────── */

int dt_signal_set_unit_impulse(dt_signal_t *sig, size_t impulse_index)
{
    if (!sig || !sig->samples) return -1;
    if (impulse_index >= sig->num_samples) return -1;

    memset(sig->samples, 0, sig->num_samples * sizeof(double));
    sig->samples[impulse_index] = 1.0;
    return 0;
}

int dt_signal_set_unit_step(dt_signal_t *sig, size_t step_index)
{
    if (!sig || !sig->samples) return -1;
    if (step_index >= sig->num_samples) return -1;

    for (size_t i = 0; i < sig->num_samples; i++) {
        sig->samples[i] = (i >= step_index) ? 1.0 : 0.0;
    }
    return 0;
}

int ct_signal_set_rectangular_pulse(ct_signal_t *sig,
                                     double start_time, double width,
                                     double amplitude)
{
    if (!sig || !sig->samples) return -1;

    for (size_t i = 0; i < sig->num_samples; i++) {
        double t = sig->t_start + (double)i *
                   (sig->t_end - sig->t_start) / (double)(sig->num_samples - 1);
        if (t >= start_time && t < start_time + width) {
            sig->samples[i] = amplitude;
        } else {
            sig->samples[i] = 0.0;
        }
    }
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 *  Utility: Norm and Comparison
 * ──────────────────────────────────────────────────────────────────────── */

double dt_signal_max_abs_diff(const dt_signal_t *a, const dt_signal_t *b)
{
    if (!a || !b || !a->samples || !b->samples) return -1.0;
    if (a->num_samples != b->num_samples) return -1.0;

    double max_diff = 0.0;
    for (size_t i = 0; i < a->num_samples; i++) {
        double diff = fabs(a->samples[i] - b->samples[i]);
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

double dt_signal_rms(const dt_signal_t *sig)
{
    if (!sig || !sig->samples || sig->num_samples == 0) return 0.0;
    return sqrt(dt_signal_power(sig));
}

int dt_signal_is_approx_equal(const dt_signal_t *a, const dt_signal_t *b,
                               double tolerance)
{
    return dt_signal_max_abs_diff(a, b) <= tolerance;
}
