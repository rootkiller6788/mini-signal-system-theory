#include "signal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

signal_ct_t *signal_ct_alloc(signal_time_t t_start, signal_time_t t_end, signal_time_t dt)
{
    if (dt <= 0.0 || t_end <= t_start) return NULL;
    signal_index_t n = (signal_index_t)((t_end - t_start) / dt) + 1;
    if (n <= 0 || n > 100000000) return NULL;
    signal_ct_t *s = (signal_ct_t *)malloc(sizeof(signal_ct_t));
    if (!s) return NULL;
    s->samples = (signal_value_t *)calloc((size_t)n, sizeof(signal_value_t));
    if (!s->samples) { free(s); return NULL; }
    s->t_start = t_start; s->t_end = t_end; s->dt = dt;
    s->num_samples = n; s->domain = SIGNAL_DOMAIN_CONTINUOUS;
    return s;
}
signal_dt_t *signal_dt_alloc(signal_index_t n_start, signal_index_t n_end)
{
    if (n_end < n_start) return NULL;
    signal_index_t len = n_end - n_start + 1;
    if (len <= 0 || len > 100000000) return NULL;
    signal_dt_t *s = (signal_dt_t *)malloc(sizeof(signal_dt_t));
    if (!s) return NULL;
    s->samples = (signal_value_t *)calloc((size_t)len, sizeof(signal_value_t));
    if (!s->samples) { free(s); return NULL; }
    s->n_start = n_start;
    s->n_end = n_end;
    s->length = len;
    s->domain = SIGNAL_DOMAIN_DISCRETE;
    return s;
}

void signal_ct_free(signal_ct_t *s)
{
    if (s) { free(s->samples); free(s); }
}

void signal_dt_free(signal_dt_t *s)
{
    if (s) { free(s->samples); free(s); }
}

signal_ct_t *signal_ct_copy(const signal_ct_t *src)
{
    if (!src) return NULL;
    signal_ct_t *d = signal_ct_alloc(src->t_start, src->t_end, src->dt);
    if (!d) return NULL;
    memcpy(d->samples, src->samples, (size_t)src->num_samples * sizeof(signal_value_t));
    return d;
}

signal_dt_t *signal_dt_copy(const signal_dt_t *src)
{
    if (!src) return NULL;
    signal_dt_t *d = signal_dt_alloc(src->n_start, src->n_end);
    if (!d) return NULL;
    memcpy(d->samples, src->samples, (size_t)src->length * sizeof(signal_value_t));
    return d;
}





signal_symmetry_t signal_ct_check_symmetry(const signal_ct_t *s, double epsilon)
{
    if (!s || s->num_samples < 2) return SIGNAL_SYMMETRY_NONE;
    signal_index_t N = s->num_samples;
    double max_even_diff = 0.0, max_odd_sum = 0.0;
    for (signal_index_t k = 0; k < N / 2; k++) {
        double v_plus = s->samples[k];
        double v_minus = s->samples[N - 1 - k];
        double even_diff = fabs(v_plus - v_minus);
        double odd_sum = fabs(v_plus + v_minus);
        if (even_diff > max_even_diff) max_even_diff = even_diff;
        if (odd_sum > max_odd_sum) max_odd_sum = odd_sum;
    }
    if (max_odd_sum < epsilon) return SIGNAL_SYMMETRY_ODD;
    if (max_even_diff < epsilon) return SIGNAL_SYMMETRY_EVEN;
    return SIGNAL_SYMMETRY_NONE;
}

/* =================================================================
 * L1: Signal Property Functions
 * Reference: Oppenheim & Willsky (1997) section 1.2
 * ================================================================= */

double signal_ct_energy(const signal_ct_t *s)
{
    if (!s || !s->samples) return 0.0;
    double E = 0.0;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double v = s->samples[k];
        E += v * v;
    }
    return E * s->dt;
}

double signal_ct_power(const signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples == 0) return 0.0;
    double E = signal_ct_energy(s);
    double duration = s->t_end - s->t_start;
    if (duration <= 0.0) return 0.0;
    return E / duration;
}

double signal_dt_energy(const signal_dt_t *s)
{
    if (!s || !s->samples) return 0.0;
    double E = 0.0;
    for (signal_index_t i = 0; i < s->length; i++) {
        double v = s->samples[i];
        E += v * v;
    }
    return E;
}

double signal_dt_power(const signal_dt_t *s)
{
    if (!s || !s->samples || s->length == 0) return 0.0;
    double E = signal_dt_energy(s);
    return E / (double)s->length;
}

signal_energy_class_t signal_ct_classify(const signal_ct_t *s)
{
    if (!s) return SIGNAL_CLASS_UNKNOWN;
    double E = signal_ct_energy(s);
    double P = signal_ct_power(s);
    double EPS = 1e-10;
    if (E < 1e300 && P < EPS) return SIGNAL_CLASS_ENERGY;
    if (P > EPS && P < 1e300) return SIGNAL_CLASS_POWER;
    return SIGNAL_CLASS_NEITHER;
}

double signal_ct_rms(const signal_ct_t *s)
{
    double P = signal_ct_power(s);
    return sqrt(P);
}

double signal_dt_rms(const signal_dt_t *s)
{
    double P = signal_dt_power(s);
    return sqrt(P);
}

double signal_ct_peak(const signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples == 0) return 0.0;
    double peak = fabs(s->samples[0]);
    for (signal_index_t k = 1; k < s->num_samples; k++) {
        double absv = fabs(s->samples[k]);
        if (absv > peak) peak = absv;
    }
    return peak;
}

/* Periodicity detection via autocorrelation peak search.
 * Computes normalized autocorrelation, finds maximum at non-zero lag.
 * If max correlation > threshold, classifies as periodic. */
signal_periodicity_t signal_ct_detect_periodicity(const signal_ct_t *s, double threshold)
{
    if (!s || s->num_samples < 4) return SIGNAL_PERIODICITY_APERIODIC;
    signal_index_t N = s->num_samples;
    double r0 = 0.0;
    for (signal_index_t k = 0; k < N; k++)
        r0 += s->samples[k] * s->samples[k];
    if (r0 < 1e-15) return SIGNAL_PERIODICITY_APERIODIC;
    double max_corr = 0.0;
    signal_index_t max_lag_idx = (N > 500) ? (N / 2) : (N - 1);
    for (signal_index_t lag = 1; lag < max_lag_idx; lag++) {
        double sum = 0.0; int count = 0;
        for (signal_index_t k = 0; k < N - lag; k++) {
            sum += s->samples[k] * s->samples[k + lag];
            count++;
        }
        if (count > 0) {
            double corr = sum / r0;
            if (corr > max_corr) max_corr = corr;
        }
    }
    return (max_corr > threshold) ? SIGNAL_PERIODICITY_PERIODIC : SIGNAL_PERIODICITY_APERIODIC;
}

/* Estimate fundamental period via autocorrelation first peak */
double signal_ct_fundamental_period(const signal_ct_t *s)
{
    if (!s || s->num_samples < 4) return 0.0;
    signal_index_t N = s->num_samples;
    double r0 = 0.0;
    for (signal_index_t k = 0; k < N; k++) r0 += s->samples[k] * s->samples[k];
    if (r0 < 1e-15) return 0.0;
    signal_index_t best_lag = 0;
    double best_corr = -1.0;
    signal_index_t max_lag_idx = (N > 500) ? (N / 2) : (N - 1);
    for (signal_index_t lag = 1; lag < max_lag_idx; lag++) {
        double sum = 0.0; int count = 0;
        for (signal_index_t k = 0; k < N - lag; k++) {
            sum += s->samples[k] * s->samples[k + lag]; count++;
        }
        if (count > 0) { double corr = sum / r0;
            if (corr > best_corr) { best_corr = corr; best_lag = lag; } }
    }
    if (best_lag == 0) return 0.0;
    return (double)best_lag * s->dt;
}

signal_causality_t signal_ct_check_causality(const signal_ct_t *s, double epsilon)
{
    if (!s) return SIGNAL_CAUSALITY_NONCAUSAL;
    int any_before = 0, any_after = 0;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        if (fabs(s->samples[k]) > epsilon) {
            if (t < 0.0) any_before = 1;
            else any_after = 1;
        }
    }
    if (any_before && !any_after) return SIGNAL_CAUSALITY_ANTICAUSAL;
    if (!any_before && any_after) return SIGNAL_CAUSALITY_CAUSAL;
    return SIGNAL_CAUSALITY_NONCAUSAL;
}


/* =================================================================
 * L1: Elementary Signal Generators
 * Each fills a pre-allocated signal with samples of one canonical signal.
 * Reference: Oppenheim & Willsky (1997) section 1.4
 * ================================================================= */

void signal_ct_fill_unit_step(signal_ct_t *s, signal_time_t step_time)
{
    if (!s) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        s->samples[k] = (t >= step_time) ? 1.0 : 0.0;
    }
}

void signal_ct_fill_dirac_delta(signal_ct_t *s, signal_time_t t0)
{
    if (!s) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        if (fabs(t - t0) < s->dt * 0.5)
            s->samples[k] = 1.0 / s->dt;
        else
            s->samples[k] = 0.0;
    }
}

void signal_ct_fill_unit_ramp(signal_ct_t *s)
{
    if (!s) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        s->samples[k] = (t >= 0.0) ? t : 0.0;
    }
}

void signal_ct_fill_sinusoid(signal_ct_t *s, const signal_sinusoid_params_t *params)
{
    if (!s || !params) return;
    double omega = 2.0 * M_PI * params->frequency;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        s->samples[k] = params->amplitude * cos(omega * t + params->phase);
    }
}

void signal_ct_fill_cisoid(signal_ct_t *s, const signal_cisoid_params_t *params)
{
    if (!s || !params) return;
    double omega = 2.0 * M_PI * params->frequency;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        double angle = omega * t + params->phase;
        s->samples[k] = params->amplitude * cos(angle);
    }
}

void signal_ct_fill_exponential(signal_ct_t *s, const signal_exponential_params_t *params)
{
    if (!s || !params) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        s->samples[k] = params->amplitude * exp(params->decay_rate * t);
    }
}

void signal_ct_fill_rect_pulse(signal_ct_t *s, const signal_rect_params_t *params)
{
    if (!s || !params) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        double x = (t - params->center) / params->width;
        if (fabs(x) < 0.5) s->samples[k] = params->amplitude;
        else if (fabs(x) < 0.5 + 1e-12) s->samples[k] = params->amplitude * 0.5;
        else s->samples[k] = 0.0;
    }
}

void signal_ct_fill_tri_pulse(signal_ct_t *s, const signal_tri_params_t *params)
{
    if (!s || !params) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        double x = fabs((t - params->center) / params->width);
        s->samples[k] = (x < 1.0) ? params->amplitude * (1.0 - x) : 0.0;
    }
}

void signal_ct_fill_gaussian(signal_ct_t *s, const signal_gaussian_params_t *params)
{
    if (!s || !params) return;
    double sigma2 = 2.0 * params->sigma * params->sigma;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        double dt0 = t - params->center;
        s->samples[k] = params->amplitude * exp(-(dt0 * dt0) / sigma2);
    }
}

void signal_ct_fill_sinc(signal_ct_t *s, const signal_sinc_params_t *params)
{
    if (!s || !params) return;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        double x = (t - params->center) / params->zero_crossing;
        s->samples[k] = params->amplitude * sinc_normalized(x);
    }
}

/* =================================================================
 * L2: Inner Product Space Operations
 * Reference: Oppenheim & Willsky (1997) section 2.2
 * ================================================================= */

double sinc_normalized(double x)
{
    if (fabs(x) < 1e-8) {
        double pix = M_PI * x;
        double pix2 = pix * pix;
        return 1.0 - pix2 / 6.0 + pix2 * pix2 / 120.0;
    }
    return sin(M_PI * x) / (M_PI * x);
}

double signal_ct_inner_product(const signal_ct_t *a, const signal_ct_t *b)
{
    if (!a || !b || !a->samples || !b->samples) return 0.0;
    signal_index_t n = (a->num_samples < b->num_samples) ? a->num_samples : b->num_samples;
    double ip = 0.0;
    for (signal_index_t k = 0; k < n; k++)
        ip += a->samples[k] * b->samples[k];
    return ip * a->dt;
}

double signal_ct_norm(const signal_ct_t *a)
{
    double ip = signal_ct_inner_product(a, a);
    return sqrt(ip);
}

double signal_ct_mean(const signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples == 0) return 0.0;
    double sum = 0.0;
    for (signal_index_t k = 0; k < s->num_samples; k++)
        sum += s->samples[k];
    return sum / (double)s->num_samples;
}

double signal_ct_variance(const signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples == 0) return 0.0;
    double mu = signal_ct_mean(s);
    double var = 0.0;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double diff = s->samples[k] - mu;
        var += diff * diff;
    }
    return var / (double)s->num_samples;
}

/* =================================================================
 * L2: Signal Differentiation and Integration
 * Numerical derivative: dx/dt ~= (x[k+1] - x[k]) / dt
 * Integral: cumulative trapezoidal rule
 * These are fundamental operations linking signals to systems.
 * ================================================================= */

int signal_ct_derivative(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    if (n < 2) return -1;
    double inv_dt = 1.0 / src->dt;
    for (signal_index_t k = 0; k < n; k++) {
        if (k == 0)
            dst->samples[k] = (src->samples[1] - src->samples[0]) * inv_dt;
        else if (k == n - 1)
            dst->samples[k] = (src->samples[k] - src->samples[k-1]) * inv_dt;
        else
            dst->samples[k] = (src->samples[k+1] - src->samples[k-1]) * 0.5 * inv_dt;
    }
    return 0;
}

int signal_ct_integral(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    signal_index_t n = (src->num_samples < dst->num_samples) ? src->num_samples : dst->num_samples;
    if (n < 2) return -1;
    dst->samples[0] = 0.0;
    for (signal_index_t k = 1; k < n; k++)
        dst->samples[k] = dst->samples[k-1] + (src->samples[k-1] + src->samples[k]) * 0.5 * src->dt;
    return 0;
}

/* =================================================================
 * L2: Signal Interpolation (Linear)
 * Resample a signal to a new sampling rate using linear interpolation.
 * ================================================================= */

int signal_ct_interpolate(const signal_ct_t *src, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        double t_dst = dst->t_start + (double)k * dst->dt;
        double idx_src = (t_dst - src->t_start) / src->dt;
        signal_index_t i0 = (signal_index_t)idx_src;
        if (i0 < 0) dst->samples[k] = src->samples[0];
        else if (i0 >= src->num_samples - 1) dst->samples[k] = src->samples[src->num_samples - 1];
        else {
            double frac = idx_src - (double)i0;
            dst->samples[k] = src->samples[i0] * (1.0 - frac) + src->samples[i0 + 1] * frac;
        }
    }
    return 0;
}

/* =================================================================
 * L2: Signal Decimation (Downsampling by integer factor)
 * Keep every M-th sample. Anti-aliasing not applied (simplified).
 * ================================================================= */

int signal_ct_decimate(const signal_ct_t *src, int factor, signal_ct_t *dst)
{
    if (!src || !dst || !src->samples || !dst->samples || factor <= 1) return -1;
    for (signal_index_t k = 0; k < dst->num_samples; k++) {
        signal_index_t src_idx = k * (signal_index_t)factor;
        if (src_idx < src->num_samples)
            dst->samples[k] = src->samples[src_idx];
        else
            dst->samples[k] = 0.0;
    }
    return 0;
}

/* =================================================================
 * L2: Zero-Crossing Rate
 * Number of sign changes per second. Used in speech/audio analysis.
 * ZCR is a simple feature for voice activity detection.
 * ================================================================= */

double signal_ct_zero_crossing_rate(const signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples < 2) return 0.0;
    int crossings = 0;
    for (signal_index_t k = 1; k < s->num_samples; k++) {
        if (s->samples[k-1] * s->samples[k] < 0.0) crossings++;
    }
    double duration = s->t_end - s->t_start;
    return (duration > 0.0) ? (double)crossings / duration : 0.0;
}

/* =================================================================
 * L2: Signal Envelope Detection via Hilbert transform approximation
 * Simple peak-envelope: track local maxima.
 * ================================================================= */

int signal_ct_envelope(const signal_ct_t *s, signal_ct_t *envelope, int window_size)
{
    if (!s || !envelope || !s->samples || !envelope->samples || window_size < 1) return -1;
    signal_index_t n = (s->num_samples < envelope->num_samples) ? s->num_samples : envelope->num_samples;
    int half_win = window_size / 2;
    for (signal_index_t k = 0; k < n; k++) {
        signal_index_t start = (k > (signal_index_t)half_win) ? k - half_win : 0;
        signal_index_t end = (k + half_win < n) ? k + half_win : n - 1;
        double max_val = fabs(s->samples[start]);
        for (signal_index_t j = start + 1; j <= end; j++)
            if (fabs(s->samples[j]) > max_val) max_val = fabs(s->samples[j]);
        envelope->samples[k] = max_val;
    }
    return 0;
}

/* L2: Signal-to-Noise Ratio in dB
 * SNR = 10 * log10(P_signal / P_noise)
 * where noise = noisy - clean */
double signal_ct_snr(const signal_ct_t *clean, const signal_ct_t *noisy)
{
    if (!clean || !noisy || !clean->samples || !noisy->samples) return -999.0;
    signal_index_t n = (clean->num_samples < noisy->num_samples) ? clean->num_samples : noisy->num_samples;
    double P_sig = 0.0, P_noise = 0.0;
    for (signal_index_t k = 0; k < n; k++) {
        P_sig += clean->samples[k] * clean->samples[k];
        double e = noisy->samples[k] - clean->samples[k];
        P_noise += e * e;
    }
    if (P_noise < 1e-15) return 999.0;
    return 10.0 * log10(P_sig / P_noise);
}


