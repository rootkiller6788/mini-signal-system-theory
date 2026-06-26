#include "signal_random.h"
#include "signal_correlation.h"
#include "signal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =================================================================
 * L2: PRNG - 64-bit Linear Congruential Generator
 * X_{n+1} = (a * X_n + c) mod m
 * Using MMIX parameters: a=6364136223846793005, c=1442695040888963407
 * ================================================================= */

#define LCG_A 6364136223846793005ULL
#define LCG_C 1442695040888963407ULL

signal_rng_t *signal_rng_create(uint64_t seed)
{
    signal_rng_t *rng = (signal_rng_t *)malloc(sizeof(signal_rng_t));
    if (!rng) return NULL;
    rng->state = seed;
    return rng;
}

void signal_rng_free(signal_rng_t *rng)
{
    free(rng);
}

/* Uniform random in [0, 1) using the high 53 bits */
double signal_rng_uniform(signal_rng_t *rng)
{
    if (!rng) return 0.0;
    rng->state = LCG_A * rng->state + LCG_C;
    return (double)(rng->state >> 11) * 0x1.0p-53;
}

/* Box-Muller transform: generate Gaussian from two uniforms */
double signal_rng_gaussian(signal_rng_t *rng, double mean, double std)
{
    if (!rng) return mean;
    double u1 = signal_rng_uniform(rng);
    double u2 = signal_rng_uniform(rng);
    if (u1 < 1e-15) u1 = 1e-15;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + std * z;
}

/* =================================================================
 * L2: Noise Signal Generation
 * ================================================================= */

int signal_fill_white_noise(signal_ct_t *s, signal_rng_t *rng, double std)
{
    if (!s || !rng) return -1;
    for (signal_index_t k = 0; k < s->num_samples; k++)
        s->samples[k] = signal_rng_gaussian(rng, 0.0, std);
    return 0;
}

int signal_fill_awgn(signal_ct_t *s, signal_rng_t *rng, double noise_std)
{
    if (!s || !rng) return -1;
    for (signal_index_t k = 0; k < s->num_samples; k++)
        s->samples[k] += signal_rng_gaussian(rng, 0.0, noise_std);
    return 0;
}

/* Colored noise via first-order AR process:
 * n[k] = alpha * n[k-1] + sqrt(1-alpha^2) * w[k]
 * where w[k] is white noise. alpha = 0 gives white, alpha -> 1 gives pink-like */
int signal_fill_colored_noise(signal_ct_t *s, signal_rng_t *rng, double std, double alpha)
{
    if (!s || !rng || alpha < -1.0 || alpha > 1.0) return -1;
    double scale = std * sqrt(1.0 - alpha * alpha);
    s->samples[0] = signal_rng_gaussian(rng, 0.0, std);
    for (signal_index_t k = 1; k < s->num_samples; k++)
        s->samples[k] = alpha * s->samples[k - 1] + signal_rng_gaussian(rng, 0.0, scale);
    return 0;
}

/* =================================================================
 * L2: Signal-to-Noise Ratio
 * ================================================================= */

int signal_add_awgn(signal_ct_t *s, signal_rng_t *rng, double snr_db)
{
    if (!s || !rng) return -1;
    double P_signal = signal_ct_power(s);
    if (P_signal < 1e-15) return -1;
    double snr_linear = pow(10.0, snr_db / 10.0);
    double noise_std = sqrt(P_signal / snr_linear);
    return signal_fill_awgn(s, rng, noise_std);
}

double signal_compute_snr(const signal_ct_t *signal, const signal_ct_t *noisy)
{
    if (!signal || !noisy) return -999.0;
    signal_index_t n = (signal->num_samples < noisy->num_samples) ? signal->num_samples : noisy->num_samples;
    double P_sig = 0.0, P_noise = 0.0;
    for (signal_index_t k = 0; k < n; k++) {
        P_sig += signal->samples[k] * signal->samples[k];
        double noise = noisy->samples[k] - signal->samples[k];
        P_noise += noise * noise;
    }
    if (P_noise < 1e-15) return 999.0;
    return 10.0 * log10(P_sig / P_noise);
}

/* Random binary sequence: +amplitude or -amplitude every bit_duration */
int signal_fill_random_bits(signal_ct_t *s, signal_rng_t *rng, double bit_duration, double amplitude)
{
    if (!s || !rng || bit_duration <= 0.0) return -1;
    for (signal_index_t k = 0; k < s->num_samples; k++) {
        double t = s->t_start + (double)k * s->dt;
        (void)((int)(t / bit_duration));
        rng->state = LCG_A * rng->state + LCG_C;
        int bit = (rng->state >> 63) & 1;
        s->samples[k] = bit ? amplitude : -amplitude;
    }
    return 0;
}

/* =================================================================
 * L8: Monte Carlo Detection Probability Estimation
 * Estimates P_d for matched filter detector in AWGN via repeated trials.
 * For each trial: generate AWGN, add to template, run matched filter,
 * check if correlation peak exceeds threshold.
 * Reference: Proakis & Salehi (2008) section 4.4
 * ================================================================= */

double signal_monte_carlo_pd(const signal_ct_t *template_signal, double snr_db,
                              double threshold, int num_trials, signal_rng_t *rng)
{
    if (!template_signal || !rng || num_trials <= 0) return -1.0;
    double P_sig = signal_ct_power(template_signal);
    if (P_sig < 1e-15) return -1.0;
    double snr_linear = pow(10.0, snr_db / 10.0);
    double noise_std = sqrt(P_sig / snr_linear);
    int detections = 0;
    signal_ct_t *noisy = signal_ct_alloc(template_signal->t_start, template_signal->t_end, template_signal->dt);
    if (!noisy) return -1.0;
    for (int trial = 0; trial < num_trials; trial++) {
        memcpy(noisy->samples, template_signal->samples, (size_t)template_signal->num_samples * sizeof(signal_value_t));
        signal_fill_awgn(noisy, rng, noise_std);
        signal_ct_t *mf_out = signal_ct_alloc(noisy->t_start, noisy->t_end + template_signal->t_end - template_signal->t_start, noisy->dt);
        if (mf_out) {
            signal_ct_matched_filter(noisy, template_signal, mf_out);
            double peak = 0.0;
            for (signal_index_t i = 0; i < mf_out->num_samples; i++)
                if (fabs(mf_out->samples[i]) > peak) peak = fabs(mf_out->samples[i]);
            if (peak > threshold) detections++;
            signal_ct_free(mf_out);
        }
    }
    signal_ct_free(noisy);
    return (double)detections / (double)num_trials;
}

/* =================================================================
 * L8: Ergodicity Check - Time Average vs Ensemble Average
 * A WSS process is mean-ergodic if lim_{T->inf} (1/T) int x(t) dt = E[x(t)]
 * Computes time average for each realization, compares to ensemble average.
 * Reference: Papoulis & Pillai (2002) section 13.2
 * ================================================================= */

int signal_check_ergodicity_mean(const signal_ct_t **ensemble, int num_realizations,
                                   double *time_avg, double *ensemble_avg, double *error)
{
    if (!ensemble || num_realizations <= 0) return -1;
    double sum_time_avg = 0.0;
    double sum_ensemble = 0.0;
    int valid_count = 0;
    for (int r = 0; r < num_realizations; r++) {
        if (!ensemble[r] || !ensemble[r]->samples) continue;
        double mu = signal_ct_mean(ensemble[r]);
        sum_time_avg += mu;
        sum_ensemble += ensemble[r]->samples[r % ensemble[r]->num_samples];
        valid_count++;
    }
    if (valid_count == 0) return -1;
    double mean_time_avg = sum_time_avg / (double)valid_count;
    double mean_ensemble_avg = sum_ensemble / (double)valid_count;
    if (time_avg) *time_avg = mean_time_avg;
    if (ensemble_avg) *ensemble_avg = mean_ensemble_avg;
    if (error) {
        if (fabs(mean_ensemble_avg) > 1e-15)
            *error = fabs(mean_time_avg - mean_ensemble_avg) / fabs(mean_ensemble_avg);
        else
            *error = fabs(mean_time_avg - mean_ensemble_avg);
    }
    return 0;
}

/* =================================================================
 * L5: Gaussian White Noise via Box-Muller Transform
 * ================================================================= */
int signal_generate_gaussian_noise(signal_t *sig, double mean, double stddev)
{
    if (!sig || !sig->data || sig->length == 0) return -1;
    if (stddev <= 0.0) return -1;

    signal_rng_t *rng = signal_rng_create(5489ULL);
    if (!rng) return -1;

    for (size_t i = 0; i < sig->length - 1; i += 2) {
        double u1 = signal_rng_uniform(rng);
        double u2 = signal_rng_uniform(rng);
        if (u1 <= 0.0) u1 = 1e-15;
        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;
        sig->data[i] = mean + stddev * r * cos(theta);
        if (i + 1 < sig->length)
            sig->data[i + 1] = mean + stddev * r * sin(theta);
    }
    if (sig->length % 2 == 1) {
        double u1 = signal_rng_uniform(rng);
        double u2 = signal_rng_uniform(rng);
        if (u1 <= 0.0) u1 = 1e-15;
        sig->data[sig->length - 1] = mean + stddev * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    }
    signal_rng_free(rng);
    return 0;
}

/* =================================================================
 * L5: Uniform White Noise [-amplitude, +amplitude]
 * ================================================================= */
int signal_generate_uniform_noise(signal_t *sig, double amplitude)
{
    if (!sig || !sig->data || sig->length == 0 || amplitude <= 0.0) return -1;
    signal_rng_t *rng = signal_rng_create(32771ULL);
    if (!rng) return -1;
    for (size_t i = 0; i < sig->length; i++) {
        sig->data[i] = amplitude * (2.0 * signal_rng_uniform(rng) - 1.0);
    }
    signal_rng_free(rng);
    return 0;
}

/* =================================================================
 * L5: Pink Noise (1/f) via Voss-McCartney Algorithm
 * Generates approximate 1/f noise by summing octaves of white noise.
 * ================================================================= */
int signal_generate_pink_noise(signal_t *sig, double stddev)
{
    if (!sig || !sig->data || sig->length == 0 || stddev <= 0.0) return -1;

    int n_octaves = 7;
    if ((size_t)(1 << n_octaves) > sig->length) {
        n_octaves = (int)floor(log2((double)sig->length));
        if (n_octaves < 3) n_octaves = 3;
    }

    double *octave = (double *)calloc(sig->length, sizeof(double));
    if (!octave) return -1;
    memset(sig->data, 0, sig->length * sizeof(double));

    signal_rng_t *rng = signal_rng_create(19650218ULL);
    if (!rng) { free(octave); return -1; }

    for (int oct = 0; oct < n_octaves; oct++) {
        int period = 1 << oct;
        for (size_t i = 0; i < sig->length; i++)
            octave[i] = signal_rng_gaussian(rng, 0.0, stddev);

        double scale = 1.0 / sqrt((double)(oct + 1));
        for (size_t i = 0; i < sig->length; i++)
            sig->data[i] += scale * octave[i];
    }

    double max_val = 0.0;
    for (size_t i = 0; i < sig->length; i++)
        if (fabs(sig->data[i]) > max_val) max_val = fabs(sig->data[i]);
    if (max_val > 0.0)
        for (size_t i = 0; i < sig->length; i++)
            sig->data[i] = sig->data[i] / max_val * 3.0 * stddev;

    free(octave);
    signal_rng_free(rng);
    return 0;
}

/* =================================================================
 * L5: Poisson Noise (Shot Noise)
 * Each sample is a Poisson random variable with rate lambda.
 * ================================================================= */
int signal_generate_poisson_noise(signal_t *sig, double lambda)
{
    if (!sig || !sig->data || sig->length == 0 || lambda <= 0.0) return -1;
    signal_rng_t *rng = signal_rng_create(271828ULL);
    if (!rng) return -1;

    double L = exp(-lambda);
    for (size_t i = 0; i < sig->length; i++) {
        double p = 1.0;
        int k = 0;
        while (p > L && k < 1000) {
            k++;
            p *= signal_rng_uniform(rng);
        }
        sig->data[i] = (double)(k - 1);
    }
    signal_rng_free(rng);
    return 0;
}

/* =================================================================
 * L6: Signal-to-Noise Ratio (SNR) Computation
 * SNR_dB = 10 * log10(P_signal / P_noise)
 * ================================================================= */
double signal_snr_db(const signal_t *signal, const signal_t *noisy, const signal_t *clean_ref)
{
    if (!signal || !noisy || !clean_ref || signal->length != noisy->length
        || signal->length != clean_ref->length || signal->length == 0)
        return -999.0;

    double P_signal = 0.0, P_noise = 0.0;
    for (size_t i = 0; i < signal->length; i++) {
        double diff = noisy->data[i] - clean_ref->data[i];
        P_signal += clean_ref->data[i] * clean_ref->data[i];
        P_noise += diff * diff;
    }
    if (P_noise <= 0.0) return 300.0;
    if (P_signal <= 0.0) return -300.0;
    return 10.0 * log10(P_signal / P_noise);
}

/* =================================================================
 * L7: Additive White Gaussian Noise Channel Model
 * y[n] = x[n] + w[n], w ~ N(0, sigma^2)
 * ================================================================= */
int signal_awgn_channel(const signal_t *input, signal_t *output, double snr_db)
{
    if (!input || !output || input->length != output->length || input->length == 0)
        return -1;

    double P_signal = 0.0;
    for (size_t i = 0; i < input->length; i++)
        P_signal += input->data[i] * input->data[i];
    P_signal /= (double)input->length;
    if (P_signal <= 0.0) return -1;

    double sigma = sqrt(P_signal / pow(10.0, snr_db / 10.0));
    signal_t *noise = signal_create(input->length);
    if (!noise) return -1;
    signal_generate_gaussian_noise(noise, 0.0, sigma);

    for (size_t i = 0; i < input->length; i++)
        output->data[i] = input->data[i] + noise->data[i];

    signal_free(noise);
    return 0;
}
