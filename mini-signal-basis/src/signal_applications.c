#include "signal_applications.h"
#include "signal_basis.h"
#include "signal_correlation.h"
#include "signal_ops.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =================================================================
 * L7: Goertzel Algorithm for single-frequency detection
 * More efficient than full DFT for detecting few frequencies.
 * Reference: Proakis & Manolakis (2007) section 7.2.2
 * ================================================================= */

double signal_goertzel_detector(const signal_ct_t *x, double target_freq)
{
    if (!x || !x->samples || x->num_samples < 2) return 0.0;
    double fs = 1.0 / x->dt;
    double omega = 2.0 * M_PI * target_freq / fs;
    double coeff = 2.0 * cos(omega);
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    signal_index_t N = x->num_samples;
    for (signal_index_t n = 0; n < N; n++) {
        s0 = x->samples[n] + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    double real = s1 - s2 * cos(omega);
    double imag = s2 * sin(omega);
    return sqrt(real * real + imag * imag) / (double)N;
}

/* =================================================================
 * L7: DTMF Tone Detection
 * Standard DTMF frequencies (ITU-T Q.23):
 * Rows: 697, 770, 852, 941 Hz
 * Cols: 1209, 1336, 1477, 1633 Hz
 * ================================================================= */

int signal_detect_dtmf(const signal_ct_t *audio, double sample_rate, dtmf_tone_t *result)
{
    (void)sample_rate;
    if (!audio || !audio->samples || !result) return -1;
    double row_freqs[4] = {697.0, 770.0, 852.0, 941.0};
    double col_freqs[4] = {1209.0, 1336.0, 1477.0, 1633.0};
    char dtmf_map[4][4] = {
        {'1', '2', '3', 'A'},
        {'4', '5', '6', 'B'},
        {'7', '8', '9', 'C'},
        {'*', '0', '#', 'D'}
    };

    int best_ri = 0, best_ci = 0;
    double max_confidence = 0.0;
    for (int ri = 0; ri < 4; ri++) {
        double rv = signal_goertzel_detector(audio, row_freqs[ri]);
        for (int ci = 0; ci < 4; ci++) {
            double cv = signal_goertzel_detector(audio, col_freqs[ci]);
            double conf = rv * cv;
            if (conf > max_confidence) {
                max_confidence = conf;
                best_ri = ri; best_ci = ci;
            }
        }
    }
    result->row_freq = row_freqs[best_ri];
    result->col_freq = col_freqs[best_ci];
    result->digit = dtmf_map[best_ri][best_ci];
    result->confidence = max_confidence;
    return 0;
}

/* =================================================================
 * L7: ECG QRS Detection (simplified Pan-Tompkins)
 * Reference: Pan & Tompkins, IEEE TBME (1985)
 * ================================================================= */

int signal_detect_ecg_qrs(const signal_ct_t *ecg, double sample_rate, ecg_detection_t *result)
{
    if (!ecg || !ecg->samples || !result) return -1;
    signal_index_t N = ecg->num_samples;
    double *filtered = (double *)malloc((size_t)N * sizeof(double));
    double *deriv = (double *)malloc((size_t)N * sizeof(double));
    double *squared = (double *)malloc((size_t)N * sizeof(double));
    double *integrated = (double *)malloc((size_t)N * sizeof(double));
    if (!filtered || !deriv || !squared || !integrated) { free(filtered); free(deriv); free(squared); free(integrated); return -1; }
    int win_bp = (int)(0.15 * sample_rate);
    if (win_bp < 2) win_bp = 2;
    for (signal_index_t i = win_bp; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < win_bp; j++) sum += ecg->samples[i - j];
        filtered[i] = ecg->samples[i] - sum / (double)win_bp;
    }
    for (signal_index_t i = 1; i < N; i++)
        deriv[i] = filtered[i] - filtered[i - 1];
    for (signal_index_t i = 0; i < N; i++)
        squared[i] = deriv[i] * deriv[i];
    int win_int = (int)(0.15 * sample_rate);
    if (win_int < 1) win_int = 1;
    for (signal_index_t i = win_int; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < win_int; j++) sum += squared[i - j];
        integrated[i] = sum / (double)win_int;
    }
    double max_val = 0.0;
    for (signal_index_t i = 0; i < N; i++)
        if (integrated[i] > max_val) max_val = integrated[i];
    double threshold = max_val * 0.3;
    int max_beats = N / 2;
    double *btimes = (double *)malloc((size_t)max_beats * sizeof(double));
    int beat_count = 0;
    int refractory = (int)(0.2 * sample_rate);
    signal_index_t last_beat = (signal_index_t)(-refractory);
    for (signal_index_t i = 1; i < N - 1; i++) {
        if (integrated[i] > threshold && integrated[i] > integrated[i-1] && integrated[i] >= integrated[i+1]) {
            if ((i - last_beat) > refractory && beat_count < max_beats) {
                btimes[beat_count] = ecg->t_start + (double)i * ecg->dt;
                beat_count++;
                last_beat = i;
            }
        }
    }
    result->num_beats = beat_count;
    result->beat_times = btimes;
    if (beat_count >= 2) {
        double total_rr = 0.0;
        for (int i = 1; i < beat_count; i++) total_rr += btimes[i] - btimes[i-1];
        double mean_rr = total_rr / (double)(beat_count - 1);
        result->mean_heart_rate = 60.0 / mean_rr;
        double var = 0.0;
        for (int i = 1; i < beat_count; i++) {
            double diff = (btimes[i] - btimes[i-1]) - mean_rr;
            var += diff * diff;
        }
        result->heart_rate_std = sqrt(var / (double)(beat_count - 1));
    } else {
        result->mean_heart_rate = 0.0;
        result->heart_rate_std = 0.0;
    }
    free(filtered); free(deriv); free(squared); free(integrated);
    return 0;
}

void ecg_detection_free(ecg_detection_t *result)
{
    if (result) { free(result->beat_times); result->beat_times = NULL; }
}

/* =================================================================
 * L7: Communication Preamble Detection via Matched Filter
 * ================================================================= */

int signal_detect_preamble(const signal_ct_t *received, const signal_ct_t *preamble, double threshold, preamble_detection_t *result)
{
    if (!received || !preamble || !result) return -1;
    result->detected = 0;
    result->detection_time = 0.0;
    result->correlation_peak = 0.0;
    result->snr_estimate = 0.0;
    signal_ct_t *mf_out = signal_ct_alloc(received->t_start, received->t_end + preamble->t_end - preamble->t_start, received->dt);
    if (!mf_out) return -1;
    signal_ct_matched_filter(received, preamble, mf_out);
    signal_index_t N = mf_out->num_samples;
    for (signal_index_t i = 0; i < N; i++) {
        if (fabs(mf_out->samples[i]) > result->correlation_peak) {
            result->correlation_peak = fabs(mf_out->samples[i]);
            result->detection_time = mf_out->t_start + (double)i * mf_out->dt;
        }
    }
    if (result->correlation_peak > threshold) {
        result->detected = 1;
        double noise_power = 0.0; int noise_count = 0;
        for (signal_index_t i = 0; i < N; i++) {
            if (fabs(mf_out->samples[i]) < threshold * 0.3) {
                noise_power += mf_out->samples[i] * mf_out->samples[i];
                noise_count++;
            }
        }
        if (noise_count > 0) {
            double noise_p = noise_power / (double)noise_count;
            if (noise_p > 1e-15)
                result->snr_estimate = 10.0 * log10(result->correlation_peak * result->correlation_peak / noise_p);
        }
    }
    signal_ct_free(mf_out);
    return 0;
}

/* =================================================================
 * L7: Radar Pulse Compression
 * ================================================================= */

int signal_radar_pulse_compression(const signal_ct_t *received, const signal_ct_t *transmitted, double carrier_freq_hz, radar_pulse_compression_t *result)
{
    (void)carrier_freq_hz;
    if (!received || !transmitted || !result) return -1;
    signal_ct_t *mf_out = signal_ct_alloc(received->t_start, received->t_end + transmitted->t_end - transmitted->t_start, received->dt);
    if (!mf_out) return -1;
    signal_ct_matched_filter(received, transmitted, mf_out);
    result->num_range_bins = (int)mf_out->num_samples;
    result->range_profile = (double *)malloc((size_t)result->num_range_bins * sizeof(double));
    if (!result->range_profile) { signal_ct_free(mf_out); return -1; }
    double c = 299792458.0;
    result->range_resolution = c / (2.0 * (1.0 / transmitted->dt));
    result->peak_bin = 0;
    result->peak_power_dB = -999.0;
    for (int i = 0; i < result->num_range_bins; i++) {
        result->range_profile[i] = fabs(mf_out->samples[i]);
        if (result->range_profile[i] > result->peak_power_dB) {
            result->peak_bin = i;
            result->peak_power_dB = result->range_profile[i];
        }
    }
    result->peak_range = (c * 0.5) * (received->t_start + (double)result->peak_bin * received->dt);
    if (result->peak_power_dB > 1e-15)
        result->peak_power_dB = 20.0 * log10(result->peak_power_dB);
    signal_ct_free(mf_out);
    return 0;
}

void radar_pulse_compression_free(radar_pulse_compression_t *result)
{
    if (result) { free(result->range_profile); result->range_profile = NULL; }
}

/* =================================================================
 * L8: Chirp Pulse Generation (Linear FM)
 * s(t) = A * cos(2*pi * (f0*t + (B/(2*T))*t^2)), 0 <= t <= T
 * Instantaneous frequency sweeps linearly from f0 to f0+B.
 * Time-bandwidth product: T*B gives pulse compression ratio.
 * Reference: Richards, Scheer & Holm (2010) section 3.2
 * ================================================================= */

int signal_generate_chirp_pulse(signal_ct_t *pulse, double f0, double B, double T, double A)
{
    if (!pulse || T <= 0.0) return -1;
    double chirp_rate = B / (2.0 * T);
    for (signal_index_t k = 0; k < pulse->num_samples; k++) {
        double t = pulse->t_start + (double)k * pulse->dt;
        if (t >= 0.0 && t <= T) {
            double phase = 2.0 * M_PI * (f0 * t + chirp_rate * t * t);
            pulse->samples[k] = A * cos(phase);
        } else {
            pulse->samples[k] = 0.0;
        }
    }
    return 0;
}

/* =================================================================
 * L8: Barker Code Pulse Generation (BPSK)
 * Barker codes have optimal aperiodic autocorrelation.
 * Known codes: N=2(+ -), 3(+ + -), 4(+ + - +), 5(+ + + - +),
 *             7(+ + + - - + -), 11(+ + + - - - + - - + -),
 *             13(+ + + + + - - + + - + - +)
 * Each chip is a rect pulse of duration chip_duration.
 * Reference: Barker (1953)
 * ================================================================= */

int signal_generate_barker_pulse(signal_ct_t *pulse, int code_length, double chip_duration, double A)
{
    if (!pulse || code_length <= 1 || chip_duration <= 0.0) return -1;
    int barker_2[] = {1, -1};
    int barker_3[] = {1, 1, -1};
    int barker_4[] = {1, 1, -1, 1};
    int barker_5[] = {1, 1, 1, -1, 1};
    int barker_7[] = {1, 1, 1, -1, -1, 1, -1};
    int barker_11[] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1};
    int barker_13[] = {1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1};
    int *code = NULL;
    switch (code_length) {
        case 2:  code = barker_2; break;
        case 3:  code = barker_3; break;
        case 4:  code = barker_4; break;
        case 5:  code = barker_5; break;
        case 7:  code = barker_7; break;
        case 11: code = barker_11; break;
        case 13: code = barker_13; break;
        default: return -1;
    }
    for (signal_index_t k = 0; k < pulse->num_samples; k++) {
        double t = pulse->t_start + (double)k * pulse->dt;
        int chip_idx = (int)(t / chip_duration);
        if (chip_idx >= 0 && chip_idx < code_length)
            pulse->samples[k] = A * (double)code[chip_idx];
        else
            pulse->samples[k] = 0.0;
    }
    return 0;
}
