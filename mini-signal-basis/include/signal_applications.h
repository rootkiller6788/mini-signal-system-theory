/*
 * signal_applications.h - Signal Applications (L7: Applications / L8: Advanced)
 *
 * Real-world applications:
 *   - DTMF audio tone detection (Goertzel algorithm)
 *   - ECG QRS complex detection (Pan-Tompkins algorithm)
 *   - Communication preamble detection (matched filter)
 *   - Radar pulse compression (chirp and Barker codes)
 *
 * Course: Stanford EE359 Ch.4 / Michigan EECS 455 Ch.4 / Georgia Tech ECE 6601 Ch.5
 * Reference: Proakis & Salehi (2008) section 4.4-4.5
 *            Richards, Scheer & Holm (2010) section 3.2
 *            Pan & Tompkins, IEEE TBME (1985)
 */
#ifndef SIGNAL_APPLICATIONS_H
#define SIGNAL_APPLICATIONS_H
#include "signal_basis.h"
#include "signal_correlation.h"
#include "signal_decomposition.h"
#include "signal_energy.h"

/* L7: DTMF tone detection */
typedef struct {
    double row_freq; double col_freq;
    char digit; double confidence;
} dtmf_tone_t;

int signal_detect_dtmf(const signal_ct_t *audio, double sample_rate, dtmf_tone_t *result);
double signal_goertzel_detector(const signal_ct_t *x, double target_freq);

/* L7: ECG QRS detection */
typedef struct {
    int num_beats; double *beat_times;
    double mean_heart_rate; double heart_rate_std;
} ecg_detection_t;

int signal_detect_ecg_qrs(const signal_ct_t *ecg, double sample_rate, ecg_detection_t *result);
void ecg_detection_free(ecg_detection_t *result);

/* L7: Communication preamble detection */
typedef struct {
    int detected; double detection_time;
    double correlation_peak; double snr_estimate;
} preamble_detection_t;

int signal_detect_preamble(const signal_ct_t *received, const signal_ct_t *preamble, double threshold, preamble_detection_t *result);

/* L7: Radar pulse compression */
typedef struct {
    double *range_profile; int num_range_bins;
    double range_resolution; int peak_bin;
    double peak_range; double peak_power_dB;
} radar_pulse_compression_t;

int signal_radar_pulse_compression(const signal_ct_t *received, const signal_ct_t *transmitted, double carrier_freq_hz, radar_pulse_compression_t *result);
void radar_pulse_compression_free(radar_pulse_compression_t *result);

/* L8: Pulse generation (advanced waveform design) */
int signal_generate_chirp_pulse(signal_ct_t *pulse, double f0, double B, double T, double A);
int signal_generate_barker_pulse(signal_ct_t *pulse, int code_length, double chip_duration, double A);

#endif /* SIGNAL_APPLICATIONS_H */
