/**
 * @file    adaptive_applications.c
 * @brief   Canonical applications of adaptive filtering
 *
 * Knowledge Coverage:
 *   L6 (Canonical Problems): System identification, noise cancellation,
 *                             echo cancellation, channel equalization,
 *                             adaptive line enhancement
 *   L7 (Applications): ECG noise cancellation, acoustic echo cancellation,
 *                      digital communications channel equalization
 *
 * Reference:
 *   Widrow et al. (1975), "Adaptive Noise Cancelling: Principles and
 *                          Applications", Proc. IEEE
 *   Haykin (2014), "Adaptive Filter Theory", 5th ed., Pearson: Chapter 16
 *   Proakis & Salehi (2008), "Digital Communications", 5th ed., McGraw-Hill
 *
 * Course Mapping:
 *   MIT 6.450 - Digital Communications (adaptive equalization)
 *   Stanford EE264 - Adaptive Filtering (noise/echo cancellation)
 *   Berkeley EE225D - Adaptive Signal Processing
 *   Georgia Tech ECE 6601 - Communications (channel estimation)
 */

#include "adaptive_filter.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =========================================================================
 * L6: System Identification - LMS-based
 * ========================================================================= */

/**
 * @brief af_app_system_identify_lms: System identification using LMS
 *
 * Problem: Given an unknown linear system with impulse response h[n],
 * observe its input x[n] and noisy output d[n] = (h * x)[n] + v[n],
 * and estimate h from the observations.
 *
 * Setup:
 *   Unknown plant:  h = [h0, h1, ..., h_{M-1}]
 *   Input:          x[n] (wideband excitation, e.g., white noise)
 *   Desired output: d[n] = sum_{k=0}^{M-1} h_k * x[n-k] + v[n]
 *   Adaptive filter: w (same order M as unknown plant)
 *
 * The LMS algorithm converges to the Wiener solution:
 *   w_opt = R^{-1} * p = h (since x is uncorrelated with v)
 *
 * Applications:
 *   - Acoustic room impulse response measurement
 *   - Control system plant identification
 *   - Geological/geophysical exploration (seismic deconvolution)
 *   - NMR spectroscopy system identification
 *
 * @param x           Excitation input signal [n_samples]
 * @param d           Observed system output [n_samples]
 * @param n_samples   Number of training samples
 * @param order       Unknown system order
 * @param mu          LMS step-size
 * @param w_out       Output: estimated impulse response [order]
 * @param mse_out     Output: final MSE (may be NULL)
 * @return            0 on success
 */
int af_app_system_identify_lms(const double *x, const double *d,
                                size_t n_samples, size_t order,
                                double mu, double *w_out, double *mse_out) {
    af_config_t cfg;
    af_filter_t filter;
    size_t n;

    if (!x || !d || !w_out || n_samples == 0 || order == 0) return -1;

    /* Configure LMS filter */
    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_LMS;
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return -1;

    /* Adaptive identification */
    for (n = 0; n < n_samples; n++) {
        af_adapt(&filter, x[n], d[n]);
    }

    /* Extract results */
    af_get_weights(&filter, w_out);
    if (mse_out) *mse_out = filter.mse;

    af_free(&filter);
    return 0;
}

/* =========================================================================
 * L6: Adaptive Noise Cancellation
 * ========================================================================= */

/**
 * @brief af_app_noise_cancel_lms: Adaptive noise cancellation
 *
 * Problem: A primary signal s[n] is corrupted by additive noise n0[n].
 * A reference sensor picks up a noise signal n1[n] that is correlated
 * with n0[n] but uncorrelated with s[n]. The adaptive filter uses n1[n]
 * to estimate n0[n], which is then subtracted from the primary input.
 *
 * Configuration:
 *   Primary input:   d[n] = s[n] + n0[n]
 *   Reference input: x[n] = n1[n]  (correlated with n0)
 *   Filter output:   y[n] = estimated n0[n]
 *   System output:   e[n] = d[n] - y[n] 「 s[n]
 *
 * The filter adapts to minimize E[|e|^2]:
 *   E[|e|^2] = E[|s + n0 - y|^2]
 *            = E[|s|^2] + E[|n0 - y|^2]  (since s uncorrelated with n0, n1)
 *
 * Minimizing E[|e|^2] is equivalent to minimizing E[|n0 - y|^2],
 * making y 「 n0 and e 「 s. This is the fundamental principle of
 * adaptive noise cancellation (Widrow et al., 1975).
 *
 * Applications:
 *   - Fetal ECG extraction from maternal ECG (Widrow's classic 1975 paper)
 *   - 60 Hz power line interference removal from ECG/EEG
 *   - Cockpit noise reduction for pilot communications (Boeing)
 *   - Engine noise cancellation in automobiles (Toyota)
 *
 * @param primary     Primary input: s[n] + noise [n_samples]
 * @param reference   Reference noise input [n_samples]
 * @param n_samples   Number of samples
 * @param order       Adaptive filter order
 * @param mu          LMS step-size
 * @param cleaned_out Output: cleaned signal [n_samples]
 * @return            0 on success
 */
int af_app_noise_cancel_lms(const double *primary, const double *reference,
                             size_t n_samples, size_t order, double mu,
                             double *cleaned_out) {
    af_config_t cfg;
    af_filter_t filter;
    size_t n;

    if (!primary || !reference || !cleaned_out || n_samples == 0 || order == 0)
        return -1;

    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_LMS;
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return -1;

    /* Adaptive noise cancellation */
    for (n = 0; n < n_samples; n++) {
        /* Reference = noise, Primary = signal+noise */
        double e = af_adapt(&filter, reference[n], primary[n]);
        cleaned_out[n] = e; /* Error = cleaned signal estimate */
    }

    af_free(&filter);
    return 0;
}

/* =========================================================================
 * L6: Adaptive Echo Cancellation
 * ========================================================================= */

/**
 * @brief af_app_echo_cancel_nlms: Acoustic echo cancellation using NLMS
 *
 * Problem: In hands-free telephony, the far-end speaker's voice is
 * played through the local loudspeaker, and the local microphone
 * picks up both the near-end speaker and the echo of the far-end
 * speaker (after room reverberation). The echo must be cancelled
 * to provide a comfortable conversation.
 *
 * Configuration:
 *   Far-end signal (reference):   x[n] (played through loudspeaker)
 *   Echo path:                    h[n] (room impulse response, unknown)
 *   Echo signal:                  y[n] = (h * x)[n]
 *   Near-end signal + echo:       d[n] = s[n] + y[n]
 *   Adaptive filter output:       y_hat[n] = (w * x)[n] 「 y[n]
 *   Residual (sent to far-end):   e[n] = d[n] - y_hat[n] 「 s[n]
 *
 * NLMS is preferred over LMS for echo cancellation because:
 *   1. Input signal power varies significantly (speech is non-stationary)
 *   2. NLMS adapts step-size to signal power automatically
 *   3. Better convergence for colored speech signals
 *
 * Double-talk detection (not implemented here) is critical to
 * prevent filter divergence when both near-end and far-end speakers
 * are active simultaneously.
 *
 * Applications:
 *   - Hands-free car kits (iPhone, Android automotive mode)
 *   - Video conferencing systems
 *   - VoIP phones
 *   - Smart speakers (Amazon Echo, Google Home)
 *
 * @param far_end     Far-end speech (reference) [n_samples]
 * @param mic_signal  Microphone signal (echo + near-end) [n_samples]
 * @param n_samples   Number of samples
 * @param order       Echo path filter length (typical: 128-1024 taps)
 * @param mu          NLMS step-size (0 < mu < 2)
 * @param residual_out Output: echo-cancelled signal [n_samples]
 * @return            0 on success
 */
int af_app_echo_cancel_nlms(const double *far_end, const double *mic_signal,
                             size_t n_samples, size_t order, double mu,
                             double *residual_out) {
    af_config_t cfg;
    af_filter_t filter;
    size_t n;

    if (!far_end || !mic_signal || !residual_out || n_samples == 0 || order == 0)
        return -1;

    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_NLMS; /* NLMS for robust echo cancellation */
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return -1;

    for (n = 0; n < n_samples; n++) {
        /* far_end = reference (played through speaker)
         * mic = desired (contains echo to cancel) */
        double e = af_adapt(&filter, far_end[n], mic_signal[n]);
        residual_out[n] = e;
    }

    af_free(&filter);
    return 0;
}

/* =========================================================================
 * L6: Adaptive Channel Equalization
 * ========================================================================= */

/**
 * @brief af_app_channel_equalize_lms: Adaptive decision-feedback equalization
 *
 * Problem: In digital communications, the transmitted symbols are
 * distorted by the channel (multipath, dispersion). An adaptive
 * equalizer at the receiver compensates for channel distortion to
 * recover the transmitted symbols.
 *
 * Channel model:
 *   r[n] = sum_{k=0}^{L-1} c_k * a[n-k] + v[n]
 * where a[n] are transmitted symbols, c_k is the channel impulse
 * response of length L, and v[n] is additive noise.
 *
 * Adaptive equalizer:
 *   y[n] = sum_{k=0}^{M-1} w_k * r[n-k]   (linear transversal equalizer)
 *   a_hat[n] = decision(y[n])             (nearest constellation point)
 *   e[n] = a_hat[n] - y[n]                (decision-directed error)
 *
 * In training mode, the true symbol a[n] is used instead of a_hat[n].
 * After convergence, the equalizer switches to decision-directed mode.
 *
 * The LMS algorithm adapts the equalizer taps to minimize MSE:
 *   w(n+1) = w(n) + mu * e(n) * r(n)
 *
 * Applications:
 *   - Digital TV (ATSC 8-VSB equalizer)
 *   - Cable modems (DOCSIS)
 *   - WiFi 802.11 baseband processors
 *   - 5G NR physical layer equalization
 *   - GPS receiver multipath mitigation
 *
 * @param received    Received signal samples [n_samples]
 * @param training    Training sequence [train_len] (known symbols)
 * @param n_samples   Total samples
 * @param train_len   Length of training sequence
 * @param order       Equalizer order (number of taps)
 * @param mu          LMS step-size
 * @param equalized_out Output: equalized soft decisions [n_samples]
 * @return            0 on success
 */
int af_app_channel_equalize_lms(const double *received,
                                 const double *training,
                                 size_t n_samples, size_t train_len,
                                 size_t order, double mu,
                                 double *equalized_out) {
    af_config_t cfg;
    af_filter_t filter;
    size_t n;

    if (!received || !equalized_out || n_samples == 0 || order == 0)
        return -1;

    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_LMS;
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_CENTER_TAP; /* Center tap init for equalizer */

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return -1;

    for (n = 0; n < n_samples; n++) {
        double desired;

        if (n < train_len && training) {
            /* Training mode: use known symbols */
            desired = training[n];
        } else {
            /* Decision-directed mode: use hard decisions */
            double output = af_filter_output(&filter, received[n]);
            /* Simple BPSK decision: sign of output */
            desired = (output >= 0.0) ? 1.0 : -1.0;
        }

        double e = af_adapt(&filter, received[n], desired);
        equalized_out[n] = e + desired; /* Reconstruct soft output */
    }

    af_free(&filter);
    return 0;
}

/* =========================================================================
 * L6: Adaptive Line Enhancement (ALE)
 * ========================================================================= */

/**
 * @brief af_app_line_enhance_lms: Adaptive line enhancement
 *
 * Problem: Detect and enhance a narrowband signal (sinusoid) buried
 * in broadband noise. This is a self-adaptive filter that uses a
 * delayed version of the input as the reference.
 *
 * Configuration:
 *   Input:     x[n] = A * cos(2*pi*f0*n + phi) + v[n]
 *   Reference: x[n - Delta]  (delayed version)
 *   Output:    y[n] 「 A * cos(2*pi*f0*n + phi)  (enhanced sinusoid)
 *
 * The decorrelation delay Delta must be chosen such that:
 *   - The noise samples become uncorrelated (v[n] vs v[n-Delta])
 *   - The sinusoidal component remains correlated
 *
 * Typically Delta = 1 is sufficient, or Delta 「 fs/(4*f0) for
 * better decorrelation of the noise.
 *
 * Applications:
 *   - Detection of weak sine waves in noise (spectrum sensing)
 *   - Carrier recovery in communication receivers
 *   - Power system harmonic detection
 *   - Rotating machinery fault diagnosis (bearing/turbine vibration)
 *
 * @param input       Noisy input signal [n_samples]
 * @param n_samples   Number of samples
 * @param order       ALE filter order
 * @param mu          LMS step-size
 * @param delay       Decorrelation delay (samples)
 * @param enhanced_out Output: enhanced narrowband signal [n_samples]
 * @return            0 on success
 */
int af_app_line_enhance_lms(const double *input, size_t n_samples,
                             size_t order, double mu, size_t delay,
                             double *enhanced_out) {
    af_config_t cfg;
    af_filter_t filter;
    size_t n;

    if (!input || !enhanced_out || n_samples == 0 || order == 0) return -1;

    memset(&cfg, 0, sizeof(cfg));
    cfg.structure = AF_STRUCTURE_FIR;
    cfg.algorithm = AF_ALGO_LMS;
    cfg.order = order;
    cfg.mu = mu;
    cfg.epsilon = 1e-6;
    cfg.lambda = 1.0;
    cfg.delta = 1.0;
    cfg.gamma_leakage = 1.0;
    cfg.init_strategy = AF_INIT_ZEROS;

    memset(&filter, 0, sizeof(filter));
    if (af_init(&cfg, &filter) != 0) return -1;

    for (n = 0; n < n_samples; n++) {
        /* Reference = delayed input (decorrelated noise, correlated signal) */
        size_t ref_idx = (n >= delay) ? (n - delay) : 0;
        af_adapt(&filter, input[ref_idx], input[n]);

        /* Enhanced output = filter output (narrowband estimate) */
        enhanced_out[n] = filter.output;
    }

    af_free(&filter);
    return 0;
}

/* =========================================================================
 * L7: ECG 60 Hz Noise Cancellation Demo
 * ========================================================================= */

/**
 * @brief af_app_ecg_noise_cancel: Simulated ECG power-line noise cancellation
 *
 * Demonstrates the noise cancellation principle on a synthetic ECG
 * signal corrupted by 60 Hz (50 Hz in Europe) power-line interference.
 * The reference is a 60 Hz tone phase-shifted through an unknown
 * transfer function (simulating wiring and grounding effects).
 *
 * This is a canonical biomedical signal processing application of
 * adaptive filtering, directly implementing Widrow's noise cancellation
 * architecture on ECG data.
 *
 * The ECG signal is modeled as a sum of Gaussian pulses (QRS complex,
 * P wave, T wave) at a heart rate of ~72 BPM.
 *
 * @param ecg_clean   Clean ECG signal [n_samples]
 * @param noise_ref   60 Hz reference tone [n_samples]
 * @param n_samples   Number of samples (fs = 360 Hz, ~10 seconds)
 * @param cleaned_out Output: denoised ECG [n_samples]
 * @return            0 on success
 */
int af_app_ecg_noise_cancel(const double *ecg_clean, const double *noise_ref,
                             size_t n_samples, double *cleaned_out) {
    /* Use NLMS with order sufficient to model the interference path.
     * The interference path between the reference sensor and the
     * ECG electrodes is typically short (order ~10-20 taps at 360 Hz). */
    return af_app_noise_cancel_lms(ecg_clean, noise_ref, n_samples,
                                    16, 0.01, cleaned_out);
}

