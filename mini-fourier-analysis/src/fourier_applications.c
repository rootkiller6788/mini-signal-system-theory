/**
 * @file    fourier_applications.c
 * @brief   Fourier Analysis Applications — L7 Applications + L8 Advanced
 *
 * @details Real-world applications of Fourier analysis in signal processing
 *          engineering domains:
 *
 *          Audio:         octave/third-octave spectrum analysis, THD
 *          Vibration:     bearing fault detection (BPFO/BPFI/BSF/FTF)
 *          Radar:         LFM chirp generation, pulse compression
 *          GPS:           C/A code acquisition via FFT correlation
 *          OFDM:          Peak-to-Average Power Ratio analysis
 *          Power quality: harmonic decomposition (IEEE 519)
 *          Speech:        formant estimation via cepstral method
 *          Sparse FT:     sub-Nyquist sparse Fourier transform
 *
 * Knowledge covered:
 *   L7: Audio, vibration, radar, GPS, OFDM, power quality, speech applications
 *   L8: Sparse Fourier transform (SFT)
 *
 * Reference: Lyons (2010), Ch.13; Misra & Enge (2011) for GPS;
 *            Haykin & Moher (2009) for OFDM; Hassanieh et al. (2012) for SFT.
 */

#include "fourier_applications.h"
#include "fourier_fft.h"
#include "fourier_spectrum.h"
#include "fourier_window.h"
#include "fourier_convolution.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "fourier_app: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── L7: Audio Octave Spectrum ─────────────────────────────────────────── */

int audio_octave_spectrum(const double *samples, size_t num_samples,
                            double sample_rate, int bands_per_octave,
                            double **band_energy, double **band_centers,
                            int *num_bands) {
    if (!samples || num_samples == 0 || sample_rate <= 0.0
        || !band_energy || !band_centers || !num_bands) return -1;

    if (bands_per_octave < 1) bands_per_octave = 1;

    /* Standard octave band edges per IEC 61260:
     * Center frequencies: 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000 (Hz)
     * For 1/3 octave: 25, 31.5, 40, 50, 63, 80, 100, 125, ... up to 20000
     */

    /* Define base center frequency (1 kHz) and work outward */
    double f_ref = 1000.0;
    double ratio = pow(2.0, 1.0 / (double)bands_per_octave);

    /* Count how many bands fit below Nyquist */
    int count = 0;
    double f_c = f_ref;
    /* Go down from reference */
    while (f_c > 20.0) { f_c /= ratio; }
    f_c *= ratio;  /* first valid band */
    while (f_c < sample_rate / 2.0 && count < 100) {
        count++;
        f_c *= ratio;
    }

    if (count == 0) return -1;

    *num_bands = count;
    *band_energy = (double*)safe_malloc(count * sizeof(double));
    *band_centers = (double*)safe_malloc(count * sizeof(double));

    /* Compute PSD via Welch for good frequency resolution */
    size_t seg_len = (num_samples > 4096) ? 4096 : num_samples;
    psd_result_t psd = welch_psd(samples, num_samples, sample_rate, seg_len, 0.5, WINDOW_HANN);

    /* Band edges: each band from f_c/√2 to f_c·√2 for octave,
     * or f_c·2^{-1/(2B)} to f_c·2^{+1/(2B)} for 1/B octave */
    f_c = f_ref;
    while (f_c / ratio >= 20.0) { f_c /= ratio; }

    for (int b = 0; b < count; b++) {
        double f_low  = f_c / pow(2.0, 1.0 / (2.0 * (double)bands_per_octave));
        double f_high = f_c * pow(2.0, 1.0 / (2.0 * (double)bands_per_octave));

        (*band_centers)[b] = f_c;

        /* Integrate PSD over band */
        double energy = 0.0;
        if (psd.psd) {
            for (size_t k = 0; k < psd.N; k++) {
                double f_k = (double)k * sample_rate / (double)psd.N;
                if (f_k >= f_low && f_k < f_high) {
                    energy += psd.psd[k];
                }
            }
        }
        (*band_energy)[b] = energy * sample_rate / (double)psd.N;  /* scale by Δf */

        f_c *= ratio;
    }

    psd_result_free(&psd);
    return 0;
}

/* ─── L7: Total Harmonic Distortion ─────────────────────────────────────── */

double total_harmonic_distortion(const double *x, size_t N, double fs,
                                   double nominal_freq, int num_harmonics) {
    if (!x || N == 0 || fs <= 0.0 || nominal_freq <= 0.0) return -1.0;

    /* Compute FFT */
    dft_result_t dft = dft_direct(x, N, fs);
    if (!dft.X) return -1.0;

    amplitude_phase_spectrum_t aps = dft_to_amplitude_phase(&dft);

    /* Find fundamental bin */
    size_t fund_bin = (size_t)(nominal_freq * (double)N / fs + 0.5);
    if (fund_bin >= N) fund_bin = N / 2;

    /* Search for peak near nominal frequency (±2 bins) */
    double V1 = 0.0;
    size_t best_bin = fund_bin;
    for (size_t k = (fund_bin > 2 ? fund_bin - 2 : 0);
         k <= fund_bin + 2 && k < N; k++) {
        if (aps.amplitude[k] > V1) {
            V1 = aps.amplitude[k];
            best_bin = k;
        }
    }

    if (V1 < 1e-20) {
        amplitude_phase_spectrum_free(&aps);
        dft_result_free(&dft);
        return -1.0;
    }

    /* Sum harmonic powers */
    double harmonic_power = 0.0;
    for (int h = 2; h <= num_harmonics; h++) {
        size_t h_bin = best_bin * (size_t)h;
        if (h_bin < N) {
            harmonic_power += aps.amplitude[h_bin] * aps.amplitude[h_bin];
        }
    }

    double thd = sqrt(harmonic_power) / V1;

    amplitude_phase_spectrum_free(&aps);
    dft_result_free(&dft);

    return thd;
}

/* ─── L7: Bearing Fault Detection ──────────────────────────────────────── */

void bearing_fault_detect(const double *spectrum, const double *freq_axis,
                            size_t num_bins, double shaft_freq,
                            int n_balls, double ball_diam, double pitch_diam,
                            double contact_angle,
                            double *bpfo_peak, double *bpfi_peak,
                            double *severity_idx) {
    if (bpfo_peak) *bpfo_peak = 0.0;
    if (bpfi_peak) *bpfi_peak = 0.0;
    if (severity_idx) *severity_idx = 0.0;

    if (!spectrum || !freq_axis || num_bins == 0 || shaft_freq <= 0.0) return;

    /* Compute characteristic fault frequencies */
    double ratio = ball_diam / pitch_diam;
    double cos_phi = cos(contact_angle);

    double bpfo_freq = ((double)n_balls / 2.0) * shaft_freq * (1.0 - ratio * cos_phi);
    double bpfi_freq = ((double)n_balls / 2.0) * shaft_freq * (1.0 + ratio * cos_phi);
    double bsf_freq  = (pitch_diam / (2.0 * ball_diam)) * shaft_freq * (1.0 - ratio * ratio * cos_phi * cos_phi);
    double ftf_freq  = (shaft_freq / 2.0) * (1.0 - ratio * cos_phi);

    /* Search for peaks near characteristic frequencies (±2% tolerance) */
    double tol = 0.02;
    double bpfo_val = 0.0, bpfi_val = 0.0, bsf_val = 0.0;
    double noise_floor = 0.0;

    (void)ftf_freq;  /* Cage frequency computed for reference */

    for (size_t i = 0; i < num_bins; i++) {
        double f = freq_axis[i];
        double amp = spectrum[i];

        /* Accumulate noise floor estimate */
        noise_floor += amp;

        /* Check near BPFO */
        if (fabs(f - bpfo_freq) / bpfo_freq < tol && amp > bpfo_val) bpfo_val = amp;
        /* Check near 2× BPFO (often more prominent) */
        if (fabs(f - 2.0 * bpfo_freq) / (2.0 * bpfo_freq) < tol && amp > bpfo_val) bpfo_val = amp;

        /* Check near BPFI */
        if (fabs(f - bpfi_freq) / bpfi_freq < tol && amp > bpfi_val) bpfi_val = amp;
        if (fabs(f - 2.0 * bpfi_freq) / (2.0 * bpfi_freq) < tol && amp > bpfi_val) bpfi_val = amp;

        /* Check near BSF */
        if (fabs(f - bsf_freq) / bsf_freq < tol && amp > bsf_val) bsf_val = amp;
    }

    noise_floor /= (double)num_bins;

    if (bpfo_peak) *bpfo_peak = bpfo_val;
    if (bpfi_peak) *bpfi_peak = bpfi_val;

    /* Severity index: ratio of fault peaks to noise floor */
    if (severity_idx) {
        double max_fault = bpfo_val;
        if (bpfi_val > max_fault) max_fault = bpfi_val;
        if (bsf_val > max_fault) max_fault = bsf_val;
        *severity_idx = (noise_floor > 1e-20) ? max_fault / noise_floor : 0.0;
    }
}

/* ─── L7: LFM Chirp Generation ─────────────────────────────────────────── */

int chirp_generate_lfm(double duration, double bandwidth, double sample_rate,
                         size_t *num_samples,
                         double **chirp_real, double **chirp_imag) {
    if (duration <= 0.0 || bandwidth <= 0.0 || sample_rate <= 0.0
        || !num_samples || !chirp_real || !chirp_imag) return -1;

    *num_samples = (size_t)(duration * sample_rate);
    if (*num_samples < 2) { *num_samples = 2; }

    *chirp_real = (double*)safe_malloc(*num_samples * sizeof(double));
    *chirp_imag = (double*)safe_malloc(*num_samples * sizeof(double));

    double chirp_rate = bandwidth / duration;  /* Hz/s */
    double dt = 1.0 / sample_rate;

    for (size_t n = 0; n < *num_samples; n++) {
        double t = (double)n * dt;
        /* Instantaneous phase: φ(t) = π·(B/T)·t² = π·chirp_rate·t² */
        double phase = M_PI * chirp_rate * t * t;
        (*chirp_real)[n] = cos(phase);
        (*chirp_imag)[n] = sin(phase);
    }

    return 0;
}

/* ─── L7: Pulse Compression (Matched Filter) ────────────────────────────── */

void pulse_compression_matched_filter(const double *received, size_t rx_len,
                                        const fourier_complex_t *tx_template,
                                        size_t tx_len,
                                        double *compressed) {
    if (!received || !tx_template || !compressed || rx_len == 0 || tx_len == 0) return;

    /* Matched filter = time-reversed complex conjugate of transmitted pulse */
    double *mf = (double*)safe_malloc(tx_len * sizeof(double));
    for (size_t i = 0; i < tx_len; i++) {
        mf[i] = creal(conj(tx_template[tx_len - 1 - i]));
    }

    /* Cross-correlation via convolution with time-reversed template */
    convolution_linear(received, rx_len, mf, tx_len, compressed);

    free(mf);
}

/* ─── L7: GPS Acquisition (FFT-based circular correlation) ──────────────── */

void gps_acquisition_fft(const double *if_signal,
                           const int *ca_code,
                           size_t sig_len, size_t code_len,
                           double doppler_hz, double if_freq,
                           double sample_rate,
                           double *correlation,
                           size_t *peak_code_phase,
                           double *peak_ratio) {
    if (!if_signal || !ca_code || !correlation || !peak_code_phase || !peak_ratio
        || sig_len == 0 || code_len == 0) return;

    /* Strip carrier: multiply by exp(-j·2π·(IF+doppler)·t) */
    fourier_complex_t *baseband = (fourier_complex_t*)safe_malloc(sig_len * sizeof(fourier_complex_t));
    double freq = if_freq + doppler_hz;
    double dt = 1.0 / sample_rate;

    for (size_t n = 0; n < sig_len; n++) {
        double phase = -2.0 * M_PI * freq * (double)n * dt;
        fourier_complex_t carrier = cos(phase) + I * sin(phase);
        baseband[n] = if_signal[n] * carrier;
    }

    /* Code replica (complex, padded to FFT size) */
    size_t fft_len = 1;
    while (fft_len < sig_len) fft_len <<= 1;

    fourier_complex_t *code_fft = (fourier_complex_t*)calloc(fft_len, sizeof(fourier_complex_t));
    fourier_complex_t *sig_fft  = (fourier_complex_t*)calloc(fft_len, sizeof(fourier_complex_t));
    if (!code_fft || !sig_fft) { free(baseband); free(code_fft); free(sig_fft); return; }

    for (size_t n = 0; n < code_len && n < sig_len; n++) code_fft[n] = (double)ca_code[n];
    for (size_t n = 0; n < sig_len; n++) sig_fft[n] = baseband[n];

    fft_plan_t *pf = fft_plan_create(fft_len, 0);
    fft_plan_t *pi = fft_plan_create(fft_len, 1);

    fft_execute_dit(pf, code_fft);
    fft_execute_dit(pf, sig_fft);

    /* Circular cross-correlation:  R = IFFT{FFT{sig}·conj(FFT{code})} */
    for (size_t i = 0; i < fft_len; i++) sig_fft[i] *= conj(code_fft[i]);

    fft_execute_dit(pi, sig_fft);

    /* Extract correlation magnitude */
    for (size_t n = 0; n < code_len && n < fft_len; n++) {
        correlation[n] = cabs(sig_fft[n]);
    }

    /* Find peak and second peak */
    double peak = 0.0, second = 0.0;
    size_t peak_idx = 0;
    for (size_t n = 0; n < code_len; n++) {
        if (correlation[n] > peak) {
            second = peak;
            peak = correlation[n];
            peak_idx = n;
        } else if (correlation[n] > second) {
            second = correlation[n];
        }
    }

    *peak_code_phase = peak_idx;
    *peak_ratio = (second > 1e-20) ? peak / second : peak;

    free(baseband);
    free(code_fft);
    free(sig_fft);
    fft_plan_destroy(pf);
    fft_plan_destroy(pi);
}

/* ─── L7: OFDM PAPR ─────────────────────────────────────────────────────── */

void ofdm_papr(const fourier_complex_t *subcarriers, size_t N,
                int oversample, double *papr_db) {
    if (!subcarriers || !papr_db || N == 0) return;

    if (oversample < 1) oversample = 4;

    size_t L = N * (size_t)oversample;
    fourier_complex_t *padded = (fourier_complex_t*)calloc(L, sizeof(fourier_complex_t));
    if (!padded) return;

    /* Place subcarriers at center of spectrum for oversampled IFFT */
    for (size_t k = 0; k < N / 2; k++) padded[k] = subcarriers[k];
    for (size_t k = 0; k < N / 2; k++) padded[L - N/2 + k] = subcarriers[N/2 + k];

    /* Oversampled IFFT: s[n] = (1/√N)·Σ d_k·exp(j·2π·k·n/L) */
    fft_plan_t *pi = fft_plan_create(L, 1);  /* inverse FFT */
    fft_execute_dit(pi, padded);

    /* Find peak and average power */
    double peak_power = 0.0, avg_power = 0.0;
    for (size_t n = 0; n < L; n++) {
        double pwr = creal(padded[n]) * creal(padded[n]) + cimag(padded[n]) * cimag(padded[n]);
        if (pwr > peak_power) peak_power = pwr;
        avg_power += pwr;
    }
    avg_power /= (double)L;

    if (papr_db) {
        *papr_db = (avg_power > 1e-20) ? 10.0 * log10(peak_power / avg_power) : 0.0;
    }

    free(padded);
    fft_plan_destroy(pi);
}

/* ─── L7: OFDM PAPR CCDF ────────────────────────────────────────────────── */

double ofdm_papr_ccdf(size_t N, double threshold_db) {
    if (N == 0) return 1.0;

    double gamma = pow(10.0, threshold_db / 10.0);  /* linear threshold */
    /* Theoretical CCDF for N independent subcarriers:
     * Pr(PAPR > γ) = 1 - (1 - e^{-γ})^N */
    double p = 1.0 - exp(-gamma);
    double ccdf = 1.0 - pow(p, (double)N);
    if (ccdf < 0.0) ccdf = 0.0;
    if (ccdf > 1.0) ccdf = 1.0;
    return ccdf;
}

/* ─── L7: Power Quality Harmonics ───────────────────────────────────────── */

void power_harmonics(const double *v, size_t N, double fs,
                       double fund_freq, size_t max_harm,
                       double *harmonics, double *thd_percent) {
    if (!v || !harmonics || !thd_percent || N == 0 || fs <= 0.0) return;

    memset(harmonics, 0, max_harm * sizeof(double));
    *thd_percent = 0.0;

    /* Compute FFT */
    dft_result_t dft = dft_direct(v, N, fs);
    if (!dft.X) return;

    amplitude_phase_spectrum_t aps = dft_to_amplitude_phase(&dft);

    /* Extract harmonic amplitudes at multiples of fundamental */
    double V1 = 0.0;
    for (size_t h = 1; h <= max_harm; h++) {
        double f_h = (double)h * fund_freq;
        size_t bin = (size_t)(f_h * (double)N / fs + 0.5);
        if (bin < N / 2) {
            harmonics[h - 1] = aps.amplitude[bin];
            if (h == 1) V1 = harmonics[h - 1];
        }
    }

    /* THD = √(Σ_{h=2}^{K} V_h²) / V₁ × 100% */
    if (V1 > 1e-20) {
        double sum_h2 = 0.0;
        for (size_t h = 1; h < max_harm; h++) {
            sum_h2 += harmonics[h] * harmonics[h];
        }
        *thd_percent = sqrt(sum_h2) / V1 * 100.0;
    }

    amplitude_phase_spectrum_free(&aps);
    dft_result_free(&dft);
}

/* ─── L7: Speech Formant Estimation ─────────────────────────────────────── */

void speech_formants(const double *speech, size_t N, double fs,
                       size_t num_formants,
                       double *formants, double *bandwidths) {
    if (!speech || !formants || !bandwidths || N == 0 || num_formants == 0) return;

    /* Compute real cepstrum */
    double *cepstrum = (double*)safe_malloc(N * sizeof(double));
    real_cepstrum(speech, N, cepstrum);

    /* Lifter: keep low quefrencies (spectral envelope).
     * Typical cutoff for vocal tract: 3 ms → quefrency = 0.003 * fs */
    size_t cutoff = (size_t)(0.003 * fs + 0.5);
    if (cutoff >= N / 2) cutoff = N / 2 - 1;
    if (cutoff < 1) cutoff = 1;

    double *cepstrum_low = (double*)safe_malloc(N * sizeof(double));
    /* Rectangular lifter: keep quefrencies 0..cutoff */
    memset(cepstrum_low, 0, N * sizeof(double));
    for (size_t n = 0; n < cutoff; n++) {
        cepstrum_low[n] = cepstrum[n];
        if (n > 0 && n < N - 1) cepstrum_low[N - n] = cepstrum[n];  /* symmetry */
    }

    /* Transform back to frequency domain → spectral envelope */
    dft_result_t envelope_dft = dft_direct(cepstrum_low, N, fs);
    if (!envelope_dft.X) { free(cepstrum); free(cepstrum_low); return; }

    amplitude_phase_spectrum_t envelope = dft_to_amplitude_phase(&envelope_dft);

    /* Find formants as peaks in the spectral envelope */
    /* Search the positive half-spectrum for peaks */
    size_t half = N / 2;
    double *peak_freqs = (double*)safe_malloc(half * sizeof(double));
    double *peak_amps  = (double*)safe_malloc(half * sizeof(double));
    size_t num_peaks = 0;

    for (size_t k = 2; k < half - 1 && num_peaks < half; k++) {
        double a0 = envelope.amplitude[k];
        double aL = envelope.amplitude[k - 1];
        double aR = envelope.amplitude[k + 1];
        /* Peak detection: local maximum */
        if (a0 > aL && a0 > aR && a0 > 1e-6) {
            /* Quadratic interpolation for more accurate peak location */
            double delta = 0.5 * (aL - aR) / (aL - 2.0 * a0 + aR);
            double f_peak = ((double)k + delta) * fs / (double)N;
            peak_freqs[num_peaks] = f_peak;
            peak_amps[num_peaks] = a0 - 0.25 * (aL - aR) * delta;
            num_peaks++;
        }
    }

    /* Select top num_formants by amplitude in valid range (50-5000 Hz) */
    size_t selected = 0;
    for (size_t attempt = 0; attempt < num_peaks && selected < num_formants; attempt++) {
        /* Find max peak not yet selected */
        double max_amp = -1.0;
        size_t max_idx = 0;
        for (size_t i = 0; i < num_peaks; i++) {
            if (peak_freqs[i] >= 50.0 && peak_freqs[i] <= 5000.0
                && peak_amps[i] > max_amp) {
                max_amp = peak_amps[i];
                max_idx = i;
            }
        }
        if (max_amp < 0.0) break;

        formants[selected] = peak_freqs[max_idx];

        /* Estimate bandwidth via 3dB width around peak */
        size_t pk_bin = (size_t)(peak_freqs[max_idx] * (double)N / fs + 0.5);
        double pk_amp = envelope.amplitude[pk_bin];
        double bw_low = peak_freqs[max_idx], bw_high = peak_freqs[max_idx];
        for (size_t k = pk_bin; k > 0; k--) {
            if (envelope.amplitude[k] < pk_amp / sqrt(2.0)) { bw_low = (double)k * fs / (double)N; break; }
        }
        for (size_t k = pk_bin; k < half; k++) {
            if (envelope.amplitude[k] < pk_amp / sqrt(2.0)) { bw_high = (double)k * fs / (double)N; break; }
        }
        bandwidths[selected] = bw_high - bw_low;

        peak_amps[max_idx] = -1.0;  /* mark as used */
        selected++;
    }

    /* For unselected formants, fill with 0 */
    for (size_t i = selected; i < num_formants; i++) {
        formants[i] = 0.0;
        bandwidths[i] = 0.0;
    }

    free(cepstrum);
    free(cepstrum_low);
    free(peak_freqs);
    free(peak_amps);
    amplitude_phase_spectrum_free(&envelope);
    dft_result_free(&envelope_dft);
}

/* ─── L8: Sparse Fourier Transform ──────────────────────────────────────── */

int sparse_fourier_transform(const double *x, size_t N, size_t k,
                               fourier_complex_t *X_sparse, size_t *freq_bins) {
    if (!x || !X_sparse || !freq_bins || N == 0 || k == 0 || k >= N) return 0;

    /* Initialize output to zero */
    memset(X_sparse, 0, N * sizeof(fourier_complex_t));

    /* Simplified SFT:
     * 1. Bin the frequencies via random permutation + flat window
     * 2. Identify "heavy" bins
     * 3. Estimate coefficients within these bins
     *
     * This basic implementation uses a spectral peak detection approach:
     * compute full FFT (or pruned if k ≪ N), identify k largest peaks. */

    /* Since we don't have a full SFT implementation (which requires randomized
     * bin permutation, aliasing filter, and iterative estimation), we provide
     * a peak-detection approach that demonstrates the principle:
     * find k largest DFT coefficients from a pruned computation.
     *
     * For sub-Nyquist SFT, see Hassanieh et al. (2012) for the full algorithm
     * with O(k·log N) complexity.  This implementation uses O(N·log N) FFT
     * with frequency-domain peak picking. */

    size_t N_pow2 = 1;
    while (N_pow2 < N) N_pow2 <<= 1;

    fourier_complex_t *fft_out = (fourier_complex_t*)calloc(N_pow2, sizeof(fourier_complex_t));
    if (!fft_out) return 0;

    for (size_t i = 0; i < N; i++) fft_out[i] = x[i];

    fft_plan_t *plan = fft_plan_create(N_pow2, 0);
    if (!plan) { free(fft_out); return 0; }
    fft_execute_dit(plan, fft_out);
    fft_plan_destroy(plan);

    /* Select k bins with largest magnitude */
    /* Simple selection via sorting indices by magnitude (partial sort) */
    size_t *indices = (size_t*)safe_malloc(N_pow2 * sizeof(size_t));
    for (size_t i = 0; i < N_pow2; i++) indices[i] = i;

    /* Simple selection sort for top k (O(k·N) is fine when k ≪ N) */
    for (size_t i = 0; i < k && i < N_pow2; i++) {
        double max_mag = -1.0;
        size_t max_idx = i;
        for (size_t j = i; j < N_pow2; j++) {
            double mag = cabs(fft_out[indices[j]]);
            if (mag > max_mag) { max_mag = mag; max_idx = j; }
        }
        /* Swap */
        size_t tmp = indices[i];
        indices[i] = indices[max_idx];
        indices[max_idx] = tmp;
    }

    /* Only keep bins within original N */
    size_t found = 0;
    for (size_t i = 0; i < k && i < N_pow2; i++) {
        if (indices[i] < N) {
            freq_bins[found] = indices[i];
            X_sparse[indices[i]] = fft_out[indices[i]];
            found++;
        }
    }

    free(fft_out);
    free(indices);
    return (int)found;
}

/* ─── L8: Orthogonal Matching Pursuit (Compressive Sensing) ───────────────── */

/**
 * OMP for sparse signal recovery from compressed measurements.
 *
 * Solves y = Φ·x where x is k-sparse using greedy atom selection
 * and least-squares refinement at each iteration.
 *
 * The core insight of compressive sensing: a k-sparse signal of dimension N
 * can be recovered from M = O(k·log(N/k)) ≪ N random measurements.
 *
 * Knowledge covered:
 *   L8: Compressive sensing, greedy recovery algorithms, OMP
 *   L5: Least-squares via normal equations, Cholesky-based LS solve
 *   L4: Restricted Isometry Property (RIP) — implicit in random Φ design
 */
int omp_sparse_recovery(const double *Phi, const double *y,
                         size_t M, size_t N, size_t k,
                         double *x_hat, size_t *support,
                         double tolerance) {
    if (!Phi || !y || !x_hat || !support || M == 0 || N == 0 || k == 0 || k > N)
        return -1;

    if (tolerance <= 0.0) tolerance = 1e-6;

    /* Initialize x̂ = 0, residual r = y */
    memset(x_hat, 0, N * sizeof(double));
    double *residual = (double*)safe_malloc(M * sizeof(double));
    memcpy(residual, y, M * sizeof(double));

    /* Track selected atoms */
    int *selected = (int*)calloc(k, sizeof(int));
    if (!selected) { free(residual); return -1; }
    int num_selected = 0;

    /* Working space for LS solve */
    double *Phi_sub = (double*)safe_malloc(M * k * sizeof(double));  /* submatrix columns */
    double *G      = (double*)safe_malloc(k * k * sizeof(double));   /* Gram matrix G = Φ_subᵀ·Φ_sub */
    double *b      = (double*)safe_malloc(k * sizeof(double));       /* rhs = Φ_subᵀ·y */
    double *w      = (double*)safe_malloc(k * sizeof(double));       /* LS solution for selected atoms */
    double *corr   = (double*)safe_malloc(N * sizeof(double));       /* correlation scores */

    for (size_t iter = 0; iter < k; iter++) {

        /* ── Step 1: Compute correlation of each atom with residual ── */
        double max_corr = -1.0;
        size_t best_atom = 0;

        for (size_t j = 0; j < N; j++) {
            /* Check if already selected */
            int already = 0;
            for (int s = 0; s < num_selected; s++) {
                if (selected[s] == (int)j) { already = 1; break; }
            }
            if (already) { corr[j] = -1.0; continue; }

            /* Correlation = |⟨φ_j, r⟩| = Σ_{i=0}^{M-1} Phi[i·N+j] · r[i] */
            double dot = 0.0;
            for (size_t i = 0; i < M; i++) {
                dot += Phi[i * N + j] * residual[i];
            }
            double abs_dot = fabs(dot);
            corr[j] = abs_dot;
            if (abs_dot > max_corr) {
                max_corr = abs_dot;
                best_atom = j;
            }
        }

        if (max_corr < tolerance) break;  /* residual too small */

        /* ── Step 2: Add best atom to support ── */
        selected[num_selected] = (int)best_atom;
        num_selected++;

        /* Copy selected atom columns into Phi_sub */
        for (int s = 0; s < num_selected; s++) {
            for (size_t i = 0; i < M; i++) {
                Phi_sub[i * k + (size_t)s] = Phi[i * N + (size_t)selected[s]];
            }
        }

        /* ── Step 3: Solve least-squares min ||y - Φ_sub·w||² ──
         * Normal equations: (Φ_subᵀ·Φ_sub)·w = Φ_subᵀ·y
         *   G = Φ_subᵀ·Φ_sub  (size m×m where m = num_selected)
         *   b = Φ_subᵀ·y
         * Solve G·w = b via Cholesky (G is symmetric positive-definite
         * for linearly independent columns).
         */

        /* Build G and b */
        size_t m = (size_t)num_selected;
        memset(G, 0, k * k * sizeof(double));
        memset(b, 0, k * sizeof(double));

        for (size_t r = 0; r < m; r++) {
            /* b[r] = Σ_i Phi_sub[i,r] · y[i] */
            double bsum = 0.0;
            for (size_t i = 0; i < M; i++) {
                bsum += Phi_sub[i * k + r] * y[i];
            }
            b[r] = bsum;

            for (size_t c = r; c < m; c++) {
                /* G[r,c] = Σ_i Phi_sub[i,r] · Phi_sub[i,c] */
                double gsum = 0.0;
                for (size_t i = 0; i < M; i++) {
                    gsum += Phi_sub[i * k + r] * Phi_sub[i * k + c];
                }
                G[r * k + c] = gsum;
                G[c * k + r] = gsum;  /* symmetric */
            }
        }

        /* Cholesky decomposition: G = L·Lᵀ (in-place)
         *   L[j,j] = √(G[j,j] - Σ_{p=0}^{j-1} L[j,p]²)
         *   L[i,j] = (G[i,j] - Σ_{p=0}^{j-1} L[i,p]·L[j,p]) / L[j,j]  for i > j
         */
        double *L = (double*)safe_malloc(k * k * sizeof(double));
        memset(L, 0, k * k * sizeof(double));
        int chol_ok = 1;

        for (size_t j = 0; j < m; j++) {
            double sdiag = G[j * k + j];
            for (size_t p = 0; p < j; p++) {
                double ljp = L[j * k + p];
                sdiag -= ljp * ljp;
            }
            if (sdiag <= 1e-20) { chol_ok = 0; break; }
            L[j * k + j] = sqrt(sdiag);

            for (size_t i = j + 1; i < m; i++) {
                double soff = G[i * k + j];
                for (size_t p = 0; p < j; p++) {
                    soff -= L[i * k + p] * L[j * k + p];
                }
                L[i * k + j] = soff / L[j * k + j];
            }
        }

        if (!chol_ok) {
            /* Fallback: use pseudoinverse via normal equations with regularization */
            double lambda = 1e-8;
            for (size_t r = 0; r < m; r++) {
                double bsum = 0.0;
                for (size_t i = 0; i < M; i++) {
                    bsum += Phi_sub[i * k + r] * y[i];
                }
                /* Match residual columns with y */
                double wr = 0.0;
                double denom = 0.0;
                for (size_t i = 0; i < M; i++) {
                    double col = Phi_sub[i * k + r];
                    denom += col * col;
                }
                /* Simple projection: w_r = ⟨φ_r, y⟩ / ⟨φ_r, φ_r⟩ */
                if (denom > 1e-20) {
                    wr = bsum / (denom + lambda);
                }
                w[r] = wr;
            }
        } else {
            /* Forward substitution: L·z = b */
            double *z = (double*)safe_malloc(k * sizeof(double));
            for (size_t i = 0; i < m; i++) {
                double sum_lz = 0.0;
                for (size_t p = 0; p < i; p++) {
                    sum_lz += L[i * k + p] * z[p];
                }
                z[i] = (b[i] - sum_lz) / L[i * k + i];
            }

            /* Back substitution: Lᵀ·w = z */
            for (size_t i = m; i > 0; i--) {
                size_t ii = i - 1;
                double sum_lw = 0.0;
                for (size_t p = ii + 1; p < m; p++) {
                    sum_lw += L[p * k + ii] * w[p];
                }
                w[ii] = (z[ii] - sum_lw) / L[ii * k + ii];
            }
            free(z);
        }

        free(L);

        /* ── Step 4: Update x̂ with current LS solution ── */
        for (int s = 0; s < num_selected; s++) {
            x_hat[selected[s]] = w[s];
        }

        /* ── Step 5: Update residual r = y - Φ·x̂ ── */
        memcpy(residual, y, M * sizeof(double));
        for (int s = 0; s < num_selected; s++) {
            double coeff = w[s];
            for (size_t i = 0; i < M; i++) {
                residual[i] -= coeff * Phi[i * N + selected[s]];
            }
        }

        /* Check convergence by residual norm */
        double r_norm = 0.0;
        for (size_t i = 0; i < M; i++) r_norm += residual[i] * residual[i];
        r_norm = sqrt(r_norm);
        if (r_norm < tolerance * sqrt((double)M)) break;
    }

    /* Copy support indices */
    for (int s = 0; s < num_selected; s++) {
        support[s] = (size_t)selected[s];
    }
    for (int s = num_selected; s < (int)k; s++) {
        support[s] = 0;
    }

    free(selected);
    free(residual);
    free(Phi_sub);
    free(G);
    free(b);
    free(w);
    free(corr);

    return num_selected;
}

/**
 * Build a random-row DFT sensing matrix for compressive spectrum sensing.
 *
 * Φ[i][j] = cos(2π·row[i]·j/N)  for real part (stored at even indices)
 * Φ[i][j] = -sin(2π·row[i]·j/N) for imag part (stored at odd indices)
 *
 * Actually stores as interleaved real/imag in a flat array:
 *   Phi[2*(i*N + j)]     = cos(-2π·row[i]·j/N) =  cos(2π·row[i]·j/N)
 *   Phi[2*(i*N + j) + 1] = sin(-2π·row[i]·j/N) = -sin(2π·row[i]·j/N)
 */
void omp_dft_sensing_matrix(size_t M, size_t N, const size_t *row_indices,
                              double *Phi) {
    if (!Phi || !row_indices || M == 0 || N == 0) return;

    for (size_t i = 0; i < M; i++) {
        size_t row = row_indices[i] % N;
        for (size_t j = 0; j < N; j++) {
            double angle = 2.0 * M_PI * (double)(row * j) / (double)N;
            Phi[2 * (i * N + j)]     = cos(angle);
            Phi[2 * (i * N + j) + 1] = -sin(angle);
        }
    }
}

/* ─── Memory Management ─────────────────────────────────────────────────── */

void app_free_vector(double *v) { free(v); }
void app_free_complex_vector(fourier_complex_t *v) { free(v); }
