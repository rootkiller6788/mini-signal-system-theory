#include "signal_energy.h"
#include "signal_basis.h"
#include "signal_correlation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* L3: Energy Spectral Density - Psi(f) = |X(f)|^2
 * Computed via DFT magnitude squared. */
int signal_energy_spectral_density(const signal_ct_t *x, double *freqs, double *esd, int N_fft)
{
    if (!x || !x->samples || !freqs || !esd || N_fft <= 0) return -1;
    double *xr = (double *)malloc((size_t)N_fft * sizeof(double));
    double *xi = (double *)calloc((size_t)N_fft, sizeof(double));
    double *Xr = (double *)malloc((size_t)N_fft * sizeof(double));
    double *Xi = (double *)malloc((size_t)N_fft * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return -1; }
    signal_index_t N = x->num_samples;
    signal_index_t n_copy = (N < N_fft) ? N : N_fft;
    for (signal_index_t i = 0; i < n_copy; i++) xr[i] = x->samples[i];
    for (int i = (int)n_copy; i < N_fft; i++) xr[i] = 0.0;
    signal_dft_forward(xr, xi, N_fft, Xr, Xi);
    double fs = 1.0 / x->dt;
    for (int k = 0; k < N_fft; k++) {
        freqs[k] = (double)k * fs / (double)N_fft;
        if (freqs[k] > fs * 0.5) freqs[k] -= fs;
        esd[k] = (Xr[k] * Xr[k] + Xi[k] * Xi[k]) * (x->dt * x->dt);
    }
    free(xr); free(xi); free(Xr); free(Xi);
    return 0;
}

/* L3: Periodogram PSD - S(f) = |X(f)|^2 / N
 * Reference: Proakis & Manolakis (2007) section 14.2.1 */
int signal_psd_periodogram(const signal_ct_t *x, double *freqs, double *psd, int N_fft)
{
    if (!x || !x->samples || !freqs || !psd || N_fft <= 0) return -1;
    double *xr = (double *)malloc((size_t)N_fft * sizeof(double));
    double *xi = (double *)calloc((size_t)N_fft, sizeof(double));
    double *Xr = (double *)malloc((size_t)N_fft * sizeof(double));
    double *Xi = (double *)malloc((size_t)N_fft * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return -1; }
    signal_index_t N = x->num_samples;
    signal_index_t n_copy = (N < N_fft) ? N : N_fft;
    for (signal_index_t i = 0; i < n_copy; i++) xr[i] = x->samples[i];
    for (int i = (int)n_copy; i < N_fft; i++) xr[i] = 0.0;
    signal_dft_forward(xr, xi, N_fft, Xr, Xi);
    double fs = 1.0 / x->dt;
    double scale = 1.0 / (double)N_fft;
    for (int k = 0; k < N_fft; k++) {
        freqs[k] = (double)k * fs / (double)N_fft;
        if (freqs[k] > fs * 0.5) freqs[k] -= fs;
        psd[k] = (Xr[k] * Xr[k] + Xi[k] * Xi[k]) * scale * x->dt;
    }
    free(xr); free(xi); free(Xr); free(Xi);
    return 0;
}

/* L3: Welch PSD Estimation
 * Reference: Welch, IEEE Trans. Audio Electroacoust. (1967) */
static void hann_window(double *w, int N)
{
    for (int i = 0; i < N; i++)
        w[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)(N - 1)));
}

int signal_psd_welch(const signal_ct_t *x, int segment_length, int overlap, double *freqs, double *psd, int N_fft)
{
    if (!x || !x->samples || segment_length <= 0 || overlap < 0 || !freqs || !psd || N_fft <= 0) return -1;
    int step = segment_length - overlap;
    if (step <= 0) return -1;
    double *window = (double *)malloc((size_t)segment_length * sizeof(double));
    double *segment = (double *)malloc((size_t)segment_length * sizeof(double));
    double *win_seg = (double *)calloc((size_t)N_fft, sizeof(double));
    double *zi = (double *)calloc((size_t)N_fft, sizeof(double));
    double *Xr = (double *)malloc((size_t)N_fft * sizeof(double));
    double *Xi = (double *)malloc((size_t)N_fft * sizeof(double));
    if (!window || !segment || !win_seg || !zi || !Xr || !Xi) {
        free(window); free(segment); free(win_seg); free(zi); free(Xr); free(Xi); return -1; }
    hann_window(window, segment_length);
    double win_power = 0.0;
    for (int i = 0; i < segment_length; i++) win_power += window[i] * window[i];
    for (int k = 0; k < N_fft; k++) psd[k] = 0.0;
    int num_segments = 0;
    signal_index_t N = x->num_samples;
    for (int start = 0; start + segment_length <= (int)N; start += step) {
        for (int i = 0; i < segment_length; i++) segment[i] = x->samples[start + i];
        for (int i = 0; i < segment_length; i++) win_seg[i] = segment[i] * window[i];
        for (int i = segment_length; i < N_fft; i++) win_seg[i] = 0.0;
        signal_dft_forward(win_seg, zi, N_fft, Xr, Xi);
        double scale = 1.0 / (win_power * x->dt);
        for (int k = 0; k < N_fft; k++)
            psd[k] += (Xr[k] * Xr[k] + Xi[k] * Xi[k]) * scale;
        num_segments++;
    }
    if (num_segments > 0) {
        double fs = 1.0 / x->dt;
        for (int k = 0; k < N_fft; k++) {
            freqs[k] = (double)k * fs / (double)N_fft;
            if (freqs[k] > fs * 0.5) freqs[k] -= fs;
            psd[k] /= (double)num_segments;
        }
    }
    free(window); free(segment); free(win_seg); free(zi); free(Xr); free(Xi);
    return 0;
}

/* =================================================================
 * L4: Parseval Theorem - integral|x(t)|^2 dt = integral|X(f)|^2 df
 * Reference: O&W (1997) section 4.3, Theorem 4.1
 * ================================================================= */

parseval_result_t signal_verify_parseval(const signal_ct_t *x, double tolerance)
{
    parseval_result_t r = {0.0, 0.0, 0.0, 0};
    if (!x || !x->samples || x->num_samples == 0) return r;
    int N = (int)x->num_samples;
    r.time_energy = signal_ct_energy(x);
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return r; }
    for (int i = 0; i < N; i++) xr[i] = x->samples[i];
    signal_dft_forward(xr, xi, N, Xr, Xi);
    r.freq_energy = 0.0;
    for (int k = 0; k < N; k++)
        r.freq_energy += (Xr[k] * Xr[k] + Xi[k] * Xi[k]) / (double)N;
    r.freq_energy *= x->dt;
    if (r.time_energy > 1e-15)
        r.relative_error = fabs(r.time_energy - r.freq_energy) / r.time_energy;
    r.verified = (r.relative_error < tolerance) ? 1 : 0;
    free(xr); free(xi); free(Xr); free(Xi);
    return r;
}

parseval_result_t signal_verify_parseval_discrete(const double *x, int N, double tolerance)
{
    parseval_result_t r = {0.0, 0.0, 0.0, 0};
    if (!x || N <= 0) return r;
    for (int i = 0; i < N; i++) r.time_energy += x[i] * x[i];
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xi || !Xr || !Xi) { free(xi); free(Xr); free(Xi); return r; }
    signal_dft_forward(x, xi, N, Xr, Xi);
    for (int k = 0; k < N; k++)
        r.freq_energy += (Xr[k] * Xr[k] + Xi[k] * Xi[k]) / (double)N;
    if (r.time_energy > 1e-15)
        r.relative_error = fabs(r.time_energy - r.freq_energy) / r.time_energy;
    r.verified = (r.relative_error < tolerance) ? 1 : 0;
    free(xi); free(Xr); free(Xi);
    return r;
}

/* =================================================================
 * L4: Wiener-Khinchin Theorem - S_xx(f) = F{R_xx(tau)}
 * ================================================================= */

wiener_khinchin_result_t signal_verify_wiener_khinchin(const signal_ct_t *x, double freq, double tolerance)
{
    wiener_khinchin_result_t r = {0.0, 0.0, 0.0, 0};
    if (!x || !x->samples || x->num_samples < 4) return r;
    int N = (int)x->num_samples;
    double fs = 1.0 / x->dt;
    int k_idx = (int)(freq * (double)N / fs + 0.5);
    if (k_idx < 0 || k_idx >= N) return r;
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return r; }
    for (int i = 0; i < N; i++) xr[i] = x->samples[i];
    signal_dft_forward(xr, xi, N, Xr, Xi);
    r.psd_from_periodogram = (Xr[k_idx] * Xr[k_idx] + Xi[k_idx] * Xi[k_idx]) * x->dt / (double)N;
    signal_ct_t *rxx = signal_ct_alloc(-x->t_end, x->t_end, x->dt);
    if (rxx) {
        signal_ct_autocorrelation(x, rxx);
        double *rr = (double *)malloc((size_t)N * sizeof(double));
        if (rr) {
            for (int i = 0; i < N; i++) rr[i] = rxx->samples[rxx->num_samples/2 - N/2 + i];
            signal_dft_forward(rr, xi, N, Xr, Xi);
            r.psd_from_autocorrelation = (Xr[k_idx] * Xr[k_idx] + Xi[k_idx] * Xi[k_idx]) * rxx->dt / (double)N;
            free(rr);
        }
        signal_ct_free(rxx);
    }
    if (r.psd_from_periodogram > 1e-15)
        r.relative_error = fabs(r.psd_from_periodogram - r.psd_from_autocorrelation) / r.psd_from_periodogram;
    r.verified = (r.relative_error < tolerance) ? 1 : 0;
    free(xr); free(xi); free(Xr); free(Xi);
    return r;
}

/* =================================================================
 * L4: Rayleigh Energy Theorem
 * ||x||^2 = sum |<x, phi_k>|^2 for any orthonormal basis
 * ================================================================= */

int signal_verify_rayleigh(const signal_ct_t *x, const signal_ct_t **basis, int K,
                           double *time_energy, double *coeff_sum, double *relative_error)
{
    if (!x || !basis || K <= 0) return -1;
    double E_time = signal_ct_energy(x);
    if (time_energy) *time_energy = E_time;
    double c_sum = 0.0;
    for (int k = 0; k < K; k++) {
        double ck = signal_ct_inner_product(x, basis[k]);
        c_sum += ck * ck;
    }
    if (coeff_sum) *coeff_sum = c_sum;
    if (relative_error) {
        if (E_time > 1e-15)
            *relative_error = fabs(E_time - c_sum) / E_time;
        else
            *relative_error = (c_sum < 1e-15) ? 0.0 : 1.0;
    }
    return 0;
}

/* Energy in a frequency band (fraction of total) */
double signal_energy_in_band(const signal_ct_t *x, double f_low, double f_high)
{
    if (!x || !x->samples || f_low >= f_high) return -1.0;
    int N = (int)x->num_samples;
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return -1.0; }
    for (int i = 0; i < N; i++) xr[i] = x->samples[i];
    (void)xi;
    signal_dft_forward(xr, xi, N, Xr, Xi);
    double fs = 1.0 / x->dt;
    double energy = 0.0, total_energy = 0.0;
    for (int k = 0; k < N; k++) {
        double freq = (double)k * fs / (double)N;
        if (freq > fs * 0.5) freq -= fs;
        double mag2 = Xr[k] * Xr[k] + Xi[k] * Xi[k];
        total_energy += mag2;
        if (fabs(freq) >= f_low && fabs(freq) <= f_high) energy += mag2;
    }
    free(xr); free(xi); free(Xr); free(Xi);
    if (total_energy < 1e-15) return 0.0;
    return energy / total_energy;
}

/* Crest Factor: x_peak / x_rms. Pure sinusoid: sqrt(2) ~= 1.414 */
double signal_ct_crest_factor(const signal_ct_t *s)
{
    if (!s) return 0.0;
    double rms = signal_ct_rms(s);
    if (rms < 1e-15) return 0.0;
    return signal_ct_peak(s) / rms;
}

/* THD: sqrt(sum of harmonic powers) / fundamental magnitude */
double signal_ct_thd(const signal_ct_t *s, double fundamental_freq, int num_harmonics)
{
    if (!s || !s->samples || fundamental_freq <= 0.0 || num_harmonics < 2) return -1.0;
    int N = (int)s->num_samples;
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return -1.0; }
    for (int i = 0; i < N; i++) xr[i] = s->samples[i];
    signal_dft_forward(xr, xi, N, Xr, Xi);
    double fs = 1.0 / s->dt;
    double fundamental_mag = 0.0, harmonic_power = 0.0;
    for (int k = 1; k <= num_harmonics; k++) {
        double freq = (double)k * fundamental_freq;
        int bin = (int)(freq * (double)N / fs + 0.5);
        if (bin >= 0 && bin < N) {
            double mag2 = Xr[bin] * Xr[bin] + Xi[bin] * Xi[bin];
            if (k == 1) fundamental_mag = sqrt(mag2);
            else harmonic_power += mag2;
        }
    }
    free(xr); free(xi); free(Xr); free(Xi);
    if (fundamental_mag < 1e-15) return -1.0;
    return sqrt(harmonic_power) / fundamental_mag;
}

/* Fractional power containment bandwidth */
double signal_ct_bandwidth(const signal_ct_t *s, double power_fraction)
{
    if (!s || !s->samples || power_fraction <= 0.0 || power_fraction > 1.0) return -1.0;
    int N = (int)s->num_samples;
    double *xr = (double *)malloc((size_t)N * sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !xi || !Xr || !Xi) { free(xr); free(xi); free(Xr); free(Xi); return -1.0; }
    for (int i = 0; i < N; i++) xr[i] = s->samples[i];
    signal_dft_forward(xr, xi, N, Xr, Xi);
    double total_power = 0.0;
    for (int k = 0; k < N; k++) total_power += Xr[k] * Xr[k] + Xi[k] * Xi[k];
    if (total_power < 1e-15) { free(xr); free(xi); free(Xr); free(Xi); return 0.0; }
    double target = total_power * power_fraction;
    int center_bin = 0; double max_mag = 0.0;
    for (int k = 0; k < N; k++) {
        double mag = Xr[k] * Xr[k] + Xi[k] * Xi[k];
        if (mag > max_mag) { max_mag = mag; center_bin = k; }
    }
    double accum = Xr[center_bin] * Xr[center_bin] + Xi[center_bin] * Xi[center_bin];
    int half_width = 0;
    while (accum < target && half_width < N / 2) {
        half_width++;
        int idx_r = (center_bin + half_width) % N;
        int idx_l = (center_bin - half_width + N) % N;
        accum += Xr[idx_r] * Xr[idx_r] + Xi[idx_r] * Xi[idx_r];
        if (idx_r != idx_l) accum += Xr[idx_l] * Xr[idx_l] + Xi[idx_l] * Xi[idx_l];
    }
    free(xr); free(xi); free(Xr); free(Xi);
    double fs = 1.0 / s->dt;
    return (double)half_width * fs / (double)N;
}
