#include "signal_correlation.h"
#include "signal_ops.h"
#include "signal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =================================================================
 * L3: Auto-Correlation
 * R_xx(tau) = integral x(t) * x(t+tau) dt
 * Properties: R_xx(0)=E, R_xx(-tau)=R_xx(tau), |R_xx(tau)|<=R_xx(0)
 * ================================================================= */

int signal_ct_autocorrelation(const signal_ct_t *x, signal_ct_t *rxx)
{
    if (!x || !rxx || !x->samples || !rxx->samples) return -1;
    signal_index_t N = x->num_samples;
    signal_index_t Nr = rxx->num_samples;
    for (signal_index_t kl = 0; kl < Nr; kl++) {
        double tau = rxx->t_start + (double)kl * rxx->dt;
        double sum = 0.0;
        int count = 0;
        for (signal_index_t k = 0; k < N; k++) {
            double t = x->t_start + (double)k * x->dt;
            double t2 = t + tau;
            if (t2 >= x->t_start && t2 <= x->t_end) {
                double idx = (t2 - x->t_start) / x->dt;
                signal_index_t i2 = (signal_index_t)idx;
                if (i2 < N - 1) {
                    double frac = idx - (double)i2;
                    double v2 = x->samples[i2] * (1.0 - frac) + x->samples[i2 + 1] * frac;
                    sum += x->samples[k] * v2;
                    count++;
                }
            }
        }
        rxx->samples[kl] = (count > 0) ? sum * x->dt : 0.0;
    }
    return 0;
}

int signal_ct_autocorrelation_normalized(const signal_ct_t *x, signal_ct_t *rhoxx)
{
    if (!x || !rhoxx) return -1;
    double r0 = signal_ct_energy(x);
    if (r0 < 1e-15) {
        if (rhoxx->samples) memset(rhoxx->samples, 0, (size_t)rhoxx->num_samples * sizeof(signal_value_t));
        return 0;
    }
    if (signal_ct_autocorrelation(x, rhoxx) != 0) return -1;
    for (signal_index_t k = 0; k < rhoxx->num_samples; k++)
        rhoxx->samples[k] /= r0;
    return 0;
}

int signal_dt_circular_autocorrelation(const signal_dt_t *x, signal_dt_t *rxx)
{
    if (!x || !rxx || !x->samples || !rxx->samples) return -1;
    signal_index_t N = x->length;
    if (rxx->length != N) return -1;
    for (signal_index_t m = 0; m < N; m++) {
        double sum = 0.0;
        for (signal_index_t n = 0; n < N; n++)
            sum += x->samples[n] * x->samples[(n + m) % N];
        rxx->samples[m] = sum;
    }
    return 0;
}

/* =================================================================
 * L3: Cross-Correlation
 * R_xy(tau) = integral x(t) * y(t+tau) dt
 * ================================================================= */

int signal_ct_crosscorrelation(const signal_ct_t *x, const signal_ct_t *y, signal_ct_t *rxy)
{
    if (!x || !y || !rxy || !x->samples || !y->samples || !rxy->samples) return -1;
    signal_index_t Nx = x->num_samples, Ny = y->num_samples;
    signal_index_t Nr = rxy->num_samples;
    for (signal_index_t kl = 0; kl < Nr; kl++) {
        double tau = rxy->t_start + (double)kl * rxy->dt;
        double sum = 0.0;
        int count = 0;
        for (signal_index_t k = 0; k < Nx; k++) {
            double t = x->t_start + (double)k * x->dt;
            double t2 = t + tau;
            if (t2 >= y->t_start && t2 <= y->t_end) {
                double idx = (t2 - y->t_start) / y->dt;
                signal_index_t i2 = (signal_index_t)idx;
                if (i2 < Ny - 1) {
                    double frac = idx - (double)i2;
                    double v2 = y->samples[i2] * (1.0 - frac) + y->samples[i2 + 1] * frac;
                    sum += x->samples[k] * v2;
                    count++;
                }
            }
        }
        rxy->samples[kl] = (count > 0) ? sum * x->dt : 0.0;
    }
    return 0;
}

int signal_ct_crosscorrelation_normalized(const signal_ct_t *x, const signal_ct_t *y, signal_ct_t *rhoxy)
{
    if (!x || !y || !rhoxy) return -1;
    double Ex = signal_ct_energy(x);
    double Ey = signal_ct_energy(y);
    if (Ex < 1e-15 || Ey < 1e-15) return -1;
    if (signal_ct_crosscorrelation(x, y, rhoxy) != 0) return -1;
    double scale = 1.0 / sqrt(Ex * Ey);
    for (signal_index_t k = 0; k < rhoxy->num_samples; k++)
        rhoxy->samples[k] *= scale;
    return 0;
}

int signal_dt_circular_crosscorrelation(const signal_dt_t *x, const signal_dt_t *y, signal_dt_t *rxy)
{
    if (!x || !y || !rxy || !x->samples || !y->samples || !rxy->samples) return -1;
    signal_index_t N = x->length;
    if (y->length != N || rxy->length != N) return -1;
    for (signal_index_t m = 0; m < N; m++) {
        double sum = 0.0;
        for (signal_index_t n = 0; n < N; n++)
            sum += x->samples[n] * y->samples[(n + m) % N];
        rxy->samples[m] = sum;
    }
    return 0;
}

/* =================================================================
 * L5: Matched Filter
 * h(t) = s(T-t) (time-reversed template), y(t) = (x * h)(t)
 * Maximizes SNR for known signal in AWGN.
 * Reference: Proakis & Salehi (2008) section 4.3
 * ================================================================= */

int signal_ct_matched_filter(const signal_ct_t *received, const signal_ct_t *template_signal, signal_ct_t *output)
{
    if (!received || !template_signal || !output) return -1;
    signal_ct_t *h = signal_ct_alloc(template_signal->t_start, template_signal->t_end, template_signal->dt);
    if (!h) return -1;
    signal_index_t N = template_signal->num_samples;
    for (signal_index_t k = 0; k < N; k++)
        h->samples[k] = template_signal->samples[N - 1 - k];
    int ret = signal_ct_convolve(received, h, output);
    signal_ct_free(h);
    return ret;
}

double signal_ct_find_best_lag(const signal_ct_t *correlation)
{
    if (!correlation || !correlation->samples || correlation->num_samples == 0) return 0.0;
    signal_index_t best_i = 0;
    double best_val = correlation->samples[0];
    for (signal_index_t k = 1; k < correlation->num_samples; k++) {
        if (correlation->samples[k] > best_val) {
            best_val = correlation->samples[k];
            best_i = k;
        }
    }
    return correlation->t_start + (double)best_i * correlation->dt;
}

/* =================================================================
 * L5: Direct DFT - O(N^2) reference implementation
 * X[k] = sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N)
 * Reference: Proakis & Manolakis (2007) section 5.1.1
 * ================================================================= */

int signal_dft_forward(const double *x_real, const double *x_imag, int N, double *X_real, double *X_imag)
{
    if (!x_real || !X_real || !X_imag || N <= 0) return -1;
    double scale = -2.0 * M_PI / (double)N;
    for (int k = 0; k < N; k++) {
        double sum_real = 0.0, sum_imag = 0.0;
        for (int n = 0; n < N; n++) {
            double angle = scale * (double)(k * n);
            double xr = x_real[n];
            double xi = (x_imag) ? x_imag[n] : 0.0;
            double c = cos(angle), s = sin(angle);
            sum_real += xr * c - xi * s;
            sum_imag += xr * s + xi * c;
        }
        X_real[k] = sum_real;
        X_imag[k] = sum_imag;
    }
    return 0;
}

int signal_dft_inverse(const double *X_real, const double *X_imag, int N, double *x_real, double *x_imag)
{
    if (!X_real || !X_imag || !x_real || !x_imag || N <= 0) return -1;
    double scale = 2.0 * M_PI / (double)N;
    for (int n = 0; n < N; n++) {
        double sum_real = 0.0, sum_imag = 0.0;
        for (int k = 0; k < N; k++) {
            double angle = scale * (double)(k * n);
            double c = cos(angle), s = sin(angle);
            sum_real += X_real[k] * c - X_imag[k] * s;
            sum_imag += X_real[k] * s + X_imag[k] * c;
        }
        x_real[n] = sum_real / (double)N;
        x_imag[n] = sum_imag / (double)N;
    }
    return 0;
}

int signal_dft_magnitude(const double *X_real, const double *X_imag, int N, double *magnitude)
{
    if (!X_real || !X_imag || !magnitude || N <= 0) return -1;
    for (int k = 0; k < N; k++)
        magnitude[k] = sqrt(X_real[k] * X_real[k] + X_imag[k] * X_imag[k]);
    return 0;
}

int signal_dft_phase(const double *X_real, const double *X_imag, int N, double *phase)
{
    if (!X_real || !X_imag || !phase || N <= 0) return -1;
    for (int k = 0; k < N; k++)
        phase[k] = atan2(X_imag[k], X_real[k]);
    return 0;
}

/* =================================================================
 * L5: Convolution via DFT
 * y = IFFT(FFT(x) .* FFT(h)), zero-padded to N = Nx + Nh - 1
 * ================================================================= */

int signal_convolution_via_dft(const double *x, int nx, const double *h, int nh, double *y)
{
    if (!x || !h || !y || nx <= 0 || nh <= 0) return -1;
    int N = nx + nh - 1;
    double *xr = (double *)calloc((size_t)N, sizeof(double));
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *hr = (double *)calloc((size_t)N, sizeof(double));
    double *hi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    double *Hr = (double *)malloc((size_t)N * sizeof(double));
    double *Hi = (double *)malloc((size_t)N * sizeof(double));
    double *Yr = (double *)malloc((size_t)N * sizeof(double));
    double *Yi = (double *)malloc((size_t)N * sizeof(double));
    if (!xr || !hr || !Xr || !Xi || !Hr || !Hi || !Yr) {
        free(xr); free(xi); free(hr); free(hi);
        free(Xr); free(Xi); free(Hr); free(Hi); free(Yr); free(Yi);
        return -1;
    }
    memcpy(xr, x, (size_t)nx * sizeof(double));
    memcpy(hr, h, (size_t)nh * sizeof(double));
    signal_dft_forward(xr, xi, N, Xr, Xi);
    signal_dft_forward(hr, hi, N, Hr, Hi);
    for (int k = 0; k < N; k++) {
        double yr = Xr[k] * Hr[k] - Xi[k] * Hi[k];
        double yi = Xr[k] * Hi[k] + Xi[k] * Hr[k];
        Yr[k] = yr; Yi[k] = yi;
    }
    signal_dft_inverse(Yr, Yi, N, xr, xi);
    memcpy(y, xr, (size_t)N * sizeof(double));
    free(xr); free(xi); free(hr); free(hi);
    free(Xr); free(Xi); free(Hr); free(Hi); free(Yr); free(Yi);
    return 0;
}

/* =================================================================
 * L5: Correlation via DFT (circular, efficient for large N)
 * R_xy = IFFT(FFT(x) .* conj(FFT(y)))
 * This is O(N^2) due to direct DFT, but demonstrates the principle.
 * For production use, an FFT library would give O(N log N).
 * ================================================================= */

int signal_dt_correlation_via_dft(const double *x, const double *y, int N, double *rxy)
{
    if (!x || !y || !rxy || N <= 0) return -1;
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *yi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    double *Yr = (double *)malloc((size_t)N * sizeof(double));
    double *Yi = (double *)malloc((size_t)N * sizeof(double));
    if (!xi || !Xr || !Xi || !Yr || !Yi) { free(xi); free(Xr); free(Xi); free(Yr); free(Yi); return -1; }
    signal_dft_forward(x, xi, N, Xr, Xi);
    signal_dft_forward(y, yi, N, Yr, Yi);
    for (int k = 0; k < N; k++) {
        double rr = Xr[k] * Yr[k] + Xi[k] * Yi[k];
        double ri = -Xr[k] * Yi[k] + Xi[k] * Yr[k];
        Xr[k] = rr; Xi[k] = ri;
    }
    signal_dft_inverse(Xr, Xi, N, rxy, xi);
    free(xi); free(Xr); free(Xi); free(Yr); free(Yi);
    return 0;
}

/* =================================================================
 * L5: Window Functions (for spectral analysis)
 * Hann, Hamming, Blackman windows for DFT preprocessing.
 * Reference: Proakis & Manolakis (2007) section 10.2
 * ================================================================= */

int signal_apply_hann_window(signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples < 2) return -1;
    signal_index_t N = s->num_samples;
    for (signal_index_t k = 0; k < N; k++)
        s->samples[k] *= 0.5 * (1.0 - cos(2.0 * M_PI * (double)k / (double)(N - 1)));
    return 0;
}

int signal_apply_hamming_window(signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples < 2) return -1;
    signal_index_t N = s->num_samples;
    for (signal_index_t k = 0; k < N; k++)
        s->samples[k] *= 0.54 - 0.46 * cos(2.0 * M_PI * (double)k / (double)(N - 1));
    return 0;
}

int signal_apply_blackman_window(signal_ct_t *s)
{
    if (!s || !s->samples || s->num_samples < 2) return -1;
    signal_index_t N = s->num_samples;
    for (signal_index_t k = 0; k < N; k++) {
        double a = 2.0 * M_PI * (double)k / (double)(N - 1);
        s->samples[k] *= 0.42 - 0.5 * cos(a) + 0.08 * cos(2.0 * a);
    }
    return 0;
}

/* =================================================================
 * L5: Circular Cross-Correlation via DFT
 * Efficient implementation using convolution theorem:
 * R_xy = IFFT(FFT(x) .* conj(FFT(y)))
 * ================================================================= */

int signal_dt_xcorr_via_dft(const signal_dt_t *x, const signal_dt_t *y, signal_dt_t *rxy)
{
    if (!x || !y || !rxy || !x->samples || !y->samples || !rxy->samples) return -1;
    signal_index_t N = x->length;
    if (y->length != N || rxy->length != N) return -1;
    double *xi = (double *)calloc((size_t)N, sizeof(double));
    double *yi = (double *)calloc((size_t)N, sizeof(double));
    double *Xr = (double *)malloc((size_t)N * sizeof(double));
    double *Xi = (double *)malloc((size_t)N * sizeof(double));
    double *Yr = (double *)malloc((size_t)N * sizeof(double));
    double *Yi = (double *)malloc((size_t)N * sizeof(double));
    if (!xi || !Xr || !Xi || !Yr || !Yi) { free(xi); free(yi); free(Xr); free(Xi); free(Yr); free(Yi); return -1; }
    signal_dft_forward(x->samples, xi, (int)N, Xr, Xi);
    signal_dft_forward(y->samples, yi, (int)N, Yr, Yi);
    for (int k = 0; k < (int)N; k++) {
        double rr = Xr[k] * Yr[k] + Xi[k] * Yi[k];
        double ri = -Xr[k] * Yi[k] + Xi[k] * Yr[k];
        Xr[k] = rr; Xi[k] = ri;
    }
    signal_dft_inverse(Xr, Xi, (int)N, rxy->samples, xi);
    free(xi); free(yi); free(Xr); free(Xi); free(Yr); free(Yi);
    return 0;
}
