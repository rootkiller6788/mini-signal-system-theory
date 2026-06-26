/**
 * @file signal_convolution.c
 * @brief Linear/circular/FFT convolution, deconvolution, step response.
 * L5: O(N^2), overlap-add via FFT.
 * L6: Impulse response estimation via Toeplitz least squares.
 */
#include "signal_basis.h"
#include "signal_ops.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int signal_linear_convolve(const signal_t *x, const signal_t *h, signal_t *y)
{
    if (!x || !h || !y || !x->data || !h->data) return -1;
    size_t N = x->length, M = h->length, L = N + M - 1;
    if (y->length != L) return -1;
    memset(y->data, 0, L * sizeof(double));
    for (size_t n = 0; n < L; n++)
        for (size_t k = 0; k <= n && k < M; k++) {
            size_t j = n - k;
            if (j < N) y->data[n] += x->data[j] * h->data[k];
        }
    return 0;
}

int signal_circular_convolve(const signal_t *x, const signal_t *h, signal_t *y)
{
    if (!x || !h || !y || !x->data || !h->data) return -1;
    if (x->length != h->length || x->length != y->length) return -1;
    size_t N = x->length;
    memset(y->data, 0, N * sizeof(double));
    for (size_t n = 0; n < N; n++)
        for (size_t k = 0; k < N; k++) {
            size_t idx = (n >= k) ? (n - k) : (n + N - k);
            y->data[n] += x->data[idx] * h->data[k];
        }
    return 0;
}

int signal_fft_convolve(const signal_t *x, const signal_t *h, signal_t *y)
{
    if (!x || !h || !y || !x->data || !h->data) return -1;
    size_t N = x->length, M = h->length, L = N + M - 1;
    if (y->length != L) return -1;
    size_t fft_len = 1;
    while (fft_len < L) fft_len <<= 1;
    double *Xr = (double *)calloc(fft_len, sizeof(double));
    double *Xi = (double *)calloc(fft_len, sizeof(double));
    double *Hr = (double *)calloc(fft_len, sizeof(double));
    double *Hi = (double *)calloc(fft_len, sizeof(double));
    if (!Xr || !Xi || !Hr || !Hi) { free(Xr); free(Xi); free(Hr); free(Hi); return -1; }
    for (size_t i = 0; i < N; i++) Xr[i] = x->data[i];
    for (size_t i = 0; i < M; i++) Hr[i] = h->data[i];
    for (size_t s = 1; s < fft_len; s <<= 1)
        for (size_t k = 0; k < fft_len; k += 2 * s)
            for (size_t j = 0; j < s; j++) {
                double a = -M_PI * j / s, wr = cos(a), wi = sin(a);
                double tr = Xr[k+j+s]*wr - Xi[k+j+s]*wi;
                double ti = Xr[k+j+s]*wi + Xi[k+j+s]*wr;
                double hr = Hr[k+j+s]*wr - Hi[k+j+s]*wi;
                double hi = Hr[k+j+s]*wi + Hi[k+j+s]*wr;
                Xr[k+j+s] = Xr[k+j] - tr; Xi[k+j+s] = Xi[k+j] - ti;
                Xr[k+j] += tr; Xi[k+j] += ti;
                Hr[k+j+s] = Hr[k+j] - hr; Hi[k+j+s] = Hi[k+j] - hi;
                Hr[k+j] += hr; Hi[k+j] += hi;
            }
    double *Yr = (double *)calloc(fft_len, sizeof(double));
    double *Yi = (double *)calloc(fft_len, sizeof(double));
    for (size_t i = 0; i < fft_len; i++) {
        Yr[i] = Xr[i]*Hr[i] - Xi[i]*Hi[i];
        Yi[i] = Xr[i]*Hi[i] + Xi[i]*Hr[i];
    }
    for (size_t s = 1; s < fft_len; s <<= 1)
        for (size_t k = 0; k < fft_len; k += 2 * s)
            for (size_t j = 0; j < s; j++) {
                double a = M_PI * j / s, wr = cos(a), wi = sin(a);
                double tr = Yr[k+j+s]*wr - Yi[k+j+s]*wi;
                double ti = Yr[k+j+s]*wi + Yi[k+j+s]*wr;
                Yr[k+j+s] = Yr[k+j] - tr; Yi[k+j+s] = Yi[k+j] - ti;
                Yr[k+j] += tr; Yi[k+j] += ti;
            }
    for (size_t i = 0; i < L; i++) y->data[i] = Yr[i] / (double)fft_len;
    free(Xr); free(Xi); free(Hr); free(Hi); free(Yr); free(Yi);
    return 0;
}

int signal_deconvolve_impulse_response(const signal_t *x, const signal_t *y,
                                        signal_t *h, size_t h_len)
{
    if (!x || !y || !h || !x->data || !y->data || h_len == 0) return -1;
    size_t M = h_len, N = y->length - M + 1;
    if (N > x->length) N = x->length;
    double *Phi = (double *)calloc(M * M, sizeof(double));
    double *rhs = (double *)calloc(M, sizeof(double));
    if (!Phi || !rhs) { free(Phi); free(rhs); return -1; }
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
            double s = 0.0;
            for (size_t n = 0; n < N; n++)
                if (n+j < x->length && n+i < x->length) s += x->data[n+j] * x->data[n+i];
            Phi[i*M+j] = s;
        }
        double s = 0.0;
        for (size_t n = 0; n < N; n++)
            if (n+i < x->length) s += x->data[n+i] * y->data[n+M-1];
        rhs[i] = s;
    }
    for (size_t k = 0; k < M; k++) {
        double mv = fabs(Phi[k*M+k]); size_t mr = k;
        for (size_t i = k+1; i < M; i++)
            if (fabs(Phi[i*M+k]) > mv) { mv = fabs(Phi[i*M+k]); mr = i; }
        if (mv < 1e-15) continue;
        if (mr != k) {
            for (size_t j = k; j < M; j++) {
                double t = Phi[k*M+j]; Phi[k*M+j] = Phi[mr*M+j]; Phi[mr*M+j] = t;
            }
            double t = rhs[k]; rhs[k] = rhs[mr]; rhs[mr] = t;
        }
        for (size_t i = k+1; i < M; i++) {
            double f = Phi[i*M+k] / Phi[k*M+k];
            for (size_t j = k; j < M; j++) Phi[i*M+j] -= f * Phi[k*M+j];
            rhs[i] -= f * rhs[k];
        }
    }
    for (int k = (int)M-1; k >= 0; k--) {
        double s = rhs[k];
        for (size_t j = k+1; j < M; j++) s -= Phi[k*M+j] * h->data[j];
        h->data[k] = (fabs(Phi[k*M+k]) > 1e-15) ? s / Phi[k*M+k] : 0.0;
    }
    free(Phi); free(rhs);
    return 0;
}

int signal_step_response(const signal_t *h, signal_t *step)
{
    if (!h || !step || !h->data || !step->data || h->length == 0) return -1;
    if (step->length != h->length) return -1;
    double a = 0.0;
    for (size_t i = 0; i < h->length; i++) { a += h->data[i]; step->data[i] = a; }
    return 0;
}