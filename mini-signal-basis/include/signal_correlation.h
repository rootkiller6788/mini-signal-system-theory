/*
 * signal_correlation.h - Correlation Analysis (L3: Math / L5: Algorithms)
 * Auto/cross-correlation, matched filter, DFT implementation.
 * Course: MIT 6.003 Ch.2 / Stanford EE102A Ch.9 / MIT 6.450 Ch.3
 * Textbook: O&W (1997) section 2.6, Proakis & Salehi (2008) section 4.3
 */
#ifndef SIGNAL_CORRELATION_H
#define SIGNAL_CORRELATION_H
#include "signal_basis.h"

/* L3: Auto-Correlation - R_xx(tau) = integral x(t)*x(t+tau) dt */
int signal_ct_autocorrelation(const signal_ct_t *x, signal_ct_t *rxx);
int signal_ct_autocorrelation_normalized(const signal_ct_t *x, signal_ct_t *rhoxx);
int signal_dt_circular_autocorrelation(const signal_dt_t *x, signal_dt_t *rxx);

/* L3: Cross-Correlation - R_xy(tau) = integral x(t)*y(t+tau) dt */
int signal_ct_crosscorrelation(const signal_ct_t *x, const signal_ct_t *y, signal_ct_t *rxy);
int signal_ct_crosscorrelation_normalized(const signal_ct_t *x, const signal_ct_t *y, signal_ct_t *rhoxy);
int signal_dt_circular_crosscorrelation(const signal_dt_t *x, const signal_dt_t *y, signal_dt_t *rxy);

/* L5: Matched Filter - optimal linear detector for known signal in AWGN */
int signal_ct_matched_filter(const signal_ct_t *received, const signal_ct_t *template_signal, signal_ct_t *output);
double signal_ct_find_best_lag(const signal_ct_t *correlation);

/* L5: Direct DFT (O(N^2) reference implementation) */
int signal_dft_forward(const double *x_real, const double *x_imag, int N, double *X_real, double *X_imag);
int signal_dft_inverse(const double *X_real, const double *X_imag, int N, double *x_real, double *x_imag);
int signal_dft_magnitude(const double *X_real, const double *X_imag, int N, double *magnitude);
int signal_dft_phase(const double *X_real, const double *X_imag, int N, double *phase);
int signal_convolution_via_dft(const double *x, int nx, const double *h, int nh, double *y);


/* L5: Window functions for spectral analysis (pre-DFT preprocessing) */
int signal_apply_hann_window(signal_ct_t *s);
int signal_apply_hamming_window(signal_ct_t *s);
int signal_apply_blackman_window(signal_ct_t *s);

/* L5: Discrete-time cross-correlation via DFT (efficient method) */
int signal_dt_correlation_via_dft(const double *x, const double *y, int N, double *rxy);
int signal_dt_xcorr_via_dft(const signal_dt_t *x, const signal_dt_t *y, signal_dt_t *rxy);

#endif /* SIGNAL_CORRELATION_H */
