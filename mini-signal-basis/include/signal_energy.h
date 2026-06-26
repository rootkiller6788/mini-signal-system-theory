/*
 * signal_energy.h - Signal Energy and Parseval's Theorem (L3-L4)
 * ESD/PSD, Parseval, Wiener-Khinchin, Rayleigh energy theorems.
 * Course: MIT 6.003 Ch.4 / Stanford EE102A Ch.9 / Berkeley EE123 Ch.4
 * Textbook: O&W (1997) section 4.3; Proakis & Manolakis (2007) section 14.1-2
 */
#ifndef SIGNAL_ENERGY_H
#define SIGNAL_ENERGY_H
#include "signal_basis.h"

/* L3: Energy Spectral Density - Psi(f) = |X(f)|^2 */
int signal_energy_spectral_density(const signal_ct_t *x, double *freqs, double *esd, int N_fft);

/* L3: Power Spectral Density estimation */
int signal_psd_periodogram(const signal_ct_t *x, double *freqs, double *psd, int N_fft);
int signal_psd_welch(const signal_ct_t *x, int segment_length, int overlap, double *freqs, double *psd, int N_fft);

/* L4: Parseval's Theorem - integral|x(t)|^2 dt = integral|X(f)|^2 df */
typedef struct {
    double time_energy; double freq_energy;
    double relative_error; int verified;
} parseval_result_t;

parseval_result_t signal_verify_parseval(const signal_ct_t *x, double tolerance);
parseval_result_t signal_verify_parseval_discrete(const double *x, int N, double tolerance);

/* L4: Wiener-Khinchin Theorem - S_xx(f) = F{R_xx(tau)} */
typedef struct {
    double psd_from_periodogram; double psd_from_autocorrelation;
    double relative_error; int verified;
} wiener_khinchin_result_t;

wiener_khinchin_result_t signal_verify_wiener_khinchin(const signal_ct_t *x, double freq, double tolerance);

/* L4: Rayleigh's Energy Theorem - ||x||^2 = sum |<x,phi_k>|^2 */
int signal_verify_rayleigh(const signal_ct_t *x, const signal_ct_t **basis, int K,
                           double *time_energy, double *coeff_sum, double *relative_error);


/* L4: Band-limited energy ratio, crest factor, THD, fractional bandwidth */
double signal_energy_in_band(const signal_ct_t *x, double f_low, double f_high);
double signal_ct_crest_factor(const signal_ct_t *s);
double signal_ct_thd(const signal_ct_t *s, double fundamental_freq, int num_harmonics);
double signal_ct_bandwidth(const signal_ct_t *s, double power_fraction);

#endif /* SIGNAL_ENERGY_H */
