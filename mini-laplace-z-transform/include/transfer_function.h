#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H
#include "laplace_z_transform_core.h"
#ifdef __cplusplus
extern "C" {
#endif

rational_func_t tf_create(const double *num_coeffs, int num_len, const double *den_coeffs, int den_len);
double tf_impulse_response(const rational_func_t *H, double t);
double tf_step_response(const rational_func_t *H, double t);
double tf_discrete_impulse_response(const rational_func_t *H, int n);
double tf_discrete_step_response(const rational_func_t *H, int n);
void tf_freq_response_mag(const rational_func_t *H, const double *freqs, double *mags, int n_freqs);
void tf_freq_response_phase(const rational_func_t *H, const double *freqs, double *phases, int n_freqs);
double tf_bode_mag_db(const rational_func_t *H, double omega);
double tf_bode_phase_deg(const rational_func_t *H, double omega);
double tf_cutoff_frequency(const rational_func_t *H, double f_low, double f_high);
double tf_dc_gain(const rational_func_t *H, int is_z_domain);
rational_func_t tf_series(const rational_func_t *H1, const rational_func_t *H2);
rational_func_t tf_parallel(const rational_func_t *H1, const rational_func_t *H2);
rational_func_t tf_feedback(const rational_func_t *G, const rational_func_t *H);
rational_func_t tf_feedback_positive(const rational_func_t *G, const rational_func_t *H);
polynomial_t tf_char_eq_feedback(const rational_func_t *G, const rational_func_t *H);
int tf_dominant_pole(const rational_func_t *H, complex_t *dominant_pole);
double tf_settling_time(const rational_func_t *H);
double tf_overshoot_percent(const rational_func_t *H);
int tf_damping_params(const rational_func_t *H, double *zeta, double *wn);

#ifdef __cplusplus
}
#endif
#endif
