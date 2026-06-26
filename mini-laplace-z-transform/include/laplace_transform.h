#ifndef LAPLACE_TRANSFORM_H
#define LAPLACE_TRANSFORM_H
#include "laplace_z_transform_core.h"
#ifdef __cplusplus
extern "C" {
#endif

rational_func_t laplace_impulse(void);
rational_func_t laplace_step(void);
rational_func_t laplace_ramp(void);
rational_func_t laplace_t_power_n(int n);
rational_func_t laplace_exp_decay(double a);
rational_func_t laplace_cos(double w0);
rational_func_t laplace_sin(double w0);
rational_func_t laplace_damped_cos(double a, double w0);
rational_func_t laplace_damped_sin(double a, double w0);
rational_func_t laplace_first_order(double K, double tau);
rational_func_t laplace_second_order(double K, double wn, double zeta);
rational_func_t laplace_linearity(double a, const rational_func_t *F1, double b, const rational_func_t *F2);
rational_func_t laplace_time_shift(const rational_func_t *F, double t0, int pade_order);
rational_func_t laplace_freq_shift(const rational_func_t *F, double a);
rational_func_t laplace_time_scale(const rational_func_t *F, double a);
rational_func_t laplace_convolution(const rational_func_t *F, const rational_func_t *G);
rational_func_t laplace_derivative(const rational_func_t *F);
rational_func_t laplace_integral(const rational_func_t *F);
double laplace_initial_value(const rational_func_t *F);
double laplace_final_value(const rational_func_t *F);
double laplace_parseval_energy(const rational_func_t *F, double w_max, int n_pts);
complex_t laplace_to_fourier(const rational_func_t *F, double omega);
int laplace_is_stable(const rational_func_t *F);
int laplace_roc_includes_jw_axis(const rational_func_t *F);

#ifdef __cplusplus
}
#endif
#endif
