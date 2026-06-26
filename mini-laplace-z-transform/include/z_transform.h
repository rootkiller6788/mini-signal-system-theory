#ifndef Z_TRANSFORM_H
#define Z_TRANSFORM_H
#include "laplace_z_transform_core.h"
#ifdef __cplusplus
extern "C" {
#endif

rational_func_t zt_impulse(void);
rational_func_t zt_step(void);
rational_func_t zt_geometric(double a);
rational_func_t zt_ramp_geometric(double a);
rational_func_t zt_cos(double w0);
rational_func_t zt_sin(double w0);
rational_func_t zt_damped_cos(double r, double w0);
rational_func_t zt_damped_sin(double r, double w0);
rational_func_t zt_first_order(double b, double a);
rational_func_t zt_second_order(double b0, double b1, double b2, double a1, double a2);
rational_func_t zt_linearity(double a, const rational_func_t *X1, double b, const rational_func_t *X2);
rational_func_t zt_time_shift_right(const rational_func_t *X, int n0);
rational_func_t zt_time_shift_left(const rational_func_t *X, int n0);
rational_func_t zt_exp_multiply(const rational_func_t *X, double a);
rational_func_t zt_n_multiply(const rational_func_t *X);
rational_func_t zt_convolution(const rational_func_t *X, const rational_func_t *H);
double zt_initial_value(const rational_func_t *X);
double zt_final_value(const rational_func_t *X);
double zt_parseval_energy(const rational_func_t *X, int n_pts);
double zt_roc_radius(const rational_func_t *X);
int zt_is_stable(const rational_func_t *X);
int zt_roc_includes_unit_circle(const rational_func_t *X);

#ifdef __cplusplus
}
#endif
#endif
