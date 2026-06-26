#ifndef STABILITY_ANALYSIS_H
#define STABILITY_ANALYSIS_H
#include "laplace_z_transform_core.h"
#ifdef __cplusplus
extern "C" {
#endif

int stability_routh_array(const double *coeffs, int n, double *routh, int *rows, int *cols);
int stability_routh_rhp_roots(const double *coeffs, int n);
int stability_routh_is_stable(const double *coeffs, int n);
int stability_routh_necessary_condition(const double *coeffs, int n);
int stability_routh_auxiliary(const double *routh, int row, int cols, double *aux_poly, int *aux_deg);
int stability_jury_table(const double *coeffs, int n, double *jury_table, int *rows, int *cols);
int stability_jury_is_stable(const double *coeffs, int n);
double stability_jury_eval_P_at_1(const double *coeffs, int n);
double stability_jury_eval_P_at_minus_1(const double *coeffs, int n);
int stability_from_poles(const complex_t *poles, int n_poles, int is_z_domain);
double stability_margin(const complex_t *poles, int n_poles, int is_z_domain);
int stability_lyapunov_continuous(const double *A, int n);
double stability_gain_margin_db(const rational_func_t *L);
double stability_phase_margin_deg(const rational_func_t *L);

#ifdef __cplusplus
}
#endif
#endif
