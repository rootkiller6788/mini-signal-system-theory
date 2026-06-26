#ifndef PARTIAL_FRACTION_H
#define PARTIAL_FRACTION_H
#include "laplace_z_transform_core.h"
#ifdef __cplusplus
extern "C" {
#endif

int pf_decompose_laplace(const rational_func_t *F, partial_fraction_t *pf);
int pf_decompose_z(const rational_func_t *X, partial_fraction_t *pf);
int pf_residues_simple(const polynomial_t *N, const polynomial_t *D, const complex_t *poles, int n, complex_t *residues);
int pf_residues_repeated(const rational_func_t *F, complex_t pole, int m, complex_t *residues);
rational_func_t pf_to_rational_func(const partial_fraction_t *pf);
int pf_find_poles(const polynomial_t *den, complex_t *poles, int *mult, int *n_distinct);

#ifdef __cplusplus
}
#endif
#endif
