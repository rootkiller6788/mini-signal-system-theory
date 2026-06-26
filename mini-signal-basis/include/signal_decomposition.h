/*
 * signal_decomposition.h - Signal Decomposition (L3/L5/L6)
 * Orthonormal basis representation, Gram-Schmidt, least-squares approximation.
 * Course: MIT 6.003 Ch.3 / Stanford EE102A Ch.4 / ETH 227-0427 Ch.4
 * Textbook: O&W (1997) section 3.2-3.3; Strang & Nguyen (1996) Ch.1
 */
#ifndef SIGNAL_DECOMPOSITION_H
#define SIGNAL_DECOMPOSITION_H
#include "signal_basis.h"

typedef signal_ct_t basis_function_t;

typedef struct {
    basis_function_t **functions;
    int count;
    int capacity;
} basis_set_t;

typedef struct {
    double *coefficients;
    signal_ct_t *approximation;
    int num_terms;
    double residual_energy;
    double energy_ratio;
} signal_approximation_t;

/* L5: Gram-Schmidt orthogonalization (Modified GS for stability) */
int signal_gram_schmidt(const signal_ct_t **vectors, int num_vectors, basis_set_t *orthonormal_basis);

/* L5: Basis projection and reconstruction */
int signal_project_onto_basis(const signal_ct_t *x, const basis_set_t *basis, double *coefficients);
int signal_reconstruct_from_basis(const double *coefficients, int K, const basis_set_t *basis, signal_ct_t *reconstructed);
int signal_truncated_approximation(const signal_ct_t *x, const basis_set_t *basis, int K, signal_approximation_t *result);
void signal_approximation_free(signal_approximation_t *approx);

/* L5: Basis set management */
basis_set_t *basis_set_alloc(int capacity);
int basis_set_add(basis_set_t *bs, const signal_ct_t *phi);
void basis_set_free(basis_set_t *bs);
int basis_set_create_fourier(basis_set_t *bs, double T, int num_harmonics, signal_time_t dt);
int basis_set_create_haar(basis_set_t *bs, signal_time_t duration, signal_time_t dt, int max_level);

/* L6: Least-squares approximation via normal equations G*c=b, G_ij=<v_i,v_j> */
int signal_least_squares_approx(const signal_ct_t *x, const signal_ct_t **templates, int K, signal_approximation_t *result);

#endif /* SIGNAL_DECOMPOSITION_H */
