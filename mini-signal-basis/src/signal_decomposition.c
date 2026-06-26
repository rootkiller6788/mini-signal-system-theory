#include "signal_decomposition.h"
#include "signal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =================================================================
 * L5: Basis Set Management
 * ================================================================= */

basis_set_t *basis_set_alloc(int capacity)
{
    if (capacity <= 0) return NULL;
    basis_set_t *bs = (basis_set_t *)malloc(sizeof(basis_set_t));
    if (!bs) return NULL;
    bs->functions = (basis_function_t **)calloc((size_t)capacity, sizeof(basis_function_t *));
    if (!bs->functions) { free(bs); return NULL; }
    bs->count = 0;
    bs->capacity = capacity;
    return bs;
}

int basis_set_add(basis_set_t *bs, const signal_ct_t *phi)
{
    if (!bs || !phi || bs->count >= bs->capacity) return -1;
    signal_ct_t *copy = signal_ct_copy(phi);
    if (!copy) return -1;
    bs->functions[bs->count] = copy;
    bs->count++;
    return bs->count - 1;
}

void basis_set_free(basis_set_t *bs)
{
    if (!bs) return;
    for (int i = 0; i < bs->count; i++)
        signal_ct_free(bs->functions[i]);
    free(bs->functions);
    free(bs);
}

void signal_approximation_free(signal_approximation_t *approx)
{
    if (!approx) return;
    free(approx->coefficients);
    signal_ct_free(approx->approximation);
    approx->coefficients = NULL;
    approx->approximation = NULL;
}

/* =================================================================
 * L5: Gram-Schmidt Orthogonalization (Modified GS for stability)
 * u_k = v_k - sum_{i=1}^{k-1} <v_k, phi_i> * phi_i
 * phi_k = u_k / ||u_k||
 * Reference: O&W (1997) section 3.2.3
 * ================================================================= */

int signal_gram_schmidt(const signal_ct_t **vectors, int num_vectors, basis_set_t *orthonormal_basis)
{
    if (!vectors || num_vectors <= 0 || !orthonormal_basis) return -1;
    if (orthonormal_basis->capacity < num_vectors) return -1;
    for (int k = 0; k < num_vectors; k++) {
        signal_ct_t *uk = signal_ct_copy(vectors[k]);
        if (!uk) return -1;
        for (int i = 0; i < k; i++) {
            double proj = signal_ct_inner_product(vectors[k], orthonormal_basis->functions[i]);
            for (signal_index_t n = 0; n < uk->num_samples; n++)
                uk->samples[n] -= proj * orthonormal_basis->functions[i]->samples[n];
        }
        double norm = signal_ct_norm(uk);
        if (norm < 1e-15) { signal_ct_free(uk); return -1; }
        for (signal_index_t n = 0; n < uk->num_samples; n++)
            uk->samples[n] /= norm;
        if (basis_set_add(orthonormal_basis, uk) < 0) { signal_ct_free(uk); return -1; }
        signal_ct_free(uk);
    }
    return 0;
}

/* =================================================================
 * L5: Basis Projection
 * c_k = <x, phi_k>
 * ================================================================= */

int signal_project_onto_basis(const signal_ct_t *x, const basis_set_t *basis, double *coefficients)
{
    if (!x || !basis || !coefficients) return -1;
    for (int k = 0; k < basis->count; k++)
        coefficients[k] = signal_ct_inner_product(x, basis->functions[k]);
    return 0;
}

int signal_reconstruct_from_basis(const double *coefficients, int K, const basis_set_t *basis, signal_ct_t *reconstructed)
{
    if (!coefficients || !basis || !reconstructed || K > basis->count) return -1;
    for (signal_index_t n = 0; n < reconstructed->num_samples; n++)
        reconstructed->samples[n] = 0.0;
    for (int k = 0; k < K; k++) {
        signal_ct_t *phi = basis->functions[k];
        for (signal_index_t n = 0; n < reconstructed->num_samples; n++)
            reconstructed->samples[n] += coefficients[k] * phi->samples[n];
    }
    return 0;
}

int signal_truncated_approximation(const signal_ct_t *x, const basis_set_t *basis, int K, signal_approximation_t *result)
{
    if (!x || !basis || !result || K <= 0 || K > basis->count) return -1;
    result->num_terms = K;
    result->coefficients = (double *)malloc((size_t)K * sizeof(double));
    if (!result->coefficients) return -1;
    signal_project_onto_basis(x, basis, result->coefficients);
    result->approximation = signal_ct_alloc(x->t_start, x->t_end, x->dt);
    if (!result->approximation) { free(result->coefficients); return -1; }
    signal_reconstruct_from_basis(result->coefficients, K, basis, result->approximation);
    signal_ct_t *diff = signal_ct_alloc(x->t_start, x->t_end, x->dt);
    if (diff) {
        for (signal_index_t n = 0; n < diff->num_samples; n++)
            diff->samples[n] = x->samples[n] - result->approximation->samples[n];
        result->residual_energy = signal_ct_energy(diff);
        double orig_energy = signal_ct_energy(x);
        result->energy_ratio = (orig_energy > 1e-15) ? 1.0 - result->residual_energy / orig_energy : 1.0;
        signal_ct_free(diff);
    }
    return 0;
}

/* =================================================================
 * L5: Fourier Basis Construction
 * phi_0(t) = 1/sqrt(T)  (DC)
 * phi_{2k-1}(t) = sqrt(2/T) * cos(2*pi*k*f0*t), k=1..K
 * phi_{2k}(t) = sqrt(2/T) * sin(2*pi*k*f0*t), k=1..K
 * Reference: O&W (1997) section 3.2
 * ================================================================= */

int basis_set_create_fourier(basis_set_t *bs, double T, int num_harmonics, signal_time_t dt)
{
    if (!bs || T <= 0.0 || num_harmonics < 0 || dt <= 0.0) return -1;
    int num_funcs = 1 + 2 * num_harmonics;
    if (bs->capacity < num_funcs) return -1;
    double f0 = 1.0 / T;
    signal_ct_t *phi_dc = signal_ct_alloc(0.0, T, dt);
    if (!phi_dc) return -1;
    double amp_dc = 1.0 / sqrt(T);
    for (signal_index_t n = 0; n < phi_dc->num_samples; n++)
        phi_dc->samples[n] = amp_dc;
    basis_set_add(bs, phi_dc); signal_ct_free(phi_dc);
    double amp_harm = sqrt(2.0 / T);
    signal_sinusoid_params_t sparams;
    sparams.amplitude = amp_harm;
    for (int k = 1; k <= num_harmonics; k++) {
        sparams.frequency = (double)k * f0;
        sparams.phase = 0.0;
        signal_ct_t *phi_cos = signal_ct_alloc(0.0, T, dt);
        if (!phi_cos) return -1;
        signal_ct_fill_sinusoid(phi_cos, &sparams);
        basis_set_add(bs, phi_cos); signal_ct_free(phi_cos);
        sparams.phase = -M_PI / 2.0;
        signal_ct_t *phi_sin = signal_ct_alloc(0.0, T, dt);
        if (!phi_sin) return -1;
        signal_ct_fill_sinusoid(phi_sin, &sparams);
        basis_set_add(bs, phi_sin); signal_ct_free(phi_sin);
    }
    return 0;
}

/* =================================================================
 * L8: Haar Wavelet Basis Construction
 * Scaling: phi(t) = 1 on [0,1)
 * Mother wavelet: psi(t) = 1 on [0,0.5), -1 on [0.5,1)
 * Dilated/translated: psi_{j,k}(t) = 2^(j/2) * psi(2^j * t - k)
 * Reference: Strang & Nguyen (1996) Ch.1
 * ================================================================= */

int basis_set_create_haar(basis_set_t *bs, signal_time_t duration, signal_time_t dt, int max_level)
{
    if (!bs || duration <= 0.0 || dt <= 0.0 || max_level < 0) return -1;
    signal_index_t N = (signal_index_t)(duration / dt) + 1;
    double norm_factor = 1.0 / sqrt(duration);
    signal_ct_t *scaling = signal_ct_alloc(0.0, duration, dt);
    if (!scaling) return -1;
    for (signal_index_t n = 0; n < N; n++)
        scaling->samples[n] = norm_factor;
    basis_set_add(bs, scaling); signal_ct_free(scaling);
    for (int level = 0; level <= max_level; level++) {
        int num_shifts = 1 << level;
        double width = duration / (double)num_shifts;
        double amp = sqrt(1.0 / width);
        for (int shift = 0; shift < num_shifts; shift++) {
            signal_ct_t *psi = signal_ct_alloc(0.0, duration, dt);
            if (!psi) return -1;
            double start = (double)shift * width;
            double mid = start + width * 0.5;
            double end = start + width;
            for (signal_index_t n = 0; n < N; n++) {
                double t = psi->t_start + (double)n * psi->dt;
                if (t >= start && t < mid) psi->samples[n] = amp;
                else if (t >= mid && t < end) psi->samples[n] = -amp;
                else psi->samples[n] = 0.0;
            }
            basis_set_add(bs, psi); signal_ct_free(psi);
        }
    }
    return 0;
}

/* =================================================================
 * L6: Least-Squares Signal Approximation
 * Minimize ||x - sum c_k * v_k||^2 via normal equations G*c = b
 * G[i][j] = <v_i, v_j> (Gram matrix), b[i] = <x, v_i>
 * Solved via Cholesky decomposition (G is SPD for LI templates).
 * ================================================================= */

static int cholesky_decompose(double *A, int n)
{
    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) {
            double sum = A[i * n + j];
            for (int k = 0; k < j; k++)
                sum -= A[i * n + k] * A[j * n + k];
            if (i == j) {
                if (sum <= 0.0) return -1;
                A[j * n + j] = sqrt(sum);
            } else {
                A[i * n + j] = sum / A[j * n + j];
            }
        }
    }
    return 0;
}

static void cholesky_solve(const double *L, int n, const double *b, double *x)
{
    double *y = (double *)malloc((size_t)n * sizeof(double));
    if (!y) return;
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = 0; j < i; j++)
            sum -= L[i * n + j] * y[j];
        y[i] = sum / L[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < n; j++)
            sum -= L[j * n + i] * x[j];
        x[i] = sum / L[i * n + i];
    }
    free(y);
}

int signal_least_squares_approx(const signal_ct_t *x, const signal_ct_t **templates, int K, signal_approximation_t *result)
{
    if (!x || !templates || K <= 0 || !result) return -1;
    double *G = (double *)calloc((size_t)(K * K), sizeof(double));
    double *b = (double *)calloc((size_t)K, sizeof(double));
    if (!G || !b) { free(G); free(b); return -1; }
    for (int i = 0; i < K; i++) {
        b[i] = signal_ct_inner_product(x, templates[i]);
        for (int j = 0; j < K; j++)
            G[i * K + j] = signal_ct_inner_product(templates[i], templates[j]);
    }
    if (cholesky_decompose(G, K) != 0) { free(G); free(b); return -1; }
    result->num_terms = K;
    result->coefficients = (double *)calloc((size_t)K, sizeof(double));
    if (!result->coefficients) { free(G); free(b); return -1; }
    cholesky_solve(G, K, b, result->coefficients);
    result->approximation = signal_ct_alloc(x->t_start, x->t_end, x->dt);
    if (!result->approximation) { free(G); free(b); free(result->coefficients); return -1; }
    for (signal_index_t n = 0; n < result->approximation->num_samples; n++)
        result->approximation->samples[n] = 0.0;
    for (int k = 0; k < K; k++)
        for (signal_index_t n = 0; n < result->approximation->num_samples; n++)
            result->approximation->samples[n] += result->coefficients[k] * templates[k]->samples[n];
    signal_ct_t *diff = signal_ct_alloc(x->t_start, x->t_end, x->dt);
    if (diff) {
        for (signal_index_t n = 0; n < diff->num_samples; n++)
            diff->samples[n] = x->samples[n] - result->approximation->samples[n];
        result->residual_energy = signal_ct_energy(diff);
        double orig_energy = signal_ct_energy(x);
        result->energy_ratio = (orig_energy > 1e-15) ? 1.0 - result->residual_energy / orig_energy : 1.0;
        signal_ct_free(diff);
    }
    free(G); free(b);
    return 0;
}
