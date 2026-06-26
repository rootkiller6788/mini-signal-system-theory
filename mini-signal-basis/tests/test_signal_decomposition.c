#include "signal_basis.h"
#include "signal_decomposition.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(void)
{
    printf("=== test_signal_decomposition ===\n");
    basis_set_t *bs = basis_set_alloc(16);
    assert(bs);
    assert(basis_set_create_fourier(bs, 1.0, 3, 0.001) == 0);
    assert(bs->count == 7);
    signal_ct_t *x = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_sinusoid_params_t p = {1.0, 2.0, 0.0};
    signal_ct_fill_sinusoid(x, &p);
    double coeffs[7];
    assert(signal_project_onto_basis(x, bs, coeffs) == 0);
    signal_approximation_t approx = {0};
    assert(signal_truncated_approximation(x, bs, 7, &approx) == 0);
    assert(approx.energy_ratio > 0.9);
    signal_approximation_free(&approx);
    signal_ct_free(x);
    basis_set_free(bs);
    printf("All decomposition tests passed\n");
    return 0;
}
