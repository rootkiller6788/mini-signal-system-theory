#include "signal_basis.h"
#include "signal_correlation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(void)
{
    printf("=== test_signal_correlation ===\n");
    signal_ct_t *x = signal_ct_alloc(0.0, 1.0, 0.01);
    assert(x);
    signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
    signal_ct_fill_sinusoid(x, &p);
    signal_ct_t *rxx = signal_ct_alloc(-0.5, 0.5, 0.01);
    assert(rxx);
    assert(signal_ct_autocorrelation(x, rxx) == 0);
    double r0 = signal_ct_energy(x);
    assert(fabs(rxx->samples[rxx->num_samples/2] - r0) < r0 * 0.1);
    signal_ct_t *mc = signal_ct_alloc(-1.0, 2.0, 0.01);
    assert(signal_ct_matched_filter(x, x, mc) == 0);
    double lag = signal_ct_find_best_lag(mc);
    printf("Best lag: %.4f\n", lag);
    assert(lag >= -2.0 && lag <= 2.0);
    double dft_in[4] = {1.0, 0.0, -1.0, 0.0};
    double Xr[4], Xi[4], xr[4], xi[4];
    assert(signal_dft_forward(dft_in, NULL, 4, Xr, Xi) == 0);
    assert(signal_dft_inverse(Xr, Xi, 4, xr, xi) == 0);
    for (int i = 0; i < 4; i++)
        assert(fabs(xr[i] - dft_in[i]) < 1e-10);
    signal_ct_free(x); signal_ct_free(rxx); signal_ct_free(mc);
    printf("All correlation tests passed\n");
    return 0;
}
