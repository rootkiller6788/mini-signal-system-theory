#include "signal_basis.h"
#include "signal_energy.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(void)
{
    printf("=== test_signal_energy ===\n");
    signal_ct_t *x = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
    signal_ct_fill_sinusoid(x, &p);
    parseval_result_t pr = signal_verify_parseval(x, 1e-6);
    assert(pr.verified);
    assert(pr.relative_error < 1e-6);
    double esd[1024], freqs[1024];
    assert(signal_energy_spectral_density(x, freqs, esd, 1024) == 0);
    double psd[1024];
    assert(signal_psd_periodogram(x, freqs, psd, 1024) == 0);
    double cf = signal_ct_crest_factor(x);
    assert(cf > 1.0 && cf < 2.0);
    signal_ct_free(x);
    printf("All energy tests passed\n");
    return 0;
}
