#include "signal_basis.h"
#include "signal_ops.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main(void)
{
    printf("=== test_signal_ops ===\n");
    signal_ct_t *a = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_ct_t *b = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_ct_t *r = signal_ct_alloc(0.0, 1.0, 0.001);
    assert(a && b && r);
    signal_sinusoid_params_t p = {1.0, 5.0, 0.0};
    signal_ct_fill_sinusoid(a, &p);
    p.amplitude = 0.5;
    signal_ct_fill_sinusoid(b, &p);
    assert(signal_ct_add(a, b, r) == 0);
    assert(signal_ct_multiply(a, b, r) == 0);
    assert(signal_ct_amplitude_scale(a, 2.0, r) == 0);
    double peak = signal_ct_peak(r);
    assert(fabs(peak - 2.0) < 0.1);
    assert(signal_ct_time_shift(a, 0.1, r) == 0);
    signal_ct_t *e = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_ct_t *o = signal_ct_alloc(0.0, 1.0, 0.001);
    signal_ct_even_odd_decompose(a, e, o);
    signal_ct_free(a); signal_ct_free(b); signal_ct_free(r);
    signal_ct_free(e); signal_ct_free(o);
    printf("All ops tests passed\n");
    return 0;
}
