#include "signal_basis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

int main(void)
{
    printf("=== test_signal_basis ===\n");

    TEST("alloc valid"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        if (s && s->num_samples > 0) PASS(); else FAIL("alloc");
        signal_ct_free(s);
    }
    TEST("alloc null checks"); {
        if (signal_ct_alloc(1.0, 0.0, 0.1) == NULL &&
            signal_ct_alloc(0.0, 1.0, 0.0) == NULL &&
            signal_ct_alloc(0.0, 1.0, -0.1) == NULL) PASS(); else FAIL("null");
    }
    TEST("unit step"); {
        signal_ct_t *s = signal_ct_alloc(-0.5, 1.0, 0.01);
        signal_ct_fill_unit_step(s, 0.0);
        if (fabs(s->samples[0]) < 1e-10 && fabs(s->samples[s->num_samples-1] - 1.0) < 1e-10) PASS(); else FAIL("step");
        signal_ct_free(s);
    }
    TEST("sinusoid"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 0.1, 0.0001);
        signal_sinusoid_params_t p = {1.0, 50.0, 0.0};
        signal_ct_fill_sinusoid(s, &p);
        if (fabs(s->samples[0] - 1.0) < 1e-10) PASS(); else FAIL("sin");
        signal_ct_free(s);
    }
    TEST("energy"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(s, &p);
        if (signal_ct_energy(s) > 0.0) PASS(); else FAIL("energy");
        signal_ct_free(s);
    }
    TEST("power"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(s, &p);
        if (fabs(signal_ct_power(s) - 0.5) < 0.05) PASS(); else FAIL("power");
        signal_ct_free(s);
    }
    TEST("RMS"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(s, &p);
        if (fabs(signal_ct_rms(s) - 0.7071) < 0.05) PASS(); else FAIL("rms");
        signal_ct_free(s);
    }
    TEST("classify POWER"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(s, &p);
        if (signal_ct_classify(s) == SIGNAL_CLASS_POWER) PASS(); else FAIL("class");
        signal_ct_free(s);
    }
    TEST("classify known"); {
        signal_ct_t *s = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t sp = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(s, &sp);
        int c = signal_ct_classify(s);
        if (c == SIGNAL_CLASS_POWER || c == SIGNAL_CLASS_ENERGY) PASS(); else FAIL("class2");
        signal_ct_free(s);
    }
    TEST("inner product"); {
        signal_ct_t *a = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_ct_t *b = signal_ct_alloc(0.0, 1.0, 0.001);
        signal_sinusoid_params_t p = {1.0, 10.0, 0.0};
        signal_ct_fill_sinusoid(a, &p);
        signal_ct_fill_sinusoid(b, &p);
        double ip = signal_ct_inner_product(a, b);
        double E = signal_ct_energy(a);
        if (fabs(ip - E) < 0.01) PASS(); else FAIL("ip");
        signal_ct_free(a); signal_ct_free(b);
    }
    printf("\nResults: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
