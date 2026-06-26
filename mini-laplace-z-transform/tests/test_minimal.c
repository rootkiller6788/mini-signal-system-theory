#include <stdio.h>
#include "laplace_z_transform_core.h"
#include "stability_analysis.h"

int main(void) {
    printf("Starting minimal test...\n");
    double coeffs[] = {1.0, 2.0, 1.0};
    int result = stability_routh_is_stable(coeffs, 2);
    printf("Routh result: %d\n", result);
    return 0;
}
