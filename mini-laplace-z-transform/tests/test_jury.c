#include <stdio.h>
#include <math.h>
#include "laplace_z_transform_core.h"
#include "stability_analysis.h"

int main(void) {
    printf("Jury test start\n"); fflush(stdout);
    double jpoly1[] = {0.06, -0.5, 1.0};
    printf("P(1) = %.6f\n", stability_jury_eval_P_at_1(jpoly1, 2)); fflush(stdout);
    printf("P(-1) check = %.6f\n", stability_jury_eval_P_at_minus_1(jpoly1, 2)); fflush(stdout);
    printf("calling jury_is_stable...\n"); fflush(stdout);
    int r = stability_jury_is_stable(jpoly1, 2);
    printf("result = %d\n", r); fflush(stdout);
    return 0;
}
