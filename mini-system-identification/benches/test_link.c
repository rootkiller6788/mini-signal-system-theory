#include <stdio.h>
#include "../include/sysid_models.h"
int main(void) {
    printf("start\n");
    sysid_arx_t model;
    int ret = sysid_arx_alloc(&model, 2, 2, 1, 0.001);
    printf("alloc ret=%d\n", ret);
    sysid_arx_free(&model);
    printf("done\n");
    return 0;
}
