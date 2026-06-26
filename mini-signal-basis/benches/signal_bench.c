#include "signal_basis.h"
#include "signal_correlation.h"
#include <stdio.h><stdlib.h><time.h>
int main(void){
    printf("=== Signal Benchmarks ===
");
    signal_ct_t *x=signal_ct_alloc(0.0,1.0,0.0001);
    signal_sinusoid_params_t p={1.0,100.0,0.0};
    signal_ct_fill_sinusoid(x,&p);
    clock_t t0=clock();
    double E=signal_ct_energy(x);
    clock_t t1=clock();
    printf("Energy (N=%lld): %.6f in %.3fms
",
        (long long)x->num_samples,E,1000.0*(double)(t1-t0)/CLOCKS_PER_SEC);
    double dft_data[4096];
    int N_dft=4096;
    double *xr=(double*)malloc(N_dft*sizeof(double));
    double *Xr=(double*)malloc(N_dft*sizeof(double));
    double *Xi=(double*)malloc(N_dft*sizeof(double));
    t0=clock();
    signal_dft_forward(xr,NULL,N_dft,Xr,Xi);
    t1=clock();
    printf("DFT (N=%d): %.3fms O(N^2)
",N_dft,
        1000.0*(double)(t1-t0)/CLOCKS_PER_SEC);
    free(xr);free(Xr);free(Xi);
    signal_ct_free(x);
    return 0;
}
