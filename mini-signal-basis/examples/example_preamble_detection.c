/* Example: Communication Preamble Detection */
#include "signal_basis.h"
#include "signal_correlation.h"
#include "signal_applications.h"
#include <stdio.h><stdlib.h><math.h><assert.h>
int main(void){
    printf("=== Preamble Detection ===
");
    double dt=0.001, dur=0.5;
    signal_ct_t *preamble=signal_ct_alloc(0.0,0.1,dt);
    assert(signal_generate_barker_pulse(preamble,7,0.01,1.0)==0);
    signal_ct_t *rx=signal_ct_alloc(0.0,dur,dt);
    double delay=0.2;
    for(signal_index_t i=0;i<rx->num_samples;i++){
        double t=(double)i*dt;
        double idx=(t-delay)/dt;
        if(idx>=0&&idx<preamble->num_samples)
            rx->samples[i]=preamble->samples[(signal_index_t)idx];
    }
    preamble_detection_t result;
    signal_detect_preamble(rx,preamble,5.0,&result);
    printf("Detected=%d time=%.3fs peak=%.3f SNR=%.1fdB
",
        result.detected,result.detection_time,result.correlation_peak,result.snr_estimate);
    assert(result.detected);
    signal_ct_free(preamble);signal_ct_free(rx);
    printf("Preamble example complete.
");
    return 0;
}
