/* Example: ECG QRS Detection */
#include "signal_basis.h"
#include "signal_applications.h"
#include <stdio.h><stdlib.h><math.h>
int main(void){
    printf("=== ECG QRS Detection ===
");
    double fs=360.0, dt=1.0/fs, dur=3.0;
    signal_ct_t *ecg=signal_ct_alloc(0.0,dur,dt);
    assert(ecg);
    for(signal_index_t i=0;i<ecg->num_samples;i++){
        double t=(double)i*dt;
        double phase=fmod(t,0.8);
        if(phase<0.08){
            ecg->samples[i]=1.0*exp(-phase*50.0)*sin(2.0*M_PI*20.0*phase);
        } else { ecg->samples[i]=0.05*sin(2.0*M_PI*1.25*t); }
    }
    ecg_detection_t result;
    signal_detect_ecg_qrs(ecg,fs,&result);
    printf("Detected %d beats, mean HR=%.1f BPM
",
        result.num_beats,result.mean_heart_rate);
    assert(result.num_beats>=2);
    ecg_detection_free(&result);
    signal_ct_free(ecg);
    printf("ECG example complete.
");
    return 0;
}
