#include "signal_basis.h"
#include "signal_correlation.h"
#include "signal_applications.h"
#include <stdio.h><stdlib.h><math.h>
int main(void){
    printf("=== Audio Tone Detection (DTMF) ===
");
    double fs=8000.0, dt=1.0/fs;
    double dur=0.5;
    signal_ct_t *audio=signal_ct_alloc(0.0,dur,dt);
    signal_sinusoid_params_t p1={0.5,697.0,0.0};
    signal_sinusoid_params_t p2={0.5,1209.0,0.0};
    signal_ct_fill_sinusoid(audio,&p1);
    signal_ct_t *tmp=signal_ct_alloc(0.0,dur,dt);
    signal_ct_fill_sinusoid(tmp,&p2);
    for(signal_index_t i=0;i<audio->num_samples;i++) audio->samples[i]+=tmp->samples[i];
    signal_ct_free(tmp);
    dtmf_tone_t result;
    signal_detect_dtmf(audio,fs,&result);
    printf("Detected: row=%.0fHz col=%.0fHz digit=%c conf=%.3f
",
        result.row_freq,result.col_freq,result.digit,result.confidence);
    assert(result.digit==1);
    signal_ct_free(audio);
    printf("DTMF example complete.
");
    return 0;
}
