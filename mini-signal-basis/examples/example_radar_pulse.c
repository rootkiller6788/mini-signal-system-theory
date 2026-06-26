/* Example: Radar Pulse Compression */
#include "signal_basis.h"
#include "signal_correlation.h"
#include "signal_applications.h"
#include <stdio.h><stdlib.h><math.h><assert.h>
int main(void){
    printf("=== Radar Pulse Compression ===
");
    double dt=1e-7, B=1e6, T=10e-6, f0=10e6;
    signal_ct_t *pulse=signal_ct_alloc(0.0,T,dt);
    assert(pulse);
    assert(signal_generate_chirp_pulse(pulse,f0,B,T,1.0)==0);
    printf("Chirp pulse: T=%.1fus B=%.1fMHz TB=%.0f
",
        T*1e6,B*1e-6,T*B);
    double E=signal_ct_energy(pulse);
    printf("Pulse energy: %.3e
",E);
    signal_ct_t *rx=signal_ct_alloc(0.0,2.0*T,dt);
    for(signal_index_t i=0;i<pulse->num_samples;i++)
        rx->samples[rx->num_samples/4+i]=pulse->samples[i]*0.5;
    radar_pulse_compression_t rpc;
    signal_radar_pulse_compression(rx,pulse,f0,&rpc);
    printf("Range profile: %d bins, res=%.3fm
",
        rpc.num_range_bins,rpc.range_resolution);
    printf("Peak at bin %d, power %.1fdB
",
        rpc.peak_bin,rpc.peak_power_dB);
    radar_pulse_compression_free(&rpc);
    signal_ct_free(pulse);signal_ct_free(rx);
    printf("Radar example complete.
");
    return 0;
}
