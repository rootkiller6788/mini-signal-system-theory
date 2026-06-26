/**
 * @file demo_bode_analysis.c
 * @brief L6 Demo: Bode Plot Analysis and Filter Characterization
 *
 * Interactive demonstration of frequency response analysis:
 *   - Bode magnitude/phase plots for lowpass/bandpass/notch filters
 *   - Resonance detection and Q factor computation
 *   - Filter type classification (lowpass/highpass/bandpass/bandstop)
 *   - Roll-off rate measurement (dB/decade and dB/octave)
 *   - Bilinear prewarping for digital filter design
 *
 * Course Mapping:
 *   MIT 6.003 Ch.6 (Frequency Response), Stanford EE102A Ch.5
 *   Berkeley EE123 Ch.5
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/system_defs.h"
#include "../include/transfer_function.h"
#include "../include/frequency_response.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void print_banner(const char *title) {
    printf("\n========================================\n  %s\n========================================\n", title);
}
static void print_tf(const ct_transfer_function_t *tf, const char *name) {
    printf("  %s: H(s) = (", name);
    for (int i = tf->num_order; i >= 0; i--) {
        if (i < tf->num_order && tf->num_coeffs[i] >= 0) printf("+");
        printf("%.4g", tf->num_coeffs[i]);
        if (i > 0) printf("s^%d", i);
    }
    printf(") / (");
    for (int i = tf->den_order; i >= 0; i--) {
        if (i < tf->den_order && tf->den_coeffs[i] >= 0) printf("+");
        printf("%.4g", tf->den_coeffs[i]);
        if (i > 0) printf("s^%d", i);
    }
    printf(")\n");
}

int main(void) {
    printf("=== Bode Plot Analysis and Filter Characterization ===\n");
    double test_freqs[6], w, f_res, Q, f_low, f_high;
    double complex H;
    int i;

    /* Demo 1: First-Order RC Lowpass, tau=1ms, fc=159 Hz */
    print_banner("Demo 1: First-Order RC Lowpass Filter");
    ct_transfer_function_t tf_lp = ct_tf_alloc(0, 1);
    tf_lp.num_coeffs[0] = 1.0;
    tf_lp.den_coeffs[0] = 1.0; tf_lp.den_coeffs[1] = 1.0e-3;
    print_tf(&tf_lp, "RC Lowpass");

    frequency_response_t fr_lp = freq_resp_alloc(100, 1.0, 100000.0, 1);
    ct_frequency_response(&tf_lp, &fr_lp);

    printf("  DC gain: %.4f (%.2f dB)\n", compute_dc_gain_ct(&tf_lp), 20.0*log10(fabs(compute_dc_gain_ct(&tf_lp))+1e-30));
    double fc = find_3db_cutoff_frequency(&fr_lp);
    printf("  -3 dB cutoff: %.2f Hz (expected fc=1/(2*pi*RC)=%.2f Hz)\n", fc, 1.0/(2.0*M_PI*1e-3));
    double rolloff = compute_rolloff_db_per_decade(&fr_lp, fc*10, fc*100);
    printf("  Roll-off: %.2f dB/decade (expected -20, 1st-order)\n", rolloff);
    printf("  Filter type: %s\n", is_lowpass(&fr_lp) ? "LOWPASS" : (is_highpass(&fr_lp) ? "HIGHPASS" : (is_bandpass(&fr_lp) ? "BANDPASS" : "OTHER")));

    test_freqs[0]=1.0; test_freqs[1]=10.0; test_freqs[2]=100.0;
    test_freqs[3]=fc; test_freqs[4]=1000.0; test_freqs[5]=10000.0;
    for (i=0;i<6;i++) { w=2.0*M_PI*test_freqs[i]; H=ct_tf_evaluate(&tf_lp,0.0+w*I);
        printf("  f=%7.1f Hz: |H|=%.4f (%6.2f dB), arg=%7.2f deg\n", test_freqs[i], cabs(H), 20.0*log10(cabs(H)+1e-30), carg(H)*180.0/M_PI); }
    freq_resp_free(&fr_lp); ct_tf_free(&tf_lp);

    /* Demo 2: RLC Bandpass, R=100, L=10mH, C=1uF */
    print_banner("Demo 2: RLC Bandpass Filter");
    double R=100.0, L=10e-3, Cv=1e-6, w0_sq=1.0/(L*Cv), w0=sqrt(w0_sq);
    ct_transfer_function_t tf_bp = ct_tf_alloc(1, 2);
    tf_bp.num_coeffs[1]=R/L; tf_bp.num_coeffs[0]=0.0;
    tf_bp.den_coeffs[0]=w0_sq; tf_bp.den_coeffs[1]=R/L; tf_bp.den_coeffs[2]=1.0;
    print_tf(&tf_bp, "RLC Bandpass");

    frequency_response_t fr_bp = freq_resp_alloc(100, 10.0, 100000.0, 1);
    ct_frequency_response(&tf_bp, &fr_bp);
    f_res=find_resonance_frequency(&fr_bp);
    printf("  Resonance: %.2f Hz (expected f0=%.2f Hz)\n", f_res, w0/(2.0*M_PI));
    Q=compute_q_factor(&fr_bp);
    printf("  Q factor: %.2f (theoretical Q=%.2f)\n", Q, w0*L/R);
    if (find_half_power_points(&fr_bp, &f_low, &f_high)==0)
        printf("  -3 dB points: %.2f - %.2f Hz, BW=%.2f Hz\n", f_low, f_high, f_high-f_low);
    printf("  Filter type: %s\n", is_bandpass(&fr_bp) ? "BANDPASS" : "OTHER");
    printf("  Passband ripple [%.0f-%.0f Hz]: %.2f dB\n", f_res*0.5, f_res*2.0, compute_passband_ripple_db(&fr_bp, f_res*0.5, f_res*2.0));
    freq_resp_free(&fr_bp); ct_tf_free(&tf_bp);

    /* Demo 3: Bilinear Prewarping */
    print_banner("Demo 3: Bilinear Frequency Prewarping");
    double Fs=8000.0, analog_fc=1000.0;
    printf("  Analog cutoff: %.2f Hz, Fs: %.2f Hz\n", analog_fc, Fs);
    printf("  Prewarped digital fc: %.2f Hz (ratio %.4f)\n", bilinear_prewarp_frequency(analog_fc,Fs), bilinear_prewarp_frequency(analog_fc,Fs)/analog_fc);
    double ta[5]={100.0,500.0,1000.0,2000.0,3500.0};
    for(i=0;i<5;i++) { double dg=bilinear_prewarp_frequency(ta[i],Fs); printf("  f_a=%6.1f -> f_d=%6.1f Hz (ratio=%.4f)\n", ta[i], dg, dg/ta[i]); }

    /* Demo 4: 60 Hz Notch Filter */
    print_banner("Demo 4: 60 Hz Notch Filter (Bandstop)");
    double f_notch=60.0, w_notch=2.0*M_PI*f_notch, Q_notch=10.0;
    ct_transfer_function_t tf_notch = ct_tf_alloc(2, 2);
    tf_notch.num_coeffs[0]=w_notch*w_notch; tf_notch.num_coeffs[1]=0.0; tf_notch.num_coeffs[2]=1.0;
    tf_notch.den_coeffs[0]=w_notch*w_notch; tf_notch.den_coeffs[1]=w_notch/Q_notch; tf_notch.den_coeffs[2]=1.0;
    print_tf(&tf_notch, "60Hz Notch");

    frequency_response_t fr_notch = freq_resp_alloc(100, 30.0, 120.0, 0);
    ct_frequency_response(&tf_notch, &fr_notch);
    printf("  Stopband attenuation [58-62 Hz]: %.2f dB\n", compute_stopband_attenuation_db(&fr_notch, 58.0, 62.0));
    printf("  Filter type: %s\n", is_bandstop(&fr_notch) ? "BANDSTOP" : "OTHER");
    printf("  Notch Profile (every 10th point):\n");
    for(i=0;i<fr_notch.num_points;i+=10) printf("    f=%8.2f Hz, |H|=%.6f, dB=%8.2f\n", fr_notch.points[i].frequency, fr_notch.points[i].magnitude, fr_notch.points[i].magnitude_db);
    freq_resp_free(&fr_notch); ct_tf_free(&tf_notch);

    /* Demo 5: Pole-Zero Map (3rd-order Butterworth) */
    print_banner("Demo 5: Pole-Zero Map and Stability");
    ct_transfer_function_t tf_butter = ct_tf_alloc(0, 3);
    tf_butter.num_coeffs[0]=1.0;
    tf_butter.den_coeffs[0]=1.0; tf_butter.den_coeffs[1]=2.0; tf_butter.den_coeffs[2]=2.0; tf_butter.den_coeffs[3]=1.0;
    pole_zero_collection_t pz = pz_collection_alloc(3, 0);
    ct_tf_pole_zero_analysis(&tf_butter, &pz);
    printf("  3rd-order Butterworth lowpass\n");
    printf("  %-4s %-20s %-8s %-8s %s\n", "Idx","Pole","Damping","Wn","Stable?");
    for(i=0;i<pz.num_poles;i++) printf("  %-4d %+.4f%+.4fj  %-8.4f %-8.4f %s\n", i, creal(pz.poles[i].value), cimag(pz.poles[i].value), pz.poles[i].damping, pz.poles[i].wn, pz.is_stable?"YES":"NO");
    pz_collection_free(&pz); ct_tf_free(&tf_butter);

    printf("\n=== Bode Analysis Demo Complete ===\n");
    return 0;
}
