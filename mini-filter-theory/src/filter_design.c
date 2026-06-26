/**
 * filter_design.c — Complete Filter Design Methods
 *
 * L5: Window function generation (10 window types)
 * L5: FIR design: window, frequency sampling, least-squares
 * L5: IIR design: bilinear transform, impulse invariance
 * L5: Order estimation for all filter families
 * L5: Coefficient quantization effects
 *
 * Reference:
 *   Harris, Proc. IEEE, 66(1):51-83, 1978 (Windows)
 *   Parks & Burrus, Digital Filter Design, Wiley, 1987
 *   Antoniou, Digital Filters, McGraw-Hill, 1993
 *   Oppenheim & Schafer (2010) Ch. 7
 */

#include "filter_design.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * L5: Window Function Generation
 * ================================================================== */

/**
 * window_generate — L5: all canonical window functions
 *
 * Generates symmetric windows of length N.
 * All windows are normalized to have maximum value = 1.0.
 *
 * Key design insight (Harris 1978):
 *   The window's frequency-domain characteristics determine the filter's
 *   transition bandwidth and stopband attenuation:
 *   - Main lobe width determines transition bandwidth
 *   - Side lobe level determines minimum stopband attenuation
 *   - Side lobe roll-off determines how attenuation improves with frequency
 *
 * Complexity: O(N) for all window types
 */
int window_generate(double *win, const window_params_t *params) {
    if (!win || !params || params->length < 1) return FILTER_ERR_NULL_PTR;

    int N = params->length;
    int n;
    double M = (N - 1) / 2.0;

    switch (params->type) {
    case WINDOW_RECTANGULAR:
        /* w[n] = 1.0 */
        for (n = 0; n < N; n++) win[n] = 1.0;
        break;

    case WINDOW_HANN:
        /* w[n] = 0.5 * (1 - cos(2*pi*n/(N-1))) */
        for (n = 0; n < N; n++) {
            win[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / (N - 1)));
        }
        break;

    case WINDOW_HAMMING:
        /* w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1))
         * Optimized to minimize the first side lobe level (-42.7 dB). */
        for (n = 0; n < N; n++) {
            win[n] = 0.54 - 0.46 * cos(2.0 * M_PI * n / (N - 1));
        }
        break;

    case WINDOW_BLACKMAN:
        /* w[n] = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1))
         * Three-term Blackman: side lobes -58 dB, 18 dB/octave roll-off. */
        for (n = 0; n < N; n++) {
            win[n] = 0.42
                     - 0.50 * cos(2.0 * M_PI * n / (N - 1))
                     + 0.08 * cos(4.0 * M_PI * n / (N - 1));
        }
        break;

    case WINDOW_BLACKMAN_HARRIS:
        /* 4-term Blackman-Harris: side lobes -92 dB.
         * Coefficients from Harris (1978), "minimum 4-term" window:
         * a0=0.35875, a1=0.48829, a2=0.14128, a3=0.01168 */
        for (n = 0; n < N; n++) {
            win[n] = 0.35875
                     - 0.48829 * cos(2.0 * M_PI * n / (N - 1))
                     + 0.14128 * cos(4.0 * M_PI * n / (N - 1))
                     - 0.01168 * cos(6.0 * M_PI * n / (N - 1));
        }
        break;

    case WINDOW_KAISER:
        /* w[n] = I0(beta*sqrt(1 - (2n/(N-1) - 1)^2)) / I0(beta)
         * where I0 is the zero-order modified Bessel function.
         *
         * The Kaiser window is near-optimal: for a given side lobe level,
         * it maximizes the main lobe energy concentration (DPSS sequence). */
        {
            double beta = params->alpha;
            if (beta < 0.0) beta = 5.0;

            /* Compute I0(beta) via series expansion */
            double I0_beta = 1.0;
            double term = 1.0;
            int k;
            for (k = 1; k < 50; k++) {
                term *= (beta / 2.0) * (beta / 2.0) / (k * k);
                I0_beta += term;
                if (term < 1e-15) break;
            }

            for (n = 0; n < N; n++) {
                double x = (2.0 * n / (N - 1) - 1.0);
                double arg = beta * sqrt(1.0 - x * x);
                /* Compute I0(arg) */
                double I0_arg = 1.0;
                term = 1.0;
                for (k = 1; k < 50; k++) {
                    term *= (arg / 2.0) * (arg / 2.0) / (k * k);
                    I0_arg += term;
                    if (term < 1e-15) break;
                }
                win[n] = I0_arg / I0_beta;
            }
        }
        break;

    case WINDOW_FLATTOP:
        /* 5-term flat-top: a0=1, a1=1.93, a2=1.29, a3=0.388, a4=0.028
         * Minimal scalloping loss (< 0.01 dB) for precision measurements.
         * Normalized to unity maximum. */
        for (n = 0; n < N; n++) {
            double a0 = 0.21557895, a1 = 0.41663158,
                   a2 = 0.277263158, a3 = 0.083578947, a4 = 0.006947368;
            win[n] = a0
                     - a1 * cos(2.0 * M_PI * n / (N - 1))
                     + a2 * cos(4.0 * M_PI * n / (N - 1))
                     - a3 * cos(6.0 * M_PI * n / (N - 1))
                     + a4 * cos(8.0 * M_PI * n / (N - 1));
        }
        break;

    case WINDOW_BARTLETT:
        /* w[n] = 1 - |2n/(N-1) - 1| (triangular window)
         * Simple, continuous, side lobes -26.5 dB. */
        for (n = 0; n < N; n++) {
            win[n] = 1.0 - fabs(2.0 * n / (N - 1) - 1.0);
        }
        break;

    case WINDOW_NUTTALL:
        /* 4-term Nuttall: continuous first derivative at boundaries.
         * a0=0.355768, a1=0.487396, a2=0.144232, a3=0.012604 */
        for (n = 0; n < N; n++) {
            win[n] = 0.355768
                     - 0.487396 * cos(2.0 * M_PI * n / (N - 1))
                     + 0.144232 * cos(4.0 * M_PI * n / (N - 1))
                     - 0.012604 * cos(6.0 * M_PI * n / (N - 1));
        }
        break;

    case WINDOW_TUKEY:
        /* w[n] = 1 for |n-M| <= alpha*N/2
         *       = 0.5*(1+cos(pi*(|n-M|-alpha*N/2)/((1-alpha)*N/2))) elsewhere
         * alpha=0 gives rectangular, alpha=1 gives Hann. */
        {
            double alpha = params->alpha;
            if (alpha < 0.0) alpha = 0.5;
            if (alpha > 1.0) alpha = 1.0;
            double limit = alpha * M;

            for (n = 0; n < N; n++) {
                double d = fabs(n - M) - limit;
                if (d <= 0.0) {
                    win[n] = 1.0;
                } else {
                    win[n] = 0.5 * (1.0 + cos(M_PI * d / (M - limit)));
                }
            }
        }
        break;

    default:
        return FILTER_ERR_NOT_IMPL;
    }

    return FILTER_OK;
}

/* ==================================================================
 * Ideal Impulse Responses
 * ================================================================== */

/**
 * ideal_lowpass_impulse — L5: sinc-based ideal LP
 *
 * h[n] = (wc/pi) * sinc(wc*(n-M)/pi) for n = 0..N-1
 *      = sin(wc*(n-M)) / (pi*(n-M)) for n != M
 *      = wc/pi for n = M (by L'Hopital's rule)
 *
 * M = (N-1)/2 is the filter center (group delay).
 * wc is normalized cutoff in rad/sample (0 to pi).
 *
 * Complexity: O(N)
 */
void ideal_lowpass_impulse(double *h, int N, double wc) {
    if (!h || N < 1) return;
    double M = (N - 1) / 2.0;
    int n;
    for (n = 0; n < N; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h[n] = wc / M_PI;
        } else {
            h[n] = sin(wc * t) / (M_PI * t);
        }
    }
}

/**
 * ideal_bandpass_impulse — L5: difference of two sincs
 *
 * h_BP[n] = h_LP(w2)[n] - h_LP(w1)[n]
 *         = (w2*sinc(w2*(n-M)) - w1*sinc(w1*(n-M))) / pi
 *
 * This follows from the linearity of the inverse DTFT:
 * H_BP(w) = rect(w1,w2) = rect(w2) - rect(w1).
 */
void ideal_bandpass_impulse(double *h, int N, double w1, double w2) {
    if (!h || N < 1 || w1 >= w2) return;
    double M = (N - 1) / 2.0;
    int n;
    for (n = 0; n < N; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h[n] = (w2 - w1) / M_PI;
        } else {
            h[n] = (sin(w2 * t) - sin(w1 * t)) / (M_PI * t);
        }
    }
}

/**
 * ideal_highpass_impulse — L5: delta minus LP
 *
 * h_HP[n] = delta[n-M] - h_LP[n]
 *         = delta[n-M] - wc/pi * sinc(wc*(n-M)/pi)
 *
 * where delta is the Kronecker delta (1 at n=M, 0 otherwise).
 */
void ideal_highpass_impulse(double *h, int N, double wc) {
    if (!h || N < 1) return;
    double M = (N - 1) / 2.0;
    int n;
    for (n = 0; n < N; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h[n] = 1.0 - wc / M_PI;
        } else {
            h[n] = -sin(wc * t) / (M_PI * t);
        }
    }
}

/* ==================================================================
 * FIR Design by Windowing
 * ================================================================== */

fir_filter_t *fir_design_by_window(const filter_spec_t *spec, int N,
                                    window_type_t win_type, double win_alpha) {
    if (!spec || !filter_spec_is_valid(spec) || N < 3) return NULL;
    if (spec->sample_rate <= 0.0) return NULL;

    double fs = spec->sample_rate;
    double wc1 = 2.0 * M_PI * spec->fc1 / fs; /* Normalized cutoff in rad/sample */

    double *h_ideal = (double *)malloc(N * sizeof(double));
    if (!h_ideal) return NULL;

    switch (spec->type) {
    case FILTER_LOWPASS:
        ideal_lowpass_impulse(h_ideal, N, wc1);
        break;
    case FILTER_HIGHPASS:
        ideal_highpass_impulse(h_ideal, N, wc1);
        break;
    case FILTER_BANDPASS:
        {
            double wc2 = 2.0 * M_PI * spec->fc2 / fs;
            ideal_bandpass_impulse(h_ideal, N, wc1, wc2);
        }
        break;
    case FILTER_BANDSTOP:
        {
            double wc2 = 2.0 * M_PI * spec->fc2 / fs;
            /* Bandstop = allpass - bandpass */
            ideal_lowpass_impulse(h_ideal, N, M_PI * 0.999); /* "Allpass" approx */
            int n;
            double M_val = (N - 1) / 2.0;
            for (n = 0; n < N; n++) {
                double t = n - M_val;
                double bp_val;
                if (fabs(t) < 1e-12) bp_val = (wc2 - wc1) / M_PI;
                else bp_val = (sin(wc2 * t) - sin(wc1 * t)) / (M_PI * t);
                h_ideal[n] -= bp_val;
            }
        }
        break;
    default:
        ideal_lowpass_impulse(h_ideal, N, wc1);
        break;
    }

    /* Generate window */
    double *window = (double *)malloc(N * sizeof(double));
    if (!window) { free(h_ideal); return NULL; }
    window_params_t wp = {win_type, win_alpha, N};
    window_generate(window, &wp);

    /* Apply window */
    double *coeff = (double *)malloc(N * sizeof(double));
    if (!coeff) { free(h_ideal); free(window); return NULL; }
    int n;
    for (n = 0; n < N; n++) coeff[n] = h_ideal[n] * window[n];

    fir_filter_t *fir = fir_filter_alloc(coeff, N);
    free(h_ideal); free(window); free(coeff);
    return fir;
}

/* ==================================================================
 * FIR Order Estimation
 * ================================================================== */

/**
 * fir_estimate_length — L5: Kaiser's formula
 *
 * N ≈ (-20*log10(sqrt(dp*ds)) - 13) / (14.6 * delta_w / (2*pi))
 *   + 1
 *
 * Simplified (Bellanger):
 *   N ≈ -10*log10(dp*ds) / (14.6 * delta_f) + 1
 *
 * where delta_f = transition bandwidth / fs (normalized)
 *
 * For typical values: -60 dB stopband, 0.1 dB ripple → N ≈ 5.5 / delta_f
 */
int fir_estimate_length(double transition_bw, double atten_db) {
    if (transition_bw <= 0.0 || transition_bw >= M_PI) return -1;
    /* Using attenuation-based estimation */
    double N = (atten_db - 8.0) / (2.285 * transition_bw) + 1.0;
    if (N < 3) N = 3;
    return (int)ceil(N);
}

/**
 * kaiser_beta_from_attenuation — L5
 *
 * Empirical formula from Kaiser (1974):
 *   beta = 0.1102*(As - 8.7) for As > 50 dB
 *   beta = 0.5842*(As-21)^0.4 + 0.07886*(As-21) for 21 <= As <= 50
 *   beta = 0 for As < 21
 */
double kaiser_beta_from_attenuation(double atten_db) {
    if (atten_db > 50.0) {
        return 0.1102 * (atten_db - 8.7);
    } else if (atten_db >= 21.0) {
        return 0.5842 * pow(atten_db - 21.0, 0.4)
               + 0.07886 * (atten_db - 21.0);
    } else {
        return 0.0;
    }
}

/**
 * kaiser_mainlobe_width — L5 approximate Kaiser window main lobe width
 *
 * The main lobe width of a Kaiser window in normalized frequency
 * (rad/sample) as a function of beta and window length N.
 *
 * Approximate formula:
 *   MLW ≈ 2*π * sqrt(beta^2 + π^2) / (N * β)   for β > 0
 *   MLW ≈ 4*π/N                                 for β = 0 (Rectangular)
 *
 * The main lobe width determines the transition bandwidth of the
 * resulting FIR filter: transition_bw ≈ MLW / (2*π) * fs.
 *
 * This is essential for predicting FIR filter performance before
 * actual design — allows quick estimation of required filter length.
 *
 * Reference: Kaiser (1974), "Nonrecursive Digital Filter Design
 *             Using the I0-sinh Window Function"
 *             Harris, Proc. IEEE, 66(1):51-83, 1978
 * Complexity: O(1)
 */
double kaiser_mainlobe_width(double beta, int N) {
    if (N < 2) return -1.0;
    if (beta < 0.0) beta = 0.0;

    if (beta < 0.001) {
        /* Rectangular window: main lobe width = 4*π/N */
        return 4.0 * M_PI / (double)N;
    }

    /* Kaiser window: the frequency response zero-crossing determines
     * the main lobe width.
     *
     * Approximate formula from Kaiser (1974):
     *   Main lobe width ≈ 2*sqrt(β^2 + π^2) / N  (normalized frequency)
     */
    double mlw = 2.0 * sqrt(beta * beta + M_PI * M_PI) / (double)N;

    /* Correction for the actual response: the Kaiser window has
     * slightly wider main lobe than the simple formula predicts.
     * Empirical correction factor ≈ 1.05 */
    mlw *= 1.05;

    return mlw;
}

/* ==================================================================
 * Frequency Sampling FIR Design
 * ================================================================== */

/**
 * fir_design_freq_sampling — L5
 *
 * Specifies desired frequency response at N equally spaced points
 * (DFT bins). The impulse response is the inverse DFT of these samples.
 *
 * For linear-phase Type I (N odd), only (N+1)/2 magnitude samples
 * are independent. The phase is constrained by the linear-phase condition.
 *
 * Complexity: O(N^2) for direct IDFT (O(N log N) with FFT)
 */
fir_filter_t *fir_design_freq_sampling(const double *desired_mag, int N,
                                        int is_linear_phase) {
    if (!desired_mag || N < 3) return NULL;

    double *h = (double *)malloc(N * sizeof(double));
    if (!h) return NULL;

    int k, n;
    int M = (N - 1) / 2;
    (void)M; /* Used in phase computation */

    /* IDFT: h[n] = (1/N) * sum_{k=0}^{N-1} H[k] * e^{j*2*pi*k*n/N}
     * For linear phase: H[k] = |H[k]| * e^{-j*pi*k*(N-1)/N} */
    for (n = 0; n < N; n++) {
        double sum = 0.0;

        /* DC bin (k=0) */
        sum += desired_mag[0];

        /* Positive frequencies k=1..(N-1)//2 */
        int K = (N - 1) / 2;
        for (k = 1; k <= K; k++) {
            double mag;
            if (k < N / 2) {
                mag = desired_mag[k];
            } else {
                mag = desired_mag[N - k]; /* Symmetry */
            }

            double phase;
            if (is_linear_phase) {
                phase = -M_PI * k * (N - 1) / N;
            } else {
                phase = 0.0;
            }

            double angle = 2.0 * M_PI * k * n / N + phase;
            sum += 2.0 * mag * cos(angle);
        }

        /* Nyquist bin (if N is even) */
        if (N % 2 == 0) {
            sum += desired_mag[N / 2] * cos(M_PI * n);
        }

        h[n] = sum / N;
    }

    fir_filter_t *fir = fir_filter_alloc(h, N);
    free(h);
    return fir;
}

/* ==================================================================
 * IIR Design via Bilinear Transform
 * ================================================================== */

iir_filter_t *iir_design_bilinear(const filter_spec_t *spec) {
    if (!spec || !filter_spec_is_valid(spec)) return NULL;
    if (spec->sample_rate <= 0.0) return NULL;

    /* Pre-warp cutoff frequencies */
    double fs = spec->sample_rate;
    double wc1_digital = 2.0 * M_PI * spec->fc1 / fs;
    double wc1_analog = bilinear_prewarp(wc1_digital, fs);

    /* Design analog prototype at pre-warped frequency */
    filter_spec_t analog_spec = *spec;
    analog_spec.sample_rate = 0.0; /* Mark as analog */
    analog_spec.fc1 = wc1_analog / (2.0 * M_PI);

    tf_continuous_t *analog_tf = analog_filter_design(&analog_spec);
    if (!analog_tf) return NULL;

    /* Apply bilinear transform */
    tf_discrete_t *digital_tf = bilinear_transform(analog_tf, fs);
    tf_continuous_free(analog_tf);
    if (!digital_tf) return NULL;

    /* Factor into biquads (simplified: single section) */
    biquad_section_t bq;
    if (digital_tf->den_len >= 3) {
        double a0 = digital_tf->den[0];
        if (fabs(a0) > 1e-15) {
            bq.b0 = digital_tf->num[0] / a0;
            bq.b1 = (digital_tf->num_len > 1) ? digital_tf->num[1] / a0 : 0.0;
            bq.b2 = (digital_tf->num_len > 2) ? digital_tf->num[2] / a0 : 0.0;
            bq.a1 = (digital_tf->den_len > 1) ? digital_tf->den[1] / a0 : 0.0;
            bq.a2 = (digital_tf->den_len > 2) ? digital_tf->den[2] / a0 : 0.0;
        } else {
            bq.b0 = 1.0; bq.b1 = 0.0; bq.b2 = 0.0;
            bq.a1 = 0.0; bq.a2 = 0.0;
        }
    } else {
        bq.b0 = 1.0; bq.b1 = 0.0; bq.b2 = 0.0;
        bq.a1 = 0.0; bq.a2 = 0.0;
    }
    bq.gain = digital_tf->gain;

    iir_filter_t *iir = iir_filter_alloc(&bq, 1, digital_tf->gain);
    tf_discrete_free(digital_tf);
    return iir;
}

/**
 * iir_design_impulse_invariance — L5 complete impulse invariance design
 *
 * Design steps:
 * 1. Design analog prototype satisfying the given specification
 * 2. Convert analog H(s) to digital H(z) via impulse_invariance()
 * 3. Factor digital TF into cascade biquads
 *
 * Important: Impulse invariance is only suitable for lowpass and
 * bandpass filters. For highpass and bandstop, aliasing from the
 * stopband folds into the passband, degrading the response.
 *
 * The method preserves the exact shape of the impulse response at
 * sample instants, making it ideal for time-domain applications.
 *
 * Reference: Oppenheim & Schafer (2010), Sec. 7.2.1
 * Complexity: O(N^2*iter + N^3) for pole finding + bilinear ops
 */
iir_filter_t *iir_design_impulse_invariance(const filter_spec_t *spec) {
    if (!spec || !filter_spec_is_valid(spec)) return NULL;
    if (spec->sample_rate <= 0.0) return NULL;

    /* Impulse invariance only works well for LP and BP */
    if (spec->type != FILTER_LOWPASS && spec->type != FILTER_BANDPASS)
        return NULL;

    double fs = spec->sample_rate;

    /* Step 1: Design analog prototype */
    filter_spec_t analog_spec = *spec;
    analog_spec.sample_rate = 0.0; /* Mark as analog */

    if (analog_spec.order <= 0) {
        double wp = 2.0 * M_PI * analog_spec.fc1 * 0.9;
        double ws = 2.0 * M_PI * analog_spec.fc1 * 1.5;
        if (analog_spec.passband_edge > 0.0) wp = 2.0 * M_PI * analog_spec.passband_edge;
        if (analog_spec.stopband_edge > 0.0) ws = 2.0 * M_PI * analog_spec.stopband_edge;
        analog_spec.order = butterworth_estimate_order(
            analog_spec.passband_ripple_db, analog_spec.stopband_atten_db, wp, ws);
        if (analog_spec.order < 1) analog_spec.order = 1;
        if (analog_spec.order > 10) return NULL;
    }

    tf_continuous_t *analog_tf = analog_filter_design(&analog_spec);
    if (!analog_tf) return NULL;

    /* Step 2: Impulse invariance conversion */
    tf_discrete_t *digital_tf = impulse_invariance(analog_tf, fs);
    tf_continuous_free(analog_tf);
    if (!digital_tf) return NULL;

    /* Step 3: Convert to cascade of biquads */
    int ns = tf_discrete_to_biquads(digital_tf, NULL, 0);
    int num_sections = (ns > 0) ? ns : 1;
    biquad_section_t *bqs =
        (biquad_section_t *)calloc(num_sections, sizeof(biquad_section_t));
    if (!bqs) { tf_discrete_free(digital_tf); return NULL; }

    if (digital_tf->den_len >= 2) {
        double a0 = digital_tf->den[0];
        if (fabs(a0) > 1e-15) {
            bqs[0].b0 = digital_tf->num[0] / a0;
            bqs[0].b1 = (digital_tf->num_len > 1) ? digital_tf->num[1] / a0 : 0.0;
            bqs[0].b2 = (digital_tf->num_len > 2) ? digital_tf->num[2] / a0 : 0.0;
            bqs[0].a1 = (digital_tf->den_len > 1) ? digital_tf->den[1] / a0 : 0.0;
            bqs[0].a2 = (digital_tf->den_len > 2) ? digital_tf->den[2] / a0 : 0.0;
        } else {
            bqs[0].b0 = 1.0; bqs[0].b1 = 0.0; bqs[0].b2 = 0.0;
            bqs[0].a1 = 0.0; bqs[0].a2 = 0.0;
        }
        bqs[0].gain = digital_tf->gain;
    }

    iir_filter_t *iir = iir_filter_alloc(bqs, num_sections, digital_tf->gain);
    free(bqs);
    tf_discrete_free(digital_tf);
    return iir;
}

iir_filter_t *iir_convert_structure(const iir_filter_t *iir,
                                     digital_structure_t target) {
    if (!iir) return NULL;
    /* For cascade form, structure change is transparent at this level */
    iir_filter_t *new_iir = iir_filter_alloc(iir->sections, iir->num_sections,
                                              iir->overall_gain);
    if (new_iir) new_iir->structure = target;
    return new_iir;
}

/* ==================================================================
 * Least-Squares FIR Design
 * ================================================================== */

/**
 * fir_design_least_squares — L5
 *
 * Minimizes: E = integral_{bands} W(w)|H_desired(w) - H_actual(w)|^2 dw
 *
 * For linear-phase FIR with N odd (Type I), the solution involves
 * solving a system of linear equations:
 *   R * a = p
 * where R is a Toeplitz matrix derived from the frequency band integrals,
 * and a contains the independent filter coefficients.
 *
 * Complexity: O(N * num_bands) for matrix construction
 */
fir_filter_t *fir_design_least_squares(const double *bands,
                                        const double *desired,
                                        const double *weights,
                                        int num_bands, int N) {
    if (!bands || !desired || !weights || num_bands < 1 || N < 3) return NULL;

    double M = (N - 1) / 2.0;
    double *h = (double *)malloc(N * sizeof(double));
    if (!h) return NULL;

    /* Simple LS approximation: weight each frequency point */
    int n;
    for (n = 0; n < N; n++) {
        double sum = 0.0;
        /* Integrate over bands (discrete approximation) */
        int b;
        for (b = 0; b < num_bands; b++) {
            double w_low = bands[2 * b];
            double w_high = bands[2 * b + 1];
            double Hd = desired[b];
            double W = weights[b];

            /* Sample the band at 10 points */
            int s;
            for (s = 0; s < 10; s++) {
                double w = w_low + (w_high - w_low) * s / 9.0;
                double ideal_phase = -w * M;
                sum += W * Hd * cos(w * (n - M) + ideal_phase) / 10.0;
            }
        }
        /* Normalize */
        h[n] = sum / (2.0 * M_PI);
    }

    fir_filter_t *fir = fir_filter_alloc(h, N);
    free(h);
    return fir;
}

/* ==================================================================
 * Coefficient Quantization
 * ================================================================== */

void coef_quantize(double *coef, int num_coefs, int frac_bits) {
    if (!coef || num_coefs < 1 || frac_bits < 1) return;
    double scale = pow(2.0, frac_bits);
    int i;
    for (i = 0; i < num_coefs; i++) {
        coef[i] = round(coef[i] * scale) / scale;
    }
}

/**
 * sos_sensitivity — L5 coefficient sensitivity analysis for SOS filters
 *
 * Measures how pole/zero locations change with small perturbations in
 * biquad coefficients. High sensitivity means small coefficient errors
 * (from quantization) cause large changes in filter response.
 *
 * For a biquad: H(z) = (b0 + b1*z^{-1} + b2*z^{-2}) / (1 + a1*z^{-1} + a2*z^{-2})
 *
 * Pole locations p satisfy: 1 + a1*p^{-1} + a2*p^{-2} = 0
 *                            → p^2 + a1*p + a2 = 0
 *                            → p = (-a1 ± sqrt(a1^2 - 4*a2))/2
 *
 * Sensitivity of pole magnitude to a2:
 *   ∂|p|/∂a2 ≈ 1/(2*|p|)  for complex poles
 *
 * Sensitivity of pole to a1:
 *   ∂p/∂a1 = ± 1/2  (for conjugate pair)
 *
 * Key result: High-Q poles (a2 close to 1, a1 close to ±2) have
 * extreme sensitivity to coefficient errors. This is WHY cascade form
 * (SOS) is essential for high-order IIR filters — direct form would
 * require impractical coefficient precision.
 *
 * The sensitivity metric reported is the maximum dB gain from coefficient
 * error to pole displacement, computed as:
 *   sens_db = 20*log10(|Δp|/|p|) / (|Δa|/|a|)
 *
 * Reference: Jackson, "On the Interaction of Roundoff Noise and
 *             Dynamic Range in Digital Filters" (1970)
 *             Oppenheim & Schafer (2010), Ch. 6.7
 * Complexity: O(num_sections)
 */
double sos_sensitivity(const sos_matrix_t *sos, double *sens) {
    if (!sos || !sos->sos_matrix) return -1.0;

    double max_sens = 0.0;
    int k;

    for (k = 0; k < sos->num_sections; k++) {
        double a1 = sos->sos_matrix[k * 6 + 4]; /* a1 coefficient */
        double a2 = sos->sos_matrix[k * 6 + 5]; /* a2 coefficient */

        /* Compute pole locations from a1, a2 */
        double discriminant = a1 * a1 - 4.0 * a2;
        double pole_mag = sqrt(fabs(a2));

        /* Pole sensitivity analysis
         *
         * For complex poles (discriminant < 0):
         *   |p| = sqrt(a2)
         *   ∂|p|/∂a2 = 1/(2*|p|)
         *
         * For real poles (discriminant >= 0):
         *   p = (-a1 ± sqrt(a1^2 - 4*a2))/2
         *   ∂p/∂a1 = (-1 ± a1/sqrt(a1^2 - 4*a2))/2
         */
        double ds_da2;
        if (fabs(pole_mag) > 1e-10) {
            ds_da2 = 1.0 / (2.0 * pole_mag);
        } else {
            ds_da2 = 0.0;
        }

        /* Sensitivity of pole magnitude to a2 coefficient */
        double sens_a2 = ds_da2 * fabs(a2) / (pole_mag + 1e-15);

        /* Sensitivity of pole magnitude to a1 coefficient
         * (only significant for real poles with large |a1|) */
        double ds_da1 = 0.0;
        if (discriminant >= 0.0) {
            double sqrt_d = sqrt(discriminant);
            /* For the dominant pole */
            double dp_da1_mag = fabs(-0.5 + a1 / (2.0 * sqrt_d + 1e-15));
            ds_da1 = dp_da1_mag / (pole_mag + 1e-15);
        }

        /* Total sensitivity as dB ratio */
        double section_sens = sens_a2 + ds_da1 * fabs(a1);
        double section_sens_db = (section_sens > 1e-15) ?
            20.0 * log10(section_sens) : -100.0;

        if (sens) sens[k] = section_sens_db;
        if (section_sens_db > max_sens) max_sens = section_sens_db;

        /* Adjustment for high-Q: Q = sqrt(a2)/(1-a2) approximation */
        if (fabs(1.0 - a2) > 1e-10 && a2 > 0.0) {
            double Q = sqrt(a2) / (1.0 - a2);
            /* High-Q poles have sensitivity proportional to Q^2 */
            double q_sens = 20.0 * log10(Q * Q + 1.0);
            if (q_sens > max_sens) max_sens = q_sens;
            if (sens) sens[k] = (q_sens > sens[k]) ? q_sens : sens[k];
        }
    }

    return max_sens;
}
