/**
 * filter_digital.c — Digital Filter Implementation (FIR and IIR)
 *
 * L2: Direct Form I/II, streaming sample processing, block processing
 * L5: FIR window design, CIC filters, half-band filters
 *
 * Reference:
 *   Oppenheim & Schafer (2010) Discrete-Time Signal Processing
 *   Hogenauer, "An Economical Class of Digital Filters..." (1981)
 *   Kaiser, "Nonrecursive Digital Filter Design Using I_0-sinh" (1974)
 */

#include "filter_digital.h"
#include "filter_design.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * FIR Filter Allocation
 * ================================================================== */

fir_filter_t *fir_filter_alloc(const double *coeff, int length) {
    if (!coeff || length < 1) return NULL;

    fir_filter_t *fir = (fir_filter_t *)malloc(sizeof(fir_filter_t));
    if (!fir) return NULL;

    fir->coeff = (double *)malloc(length * sizeof(double));
    if (!fir->coeff) { free(fir); return NULL; }

    int i;
    double sum = 0.0;
    for (i = 0; i < length; i++) {
        fir->coeff[i] = coeff[i];
        sum += coeff[i];
    }

    fir->length = length;
    fir->gain_dc = sum;
    fir->group_delay = (length - 1) / 2.0;
    fir->lp_type = fir_determine_linear_phase_type(coeff, length, 1e-10);

    return fir;
}

void fir_filter_free(fir_filter_t *fir) {
    if (fir) {
        free(fir->coeff);
        free(fir);
    }
}

/* ==================================================================
 * FIR State Management
 * ================================================================== */

filter_state_t *fir_state_alloc(int length) {
    if (length < 1) return NULL;

    filter_state_t *state = (filter_state_t *)malloc(sizeof(filter_state_t));
    if (!state) return NULL;

    state->dline = (double *)calloc(length, sizeof(double));
    state->dline_len = length;
    state->dline_idx = 0;
    state->w_buf = NULL;
    state->w_buf_len = 0;
    state->v_buf = NULL;
    state->x_history = NULL;
    state->hist_len = 0;

    if (!state->dline) { free(state); return NULL; }
    return state;
}

void filter_state_free(filter_state_t *state) {
    if (state) {
        free(state->dline);
        free(state->w_buf);
        free(state->v_buf);
        free(state->x_history);
        free(state);
    }
}

void filter_state_reset(filter_state_t *state) {
    if (!state) return;
    int i;
    if (state->dline) {
        for (i = 0; i < state->dline_len; i++) state->dline[i] = 0.0;
        state->dline_idx = 0;
    }
    if (state->w_buf) {
        for (i = 0; i < state->w_buf_len; i++) state->w_buf[i] = 0.0;
    }
    if (state->v_buf) {
        for (i = 0; i < state->w_buf_len; i++) state->v_buf[i] = 0.0;
    }
}

/* ==================================================================
 * FIR Streaming Processing
 * ================================================================== */

/**
 * fir_process_sample — L2 streaming FIR
 *
 * Direct form FIR with circular buffer:
 *   y[n] = sum_{k=0}^{N-1} h[k] * x[n-k]
 *
 * The circular buffer stores the most recent N input samples.
 * dline_idx points to the oldest sample (next to be overwritten).
 *
 * Complexity: O(N) per sample
 */
double fir_process_sample(const fir_filter_t *fir, filter_state_t *state,
                           double x) {
    if (!fir || !state || !state->dline) return 0.0;

    /* Write new sample to circular buffer (overwrite oldest) */
    state->dline[state->dline_idx] = x;

    /* Compute dot product h * dline, wrapping around */
    double y = 0.0;
    int i;
    int idx = state->dline_idx;
    for (i = 0; i < fir->length; i++) {
        y += fir->coeff[i] * state->dline[idx];
        idx--;
        if (idx < 0) idx = state->dline_len - 1;
    }

    /* Advance write pointer */
    state->dline_idx++;
    if (state->dline_idx >= state->dline_len)
        state->dline_idx = 0;

    return y;
}

/**
 * fir_process_block — L2 block FIR processing
 *
 * Processes a contiguous block of samples.
 * Complexity: O(N * block_len)
 */
void fir_process_block(const fir_filter_t *fir, filter_state_t *state,
                        const double *x, double *y, int block_len) {
    if (!fir || !state || !x || !y) return;
    int i;
    for (i = 0; i < block_len; i++) {
        y[i] = fir_process_sample(fir, state, x[i]);
    }
}

/* ==================================================================
 * IIR Filter Allocation
 * ================================================================== */

iir_filter_t *iir_filter_alloc(const biquad_section_t *sections,
                                int num_sections, double overall_gain) {
    if (!sections || num_sections < 1) return NULL;

    iir_filter_t *iir = (iir_filter_t *)malloc(sizeof(iir_filter_t));
    if (!iir) return NULL;

    iir->sections = (biquad_section_t *)malloc(num_sections * sizeof(biquad_section_t));
    if (!iir->sections) { free(iir); return NULL; }

    int i;
    for (i = 0; i < num_sections; i++) iir->sections[i] = sections[i];

    iir->num_sections = num_sections;
    iir->structure = CASCADE_BQ;
    iir->overall_gain = overall_gain;

    return iir;
}

void iir_filter_free(iir_filter_t *iir) {
    if (iir) {
        free(iir->sections);
        free(iir);
    }
}

/* ==================================================================
 * IIR Stream Processing — Direct Form II (Transposed)
 * ================================================================== */

/**
 * iir_process_sample_df2 — L2 DF-II biquad cascade
 *
 * Direct Form II structure for each biquad section:
 *   w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
 *   y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
 *
 * DF-II uses only 2 state variables per section (w[n-1], w[n-2]),
 * making it the most memory-efficient structure.
 *
 * The state buffer layout for ns sections:
 *   w_buf[2*k]   = w_k[n-1]
 *   w_buf[2*k+1] = w_k[n-2]
 *
 * Complexity: O(num_sections) per sample
 */
double iir_process_sample_df2(const iir_filter_t *iir, filter_state_t *state,
                                double x) {
    if (!iir || !state) return 0.0;

    /* Allocate state buffer on first call */
    if (!state->w_buf) {
        state->w_buf_len = iir->num_sections * 2;
        state->w_buf = (double *)calloc(state->w_buf_len, sizeof(double));
        if (!state->w_buf) return 0.0;
    }

    double y = x * iir->overall_gain;
    int k;
    for (k = 0; k < iir->num_sections; k++) {
        biquad_section_t *bq = &iir->sections[k];
        double w_n1 = state->w_buf[2 * k];
        double w_n2 = state->w_buf[2 * k + 1];

        /* Compute intermediate node w[n] */
        double w_n = y - bq->a1 * w_n1 - bq->a2 * w_n2;

        /* Compute output of this section */
        y = bq->b0 * w_n + bq->b1 * w_n1 + bq->b2 * w_n2;

        /* Update state */
        state->w_buf[2 * k + 1] = w_n1;  /* w[n-2] = w[n-1] */
        state->w_buf[2 * k]     = w_n;    /* w[n-1] = w[n] */
    }

    return y;
}

/**
 * iir_process_sample_df1 — L2 DF-I biquad cascade
 *
 * Direct Form I:
 *   y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
 *        - a1*y[n-1] - a2*y[n-2]
 *
 * DF-I requires 4 state variables per section but is less susceptible
 * to limit cycle oscillations in fixed-point arithmetic.
 *
 * Complexity: O(num_sections) per sample
 */
double iir_process_sample_df1(const iir_filter_t *iir, filter_state_t *state,
                                double x) {
    if (!iir || !state) return 0.0;

    if (!state->v_buf) {
        state->w_buf_len = iir->num_sections * 4;
        state->v_buf = (double *)calloc(state->w_buf_len, sizeof(double));
        if (!state->v_buf) return 0.0;
    }

    double y = x * iir->overall_gain;
    int k;
    for (k = 0; k < iir->num_sections; k++) {
        biquad_section_t *bq = &iir->sections[k];
        double *st = &state->v_buf[4 * k];
        double xn1 = st[0], xn2 = st[1], yn1 = st[2], yn2 = st[3];

        double yn = bq->b0 * y + bq->b1 * xn1 + bq->b2 * xn2
                    - bq->a1 * yn1 - bq->a2 * yn2;

        /* Shift state */
        st[1] = xn1; st[0] = y;
        st[3] = yn1; st[2] = yn;
        y = yn;
    }
    return y;
}

void iir_process_block(const iir_filter_t *iir, filter_state_t *state,
                        const double *x, double *y, int block_len) {
    if (!iir || !state || !x || !y) return;
    int i;
    for (i = 0; i < block_len; i++) {
        y[i] = iir_process_sample_df2(iir, state, x[i]);
    }
}

/* ==================================================================
 * FIR Window Design
 * ================================================================== */

/**
 * fir_design_window — L5: Window-based FIR design
 *
 * h[n] = h_ideal[n] * w[n], for n = 0..N-1
 *
 * Ideal lowpass impulse response (sinc function):
 *   h_ideal[n] = (wc/pi) * sinc(wc*(n - (N-1)/2) / pi)
 *   where sinc(x) = sin(pi*x)/(pi*x)
 *
 * At the center (n = (N-1)/2), h_ideal = wc/pi (the l'Hopital limit).
 *
 * The window tapers the ideal response to reduce Gibbs phenomenon
 * (overshoot and ringing at discontinuities).
 */
fir_filter_t *fir_design_window(double cutoff, int num_taps,
                                 window_type_t win_type, double win_alpha) {
    if (cutoff <= 0.0 || cutoff >= 0.5 || num_taps < 3) return NULL;

    double wc = 2.0 * M_PI * cutoff; /* Normalized cutoff in rad/sample */
    double M = (num_taps - 1) / 2.0;

    /* Allocate temporary arrays */
    double *h_ideal = (double *)malloc(num_taps * sizeof(double));
    double *window  = (double *)malloc(num_taps * sizeof(double));
    if (!h_ideal || !window) {
        free(h_ideal); free(window); return NULL;
    }

    /* Compute h_ideal[n] */
    int n;
    for (n = 0; n < num_taps; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h_ideal[n] = wc / M_PI; /* Limit at center */
        } else {
            h_ideal[n] = sin(wc * t) / (M_PI * t);
        }
    }

    /* Generate window */
    window_params_t wp = {win_type, win_alpha, num_taps};
    window_generate(window, &wp);

    /* Apply window: h[n] = h_ideal[n] * w[n] */
    double *coeff = (double *)malloc(num_taps * sizeof(double));
    if (!coeff) { free(h_ideal); free(window); return NULL; }

    for (n = 0; n < num_taps; n++) {
        coeff[n] = h_ideal[n] * window[n];
    }

    fir_filter_t *fir = fir_filter_alloc(coeff, num_taps);
    free(h_ideal); free(window); free(coeff);
    return fir;
}

/* ==================================================================
 * Specialized FIR Designs
 * ================================================================== */

/**
 * fir_design_differentiator — L5
 *
 * Ideal differentiator frequency response: H(w) = j*w
 * Impulse response: h[n] = cos(pi*n)/n - sin(pi*n)/(pi*n^2)
 * for n != 0, h[0] = 0.
 *
 * Type III (odd length) or IV (even length), odd symmetry.
 * Amplitude response proportional to frequency (up to Nyquist).
 */
fir_filter_t *fir_design_differentiator(int num_taps, window_type_t win_type,
                                         double win_alpha) {
    if (num_taps < 3) return NULL;

    double M = (num_taps - 1) / 2.0;
    double *h = (double *)malloc(num_taps * sizeof(double));
    double *w = (double *)malloc(num_taps * sizeof(double));
    if (!h || !w) { free(h); free(w); return NULL; }

    int n;
    for (n = 0; n < num_taps; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h[n] = 0.0;
        } else {
            h[n] = cos(M_PI * t) / t - sin(M_PI * t) / (M_PI * t * t);
        }
    }

    window_params_t wp = {win_type, win_alpha, num_taps};
    window_generate(w, &wp);

    double *coeff = (double *)malloc(num_taps * sizeof(double));
    if (!coeff) { free(h); free(w); return NULL; }
    for (n = 0; n < num_taps; n++) coeff[n] = h[n] * w[n];

    fir_filter_t *fir = fir_filter_alloc(coeff, num_taps);
    free(h); free(w); free(coeff);
    return fir;
}

/**
 * fir_design_hilbert — L5
 *
 * Ideal Hilbert transformer: H(w) = -j*sign(w) for 0 < |w| < pi
 * Impulse response: h[n] = 0 for n even (relative to center),
 *   h[n] = 2/(pi*n) for n odd.
 *
 * Produces a 90-degree phase shift across the band.
 * Type III (odd length, odd symmetry).
 */
fir_filter_t *fir_design_hilbert(int num_taps, window_type_t win_type,
                                  double win_alpha) {
    if (num_taps < 5 || num_taps % 2 == 0) return NULL; /* Must be odd for Type III */

    double M = (num_taps - 1) / 2.0;
    double *h = (double *)malloc(num_taps * sizeof(double));
    double *w = (double *)malloc(num_taps * sizeof(double));
    if (!h || !w) { free(h); free(w); return NULL; }

    int n;
    for (n = 0; n < num_taps; n++) {
        double t = n - M;
        if (fabs(t) < 1e-12) {
            h[n] = 0.0;
        } else {
            int k = (int)t;
            if (k < 0) k = -k;
            if (k % 2 == 0) {
                h[n] = 0.0;
            } else {
                h[n] = 2.0 / (M_PI * t);
            }
        }
    }

    window_params_t wp = {win_type, win_alpha, num_taps};
    window_generate(w, &wp);

    double *coeff = (double *)malloc(num_taps * sizeof(double));
    if (!coeff) { free(h); free(w); return NULL; }
    for (n = 0; n < num_taps; n++) coeff[n] = h[n] * w[n];

    fir_filter_t *fir = fir_filter_alloc(coeff, num_taps);
    free(h); free(w); free(coeff);
    return fir;
}

/* ==================================================================
 * CIC Filter
 * ================================================================== */

/**
 * cic_magnitude — L5
 *
 * Hogenauer CIC filter magnitude response:
 *   |H(f)| = |sin(pi*R*M*f) / sin(pi*f)|^N
 * where R = rate change, M = differential delay, N = stages.
 *
 * The response is a sinc^N function — a lowpass filter with
 * periodic notches at multiples of fs/R (used for anti-aliasing).
 */
double cic_magnitude(double freq_norm, int rate_change, int num_stages,
                      int diff_delay) {
    if (freq_norm <= 0.0 || rate_change < 1 || num_stages < 1 || diff_delay < 1)
        return -1.0;

    double w = 2.0 * M_PI * freq_norm;
    double num = sin((double)rate_change * diff_delay * w / 2.0);
    double den = sin(w / 2.0);

    if (fabs(den) < 1e-15) {
        /* L'Hopital: limit as w->0 gives R*M */
        return pow((double)(rate_change * diff_delay), num_stages);
    }

    return pow(fabs(num / den), num_stages);
}

/**
 * cic_compensation_filter — L5
 *
 * The CIC filter has sinc^N passband droop. The compensation FIR
 * approximates the inverse response: ~ 1/sinc^N(f) in the passband.
 *
 * Typically a short FIR (3-7 taps) sufficient for most applications.
 * Reference: Xilinx CIC Compiler documentation, Altera AN455.
 */
fir_filter_t *cic_compensation_filter(int rate_change, int num_stages,
                                       int diff_delay, int comp_taps) {
    if (rate_change < 1 || num_stages < 1 || diff_delay < 1 || comp_taps < 1)
        return NULL;

    double M = (comp_taps - 1) / 2.0;
    double *coeff = (double *)malloc(comp_taps * sizeof(double));
    if (!coeff) return NULL;

    /* Simple inverse-sinc approximation:
     * Design at a few frequency points, inverse-DFT */
    int n;
    for (n = 0; n < comp_taps; n++) {
        double w = M_PI * (n - M) / (double)comp_taps;
        double freq = fabs(w) / (2.0 * M_PI);
        if (freq > 0.48) freq = 0.48;
        double mag = cic_magnitude(freq, rate_change, num_stages, diff_delay);
        if (mag > 1e-10) {
            coeff[n] = 1.0 / mag;
        } else {
            coeff[n] = 1.0;
        }
    }

    fir_filter_t *fir = fir_filter_alloc(coeff, comp_taps);
    free(coeff);
    return fir;
}

/* ==================================================================
 * Half-Band Filter
 * ================================================================== */

fir_filter_t *fir_design_halfband(int num_taps, window_type_t win_type) {
    if (num_taps < 3 || num_taps % 2 == 0) return NULL;
    return fir_design_window(0.25, num_taps, win_type, 5.0);
}

int fir_is_halfband(const fir_filter_t *fir, double tolerance) {
    if (!fir) return 0;
    if (fir->length % 2 == 0) return 0;

    int center = fir->length / 2;
    int n;
    for (n = 0; n < fir->length; n++) {
        if (n != center && (n - center) % 2 == 0) {
            if (fabs(fir->coeff[n]) > tolerance) return 0;
        }
    }
    return 1;
}

/* ==================================================================
 * Linear Phase Type Detection
 * ================================================================== */

/**
 * fir_determine_linear_phase_type — L2
 *
 * Checks coefficient symmetry/antisymmetry:
 * - Even symmetry + odd length  = Type I
 * - Even symmetry + even length = Type II
 * - Odd symmetry + odd length   = Type III
 * - Odd symmetry + even length  = Type IV
 *
 * The symmetry implies h[n] = +- h[N-1-n].
 *
 * Complexity: O(N)
 */
int fir_determine_linear_phase_type(const double *coeff, int length,
                                     double tol) {
    if (!coeff || length < 1) return -1;

    int is_even = 1, is_odd = 1;
    int n;
    int half = length / 2;
    for (n = 0; n < half; n++) {
        double diff_plus  = coeff[n] - coeff[length - 1 - n];
        double diff_minus = coeff[n] + coeff[length - 1 - n];
        if (fabs(diff_plus) > tol)  is_even = 0;
        if (fabs(diff_minus) > tol) is_odd  = 0;
    }

    if (is_even && length % 2 == 1) return FIR_TYPE_I;
    if (is_even && length % 2 == 0) return FIR_TYPE_II;
    if (is_odd  && length % 2 == 1) return FIR_TYPE_III;
    if (is_odd  && length % 2 == 0) return FIR_TYPE_IV;
    return -1;
}
