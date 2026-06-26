/**
 * filter_defs.h — Core Type Definitions for Filter Theory
 *
 * L1 Definitions: filter types, coefficient structures,
 * specification parameters, pole-zero representation.
 *
 * Courses: MIT 6.003, Stanford EE102A, Berkeley EE123
 * Reference: Oppenheim & Schafer (2010) Discrete-Time Signal Processing
 */

#ifndef FILTER_DEFS_H
#define FILTER_DEFS_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L1: Filter Type Classifications — 10 canonical filter types
 * ================================================================== */

typedef enum {
    FILTER_LOWPASS   = 0,  /* Pass low frequencies, attenuate high */
    FILTER_HIGHPASS  = 1,  /* Pass high frequencies, attenuate low */
    FILTER_BANDPASS  = 2,  /* Pass a band, attenuate outside */
    FILTER_BANDSTOP  = 3,  /* Attenuate a band (notch), pass elsewhere */
    FILTER_ALLPASS   = 4,  /* Pass all frequencies, modify phase only */
    FILTER_NOTCH     = 5,  /* Narrow band-reject, high Q */
    FILTER_PEAKING   = 6,  /* Narrow band-boost (resonant) */
    FILTER_LOWSHELF  = 7,  /* Boost/cut below a corner frequency */
    FILTER_HIGHSHELF = 8,  /* Boost/cut above a corner frequency */
    FILTER_CUSTOM    = 9   /* User-defined arbitrary response */
} filter_type_t;

/** Approximation type — each represents a distinct optimization criterion:
 *  Butterworth:   maximally flat magnitude at w=0
 *  Chebyshev I:   equiripple in passband, monotonic stopband
 *  Chebyshev II:  monotonic passband, equiripple in stopband
 *  Elliptic:      equiripple in both passband and stopband (Cauer)
 *  Bessel:        maximally flat group delay (best phase linearity)
 *  Gaussian:      Gaussian magnitude shape, minimal time-bandwidth
 *  Legendre:      monotonic with maximum slope at cutoff
 *  Synchronous:   cascaded identical stages, simple tuning */
typedef enum {
    APPROX_BUTTERWORTH   = 0,
    APPROX_CHEBYSHEV_I   = 1,
    APPROX_CHEBYSHEV_II  = 2,
    APPROX_ELLIPTIC      = 3,
    APPROX_BESSEL        = 4,
    APPROX_GAUSSIAN      = 5,
    APPROX_LEGENDRE      = 6,
    APPROX_SYNCHRONOUS   = 7
} approx_type_t;

/** Analog circuit topology for active filter realization.
 *  Reference: Sedra & Smith (2020) Microelectronic Circuits, Ch. 16. */
typedef enum {
    TOPO_SALLEN_KEY      = 0,  /* Sallen-Key, single op-amp, low sensitivity */
    TOPO_MFB             = 1,  /* Multiple Feedback, inverting, stable */
    TOPO_STATE_VARIABLE  = 2,  /* State-variable, simultaneous LP/BP/HP outputs */
    TOPO_BIQUAD          = 3,  /* Two-integrator-loop (Tow-Thomas) biquad */
    TOPO_GIC             = 4,  /* Generalized Immittance Converter, inductor sim */
    TOPO_FDNR            = 5,  /* Frequency-Dependent Negative Resistance */
    TOPO_GYRATOR         = 6,  /* Gyrator-based LC simulation */
    TOPO_LADDER          = 7   /* Doubly-terminated LC ladder (passive prototype) */
} filter_topology_t;

/** Digital filter implementation structure.
 *  DF_I:            non-canonical, separate delay lines
 *  DF_II:           canonical, shared delay line, minimum memory
 *  DF_I_TRANSPOSED: transposed DF_I, zero-input limit cycle free
 *  DF_II_TRANSPOSED: transposed DF_II, preferred for fixed-point
 *  CASCADE_BQ:      cascade of biquads — lowest coefficient sensitivity
 *  PARALLEL_BQ:     parallel biquads — lowest roundoff noise
 *  LATTICE:         orthogonal structure for adaptive filtering
 *  LADDER_WAVE:     wave digital ladder — preserves passivity */
typedef enum {
    DF_I              = 0,
    DF_II             = 1,
    DF_I_TRANSPOSED   = 2,
    DF_II_TRANSPOSED  = 3,
    CASCADE_BQ        = 4,
    PARALLEL_BQ       = 5,
    LATTICE           = 6,
    LADDER_WAVE       = 7
} digital_structure_t;

/** FIR linear-phase type classification.
 *  L2: Symmetry conditions guarantee exact linear phase.
 *  Type I:   N odd,  h[n] = h[N-1-n] — universal (LP, HP, BP, BS)
 *  Type II:  N even, h[n] = h[N-1-n] — no HP/BS (zero at z=-1)
 *  Type III: N odd,  h[n] = -h[N-1-n] — differentiator, Hilbert (zeros at z=1,-1)
 *  Type IV:  N even, h[n] = -h[N-1-n] — differentiator, Hilbert (zero at z=1) */
typedef enum {
    FIR_TYPE_I   = 0,
    FIR_TYPE_II  = 1,
    FIR_TYPE_III = 2,
    FIR_TYPE_IV  = 3
} fir_linear_phase_type_t;

/* ==================================================================
 * L1: Core Data Structures — 7 canonical filter representations
 * ================================================================== */

/** Filter design specification — translation of application needs to parameters.
 *  L1 Definition: All user-facing requirements for a filter design task. */
typedef struct {
    filter_type_t   type;               /* Filter type (LP/HP/BP/BS/...) */
    approx_type_t   approx;             /* Approximation family */
    double          sample_rate;        /* Sampling frequency (Hz), 0 = analog */
    double          fc1;                /* Lower cutoff / corner freq (Hz) */
    double          fc2;                /* Upper cutoff freq (Hz), 0 for LP/HP */
    double          passband_ripple_db; /* Max passband ripple (dB) */
    double          stopband_atten_db;  /* Min stopband attenuation (dB) */
    double          passband_edge;      /* Passband edge frequency (Hz) */
    double          stopband_edge;      /* Stopband edge frequency (Hz) */
    int             order;              /* Filter order, 0 = auto-estimate */
    double          q_factor;           /* Quality factor (for BP/BS), 0 = auto */
    double          gain_db;            /* Passband gain in dB */
} filter_spec_t;

/** Second-order section (SOS / biquad) — canonical building block.
 *  L1: H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
 *  High-order filters are decomposed into cascaded SOS for numerical robustness. */
typedef struct {
    double b0, b1, b2;  /* Numerator (zero) coefficients */
    double a1, a2;      /* Denominator (pole) coefficients, a0 = 1 */
    double gain;        /* Overall section gain */
} biquad_section_t;

/** Analog second-order section — continuous-time biquad.
 *  L1: H(s) = (b2*s^2 + b1*s + b0) / (s^2 + a1*s + a0)
 *  Used as building block for cascaded active filter design. */
typedef struct {
    double b2, b1, b0;  /* Numerator coefficients */
    double a1, a0;      /* Denominator coefficients (a2 = 1 for normalized) */
    double gain;        /* Section gain factor */
} analog_biquad_t;

/** Pole-zero pair — fundamental dynamic characterization.
 *  L1: Complex-conjugate pole or zero in s-plane or z-plane.
 *  Pole location determines natural frequency wn and quality factor Q.
 *  Q = wn / (2*|Re(pole)|), or Q = 1/(2*sin(angle)) for z-plane. */
typedef struct {
    double complex pole;       /* Pole location in complex plane */
    double complex zero;       /* Zero location in complex plane */
    double         q_pole;     /* Q-factor of pole pair: wn/(2*|sigma|) */
    double         q_zero;     /* Q-factor of zero pair */
    double         freq_pole;  /* Natural frequency wn (rad/s or rad) */
    double         freq_zero;  /* Natural frequency of zero (rad/s or rad) */
} pz_pair_t;

/** Complete pole-zero map — full characterization of a linear filter.
 *  L1/L3: Determines stability, frequency response, transient behavior. */
typedef struct {
    pz_pair_t *pz_pairs;       /* Complex-conjugate pole-zero pairs */
    int        num_pairs;      /* Number of second-order sections */
    int        num_real_poles; /* Additional real-axis poles */
    int        num_real_zeros; /* Additional real-axis zeros */
    double    *real_poles;     /* Array of real pole locations */
    double    *real_zeros;     /* Array of real zero locations */
    double     gain;           /* Overall gain factor K */
} pz_map_t;

/** FIR filter descriptor — all-zero filter, inherently stable.
 *  L1: y[n] = sum_{k=0}^{N-1} h[k] * x[n-k]
 *  Finite impulse response means the output depends only on N past inputs. */
typedef struct {
    double       *coeff;        /* Filter coefficients h[0..N-1] */
    int           length;       /* Number of taps N */
    fir_linear_phase_type_t lp_type; /* Linear phase type */
    double        gain_dc;      /* DC gain = sum of coefficients */
    double        group_delay;  /* Constant group delay = (N-1)/2 samples */
} fir_filter_t;

/** IIR filter descriptor — feedback creates poles.
 *  L1: y[n] = sum b_k*x[n-k] - sum a_k*y[n-k]
 *  Infinite impulse response, can achieve sharper transitions with lower order. */
typedef struct {
    biquad_section_t *sections;     /* Array of biquad sections */
    int                num_sections; /* Number of second-order sections */
    digital_structure_t structure;  /* Implementation structure */
    double             overall_gain; /* Overall gain factor */
} iir_filter_t;

/** Filter state for streaming sample-by-sample processing.
 *  L2: Delay-line memory for real-time filter operation.
 *  Supports both FIR (circular delay line) and IIR (biquad state buffers). */
typedef struct {
    double *dline;          /* FIR delay line of length N (circular buffer) */
    int     dline_len;      /* Length of FIR delay line */
    int     dline_idx;      /* Current write position (circular) */
    double *w_buf;          /* IIR state buffer (for DF-II: 2*nsections) */
    int     w_buf_len;      /* Length of IIR state buffer */
    double *v_buf;          /* IIR auxiliary state buffer (for DF-I) */
    double *x_history;      /* Input history for block processing */
    int     hist_len;       /* History buffer length */
} filter_state_t;

/** Frequency response at a single frequency point.
 *  L1: Fundamental filter characterization — H(jw) or H(e^{jw}).
 *  Magnitude shows gain; phase shows delay; group delay shows envelope delay. */
typedef struct {
    double freq_hz;         /* Frequency in Hz */
    double freq_rad;        /* Frequency in rad/s (analog) or rad/sample (digital) */
    double magnitude;       /* |H(jw)| or |H(e^{jw})| (linear) */
    double magnitude_db;    /* 20*log10(|H|) in dB */
    double phase_rad;       /* arg(H) in radians */
    double phase_deg;       /* arg(H) in degrees */
    double group_delay_s;   /* -d(phase)/dw in seconds (or samples for digital) */
    double real_part;       /* Re{H} */
    double imag_part;       /* Im{H} */
} freq_resp_point_t;

/** Complete frequency response vector.
 *  L1/L3: Discrete evaluation of H(w) across a frequency range.
 *  Used for Bode plots, Nyquist plots, Nichols charts. */
typedef struct {
    freq_resp_point_t *points;      /* Array of frequency response points */
    int                num_points;   /* Number of frequency points */
    double             f_start;      /* Start frequency (Hz) */
    double             f_stop;       /* Stop frequency (Hz) */
    double             f_resolution; /* Frequency step (Hz) */
} freq_resp_t;

/** Window function type for FIR filter design.
 *  L2/L5: Different windows trade off main-lobe width vs side-lobe level.
 *  Rectangular:     main lobe narrowest, side lobes -13 dB
 *  Hann:            side lobes -31.5 dB, 18 dB/octave roll-off
 *  Hamming:         side lobes -42.7 dB, optimized first side lobe
 *  Blackman:        side lobes -58 dB, 18 dB/octave roll-off
 *  Blackman-Harris: side lobes -92 dB (4-term), very wide main lobe
 *  Kaiser:          adjustable beta, near-optimal energy concentration
 *  Flat-top:        minimal scalloping loss, for precision spectral analysis
 *  Bartlett:        triangular, side lobes -26.5 dB
 *  Nuttall:         4-term, continuous derivatives at boundary
 *  Tukey:           tapered cosine, alpha controls cosine fraction
 *  Reference: Harris, Proc. IEEE, vol. 66, no. 1, pp. 51-83, 1978. */
typedef enum {
    WINDOW_RECTANGULAR     = 0,
    WINDOW_HANN            = 1,
    WINDOW_HAMMING         = 2,
    WINDOW_BLACKMAN        = 3,
    WINDOW_BLACKMAN_HARRIS = 4,
    WINDOW_KAISER          = 5,
    WINDOW_FLATTOP         = 6,
    WINDOW_BARTLETT        = 7,
    WINDOW_NUTTALL         = 8,
    WINDOW_TUKEY           = 9
} window_type_t;

/** Window function parameters. */
typedef struct {
    window_type_t type;   /* Window family */
    double        alpha;  /* Tukey: taper ratio [0,1]; Kaiser: beta parameter */
    int           length; /* Window length in samples */
} window_params_t;

/* ==================================================================
 * L1: Error Codes — 12 canonical filter error conditions
 * ================================================================== */

typedef enum {
    FILTER_OK              =  0,  /* Success */
    FILTER_ERR_NULL_PTR    = -1,  /* NULL pointer argument */
    FILTER_ERR_INVALID_ORDER = -2,  /* Order out of valid range */
    FILTER_ERR_INVALID_FREQ  = -3,  /* Frequency specification invalid */
    FILTER_ERR_INVALID_RIPPLE = -4, /* Ripple/attenuation spec invalid */
    FILTER_ERR_UNSTABLE    = -5,  /* Resulting filter is unstable */
    FILTER_ERR_NO_CONVERGE = -6,  /* Design algorithm did not converge */
    FILTER_ERR_MEMORY      = -7,  /* Memory allocation failure */
    FILTER_ERR_NOT_IMPL    = -8,  /* Feature not yet implemented */
    FILTER_ERR_NUMERICAL   = -9,  /* Numerical precision issue */
    FILTER_ERR_ZERO_DIV    = -10, /* Division by zero encountered */
    FILTER_ERR_OVERFLOW    = -11  /* Numeric overflow */
} filter_error_t;

/* ==================================================================
 * L1: Convenience inline functions
 * ================================================================== */

/** Initialize a filter_spec_t with safe defaults (1 kHz LP Butterworth). */
static inline filter_spec_t filter_spec_init(void) {
    return (filter_spec_t){
        .type = FILTER_LOWPASS, .approx = APPROX_BUTTERWORTH,
        .sample_rate = 0.0, .fc1 = 1000.0, .fc2 = 0.0,
        .passband_ripple_db = 0.1, .stopband_atten_db = 60.0,
        .passband_edge = 0.0, .stopband_edge = 0.0,
        .order = 0, .q_factor = 0.0, .gain_db = 0.0
    };
}

/** Validate filter specification parameters.
 *  Checks boundary conditions: non-negative frequencies, valid ripple,
 *  consistent band edges for BP/BS filters. Returns non-zero if valid. */
static inline int filter_spec_is_valid(const filter_spec_t *spec) {
    if (!spec) return 0;
    if (spec->fc1 < 0.0) return 0;
    if (spec->passband_ripple_db < 0.0) return 0;
    if (spec->stopband_atten_db < 0.0) return 0;
    if (spec->order < 0) return 0;
    if (spec->type == FILTER_BANDPASS || spec->type == FILTER_BANDSTOP) {
        if (spec->fc2 <= spec->fc1 && spec->fc2 > 0.0) return 0;
    }
    return 1;
}

#ifdef __cplusplus
}
#endif

#endif /* FILTER_DEFS_H */
