/**
 * filter_transfer.h — Transfer Function Representations
 *
 * L1: Polynomial transfer functions, partial fraction expansion
 * L3: Mathematical structures for s-domain and z-domain
 * L4: Stability criterion (pole locations)
 *
 * Courses: MIT 6.003, Berkeley EE123, Cambridge 3F3
 * Reference: Oppenheim & Willsky (1997) Signals and Systems, Ch. 9
 */

#ifndef FILTER_TRANSFER_H
#define FILTER_TRANSFER_H

#include "filter_defs.h"
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L1/L3: Polynomial Transfer Function Representations
 * ================================================================== */

/** Continuous-time transfer function.
 *  L3: H(s) = N(s)/D(s) where s = sigma + j*omega
 *  num[0..num_len-1] defines numerator polynomial b0 + b1*s + ... + bN*s^N
 *  den[0..den_len-1] defines denominator polynomial a0 + a1*s + ... + aM*s^M
 *  Poles of H(s) must be in left half-plane for BIBO stability. */
typedef struct {
    double *num;        /* Numerator coefficients, length = num_len */
    double *den;        /* Denominator coefficients, length = den_len */
    int     num_len;    /* Number of numerator coefficients */
    int     den_len;    /* Number of denominator coefficients */
    double  gain;       /* Overall gain K: H(s) = K * num(s)/den(s) */
} tf_continuous_t;

/** Discrete-time transfer function.
 *  L3: H(z) = N(z)/D(z) where z = e^{j*omega}
 *  num[0..num_len-1] = b0, b1, ..., bN for sum b_k*z^{-k}
 *  den[0..den_len-1] = a0, a1, ..., aM for sum a_k*z^{-k} (a0 typically = 1)
 *  Poles must be inside unit circle |z| < 1 for BIBO stability.
 *  L4: Jury's stability test verifies this condition algebraically. */
typedef struct {
    double *num;        /* Numerator coefficients (b), length = num_len */
    double *den;        /* Denominator coefficients (a), length = den_len */
    int     num_len;    /* Number of numerator coefficients */
    int     den_len;    /* Number of denominator coefficients */
    double  gain;       /* Overall gain */
} tf_discrete_t;

/** Second-order section (SOS) matrix form.
 *  L3: H(z) = prod_{k} (b0k + b1k*z^{-1} + b2k*z^{-2}) / (1 + a1k*z^{-1} + a2k*z^{-2})
 *  This is the preferred representation for high-order IIR filters.
 *  Stored as a flat array of [b0, b1, b2, 1, a1, a2] per section. */
typedef struct {
    double *sos_matrix; /* Flat SOS matrix: 6 doubles per section */
    int     num_sections; /* Number of second-order sections */
    double  gain;       /* Overall gain */
} sos_matrix_t;

/* ==================================================================
 * L3: Partial Fraction Expansion
 * ================================================================== */

/** Single term in partial fraction expansion.
 *  L3: For real pole p: term = r/(s-p) or r/(1-p*z^{-1})
 *  For complex pair: two conjugate terms combined into one biquad.
 *  Essential for parallel filter realization and inverse Laplace/Z transform. */
typedef struct {
    double complex residue;  /* Residue at the pole */
    double complex pole_loc; /* Pole location */
    int            is_complex; /* 1 if part of complex-conjugate pair */
} pfe_term_t;

/** Partial fraction expansion of a transfer function.
 *  L3: H(s) = sum r_k/(s-p_k) + direct_term
 *  Used for parallel-form IIR realization and analog circuit synthesis. */
typedef struct {
    pfe_term_t *terms;       /* Array of PFE terms */
    int         num_terms;   /* Number of terms */
    double      direct_term; /* Direct feedthrough term (if proper) */
} pfe_expansion_t;

/* ==================================================================
 * L3: Transfer Function API
 * ================================================================== */

/** Allocate and initialize a continuous-time transfer function.
 *  @param num_len  Number of numerator coefficients
 *  @param den_len  Number of denominator coefficients
 *  @return Allocated tf_continuous_t, or NULL on failure
 *  Complexity: O(num_len + den_len) */
tf_continuous_t *tf_continuous_alloc(int num_len, int den_len);

/** Free a continuous-time transfer function. */
void tf_continuous_free(tf_continuous_t *tf);

/** Allocate and initialize a discrete-time transfer function. */
tf_discrete_t *tf_discrete_alloc(int num_len, int den_len);

/** Free a discrete-time transfer function. */
void tf_discrete_free(tf_discrete_t *tf);

/** Evaluate H(s) at a complex frequency s = sigma + j*omega.
 *  L3: Direct polynomial evaluation using Horner's method.
 *  @param tf   Transfer function
 *  @param s    Complex frequency point
 *  @return     H(s) as complex double
 *  Complexity: O(max(num_len, den_len)) */
double complex tf_continuous_eval(const tf_continuous_t *tf, double complex s);

/** Evaluate H(z) at a complex frequency z = e^{j*omega}.
 *  L3: Direct polynomial evaluation using Horner's method.
 *  @param tf   Transfer function
 *  @param z    Complex point on unit circle or elsewhere
 *  @return     H(z) as complex double
 *  Complexity: O(max(num_len, den_len)) */
double complex tf_discrete_eval(const tf_discrete_t *tf, double complex z);

/** Compute magnitude of H(s) at real frequency omega.
 *  Wrapper around tf_continuous_eval for convenience.
 *  @return |H(j*omega)| */
double tf_continuous_mag_at(const tf_continuous_t *tf, double omega);

/** Compute magnitude of H(z) at normalized frequency omega.
 *  @param omega  Normalized frequency in radians/sample [0, pi]
 *  @return |H(e^{j*omega})| */
double tf_discrete_mag_at(const tf_discrete_t *tf, double omega);

/* ==================================================================
 * L3: Transfer Function Manipulation
 * ================================================================== */

/** Multiply two continuous-time transfer functions: H(s) = H1(s) * H2(s).
 *  L3: Polynomial convolution of numerator and denominator.
 *  Used for cascading filter stages.
 *  Complexity: O(N1*N2 + M1*M2) where N=num_len, M=den_len */
tf_continuous_t *tf_continuous_multiply(const tf_continuous_t *tf1,
                                         const tf_continuous_t *tf2);

/** Multiply two discrete-time transfer functions: H(z) = H1(z) * H2(z). */
tf_discrete_t *tf_discrete_multiply(const tf_discrete_t *tf1,
                                     const tf_discrete_t *tf2);

/** Convert continuous-time transfer function to cascade of analog biquads.
 *  L3: Factoring high-order polynomials into second-order sections.
 *  @param tf     Input transfer function
 *  @param bqs    Output array of analog biquads (caller allocates)
 *  @param max_bq Maximum number of biquads
 *  @return       Number of biquads filled, or -1 on error */
int tf_continuous_to_biquads(const tf_continuous_t *tf,
                              analog_biquad_t *bqs, int max_bq);

/** Convert discrete-time transfer function to cascade of digital biquads. */
int tf_discrete_to_biquads(const tf_discrete_t *tf,
                            biquad_section_t *bqs, int max_bq);

/* ==================================================================
 * L3: Pole-Zero Map Conversion
 * ================================================================== */

/** Build a pole-zero map from a continuous-time transfer function.
 *  Computes roots of numerator and denominator polynomials.
 *  Complexity: O(N^3) for polynomial root-finding (eigenvalue method) */
pz_map_t *tf_continuous_to_pzmap(const tf_continuous_t *tf);

/** Build a pole-zero map from a discrete-time transfer function. */
pz_map_t *tf_discrete_to_pzmap(const tf_discrete_t *tf);

/** Free a pole-zero map. */
void pz_map_free(pz_map_t *pz);

/** Build a continuous-time transfer function from a pole-zero map.
 *  L3: H(s) = K * prod(s - z_k) / prod(s - p_k)
 *  Polynomial expansion from roots. */
tf_continuous_t *pzmap_to_tf_continuous(const pz_map_t *pz);

/** Build a discrete-time transfer function from a pole-zero map. */
tf_discrete_t *pzmap_to_tf_discrete(const pz_map_t *pz);

/* ==================================================================
 * L2: Frequency Transformations
 * ================================================================== */

/** Lowpass-to-Lowpass frequency transformation (analog).
 *  L2: s -> s/omega_c — scales the cutoff frequency.
 *  Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 7 */
void tf_analog_lp_to_lp(tf_continuous_t *tf, double new_wc);

/** Lowpass-to-Highpass frequency transformation (analog).
 *  L2: s -> omega_c/s — converts LP prototype to HP.
 *  Poles at origin become zeros, zeros become poles. */
void tf_analog_lp_to_hp(tf_continuous_t *tf, double new_wc);

/** Lowpass-to-Bandpass frequency transformation (analog).
 *  L2: s -> (s^2 + w0^2)/(s*B) — doubles the order.
 *  Each pole becomes a pair; w0 = sqrt(wL*wH), B = wH - wL. */
tf_continuous_t *tf_analog_lp_to_bp(const tf_continuous_t *tf,
                                     double w0, double bw);

/** Lowpass-to-Bandstop frequency transformation (analog).
 *  L2: s -> (s*B)/(s^2 + w0^2) — doubles the order. */
tf_continuous_t *tf_analog_lp_to_bs(const tf_continuous_t *tf,
                                     double w0, double bw);

/* ==================================================================
 * L4: Stability Tests
 * ================================================================== */

/** Routh-Hurwitz stability criterion for continuous-time systems.
 *  L4: All poles must have negative real parts for BIBO stability.
 *  Constructs the Routh array and counts sign changes in first column.
 *  A system is stable iff all elements in the first column have the same sign.
 *  Reference: Dorf & Bishop, Modern Control Systems (2017), Ch. 6
 *  @param tf  Transfer function to test
 *  @return    1 if stable, 0 if unstable, -1 on error */
int tf_continuous_is_stable(const tf_continuous_t *tf);

/** Jury's stability test for discrete-time systems.
 *  L4: All poles must have magnitude < 1 for BIBO stability.
 *  Constructs the Jury table and checks the necessary and sufficient conditions.
 *  Reference: Jury, "Theory and Application of the z-Transform Method" (1964)
 *  @param tf  Transfer function to test
 *  @return    1 if stable, 0 if unstable, -1 on error */
int tf_discrete_is_stable(const tf_discrete_t *tf);

/** Check if an SOS matrix represents a stable filter.
 *  Verifies that each section's poles are inside the unit circle.
 *  For biquad (1 + a1*z^{-1} + a2*z^{-2}):
 *    Stability requires |a2| < 1 and |a1| < 1 + a2.
 *  @return 1 if stable, 0 if any section is unstable */
int sos_is_stable(const sos_matrix_t *sos);

/* ==================================================================
 * L4: Bilinear Transform (s <-> z mapping)
 * ================================================================== */

/** Bilinear transform: s = 2*fs * (z-1)/(z+1)
 *  L4: Maps analog s-plane to digital z-plane.
 *  Maps jw-axis to unit circle (one-to-one).
 *  Introduces frequency warping: w_digital = 2*fs*tan(w_analog/(2*fs)).
 *  Reference: Oppenheim & Schafer (2010), Sec. 7.1
 *  @param s     Analog frequency variable
 *  @param fs    Sampling frequency (Hz)
 *  @return      Corresponding z-plane value */
double complex bilinear_s_to_z(double complex s, double fs);

/** Inverse bilinear transform: z -> s = 2*fs * (z-1)/(z+1) */
double complex bilinear_z_to_s(double complex z, double fs);

/** Pre-warp an analog cutoff frequency for bilinear transform.
 *  L4: w_a = 2*fs*tan(w_d/(2*fs))
 *  Essential for accurate digital filter design.
 *  @param w_d  Desired digital cutoff (rad/sample)
 *  @param fs   Sampling frequency (Hz)
 *  @return     Corresponding analog frequency (rad/s) */
double bilinear_prewarp(double w_d, double fs);

/** Convert an analog transfer function to digital via bilinear transform.
 *  L4/L5: Maps H(s) -> H(z) by substituting s = 2*fs*(z-1)/(z+1).
 *  Preserves stability: LHP poles -> inside unit circle.
 *  Preserves the number of poles and zeros.
 *  Does NOT preserve impulse response or phase linearity.
 *  Complexity: O((max(N,M))^3) due to polynomial operations */
tf_discrete_t *bilinear_transform(const tf_continuous_t *tf_analog, double fs);

/* ==================================================================
 * L5: Impulse Invariance Method
 * ================================================================== */

/** Impulse Invariance: sample the continuous-time impulse response.
 *  L5: h[n] = T * h_c(n*T) where h_c(t) is the analog impulse response.
 *  H(z) = sum Residues/(1 - e^{p_k*T} * z^{-1})
 *  Preserves impulse response shape but may cause aliasing.
 *  Only suitable for bandlimited (lowpass/bp) filters.
 *  @param tf_analog  Analog transfer function
 *  @param fs         Sampling frequency (Hz)
 *  @return           Digital transfer function, or NULL on error */
tf_discrete_t *impulse_invariance(const tf_continuous_t *tf_analog, double fs);

/* ==================================================================
 * L3: Math Utility Functions (implemented in filter_math.c)
 * ================================================================== */

double binomial(int n, int k);
double bessel_i0(double x);
double bessel_i0_ratio(double x, double x0);
int polynomial_multiply(const double *A, int lenA,
                        const double *B, int lenB,
                        double *C, int max_lenC);
int polynomial_add(const double *A, int lenA,
                   const double *B, int lenB,
                   double *C, int max_lenC);
int polynomial_derivative(const double *coeff, int degree,
                           double *deriv);
double polynomial_normalize(double *coeff, int degree);
double complex polynomial_eval_complex(const double *coeff, int degree,
                                        double complex x);
int binomial_expand(int n, double *coeff);
int binomial_expand_signed(int n, double *coeff);
int durand_kerner(const double *coeff, int degree,
                   double complex *roots, int max_iter);
double trapezoidal_integrate(double (*f)(double, void*), void *params,
                              double a, double b, int n);
double simpson_integrate(double (*f)(double, void*), void *params,
                          double a, double b, int n);
double bisection_find_root(double (*f)(double, void*), void *params,
                            double a, double b, double tol);
double newton_raphson(double (*f)(double, void*),
                       double (*fprime)(double, void*),
                       void *params, double x0, double tol, int max_iter);
int gaussian_elimination(double *A, double *b, int n);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_TRANSFER_H */
