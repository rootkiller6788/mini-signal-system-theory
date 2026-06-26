/**
 * filter_analog.h — Analog Filter Prototypes and Design
 *
 * L1: Butterworth, Chebyshev I/II, Elliptic, Bessel prototype polynomials
 * L5: Analog filter design algorithms — pole placement, coefficient computation
 * L6: Canonical problem: design a normalized LP prototype, transform to desired spec
 *
 * Courses: MIT 6.003, Stanford EE102A, Berkeley EE105, ETH 227-0427
 * Reference: Van Valkenburg, Analog Filter Design (1982)
 *             Zverev, Handbook of Filter Synthesis (1967)
 */

#ifndef FILTER_ANALOG_H
#define FILTER_ANALOG_H

#include "filter_defs.h"
#include "filter_transfer.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L5: Butterworth Filter — Maximally Flat Magnitude
 * ================================================================== */

/** Compute Butterworth prototype poles for order N.
 *  L5: Poles lie on a circle of radius wc in the left half-plane.
 *  p_k = wc * exp(j*pi*(0.5 + (2k+N-1)/(2N))) for k = 0, 1, ..., N-1
 *  Magnitude: |H(jw)| = 1/sqrt(1 + (w/wc)^{2N})
 *  First 2N-1 derivatives of |H(jw)| are zero at w=0 — maximally flat.
 *  Reference: Butterworth, "On the Theory of Filter Amplifiers" (1930)
 *  @param order  Filter order (>= 1)
 *  @param wc     Cutoff frequency (rad/s)
 *  @param poles  Output array of pole locations (length >= order)
 *  @return       FILTER_OK or error code
 *  Complexity: O(N) */
int butterworth_poles(int order, double wc, double complex *poles);

/** Design a normalized Butterworth lowpass prototype.
 *  L5: Returns transfer function H(s) with wc = 1 rad/s, passband gain = 1.
 *  @param order  Filter order
 *  @return       Normalized prototype, or NULL on error */
tf_continuous_t *butterworth_prototype(int order);

/** Denormalize a Butterworth prototype to desired cutoff and gain.
 *  L2: frequency scaling: s -> s/wc, then adjust gain.
 *  @param proto  Normalized prototype (wc=1)
 *  @param wc     Desired cutoff (rad/s)
 *  @param gain   Passband gain (linear, not dB) */
void butterworth_denormalize(tf_continuous_t *proto, double wc, double gain);

/* ==================================================================
 * L5: Chebyshev Type I — Equiripple in Passband
 * ================================================================== */

/** Compute Chebyshev Type I poles for order N.
 *  L5: Poles lie on an ellipse in the left half-plane.
 *  The ellipse intersects the imaginary axis at +-j*wc.
 *  Magnitude: |H(jw)| = 1/sqrt(1 + epsilon^2 * T_N^2(w/wc))
 *  where T_N is the Chebyshev polynomial of the first kind.
 *  epsilon = sqrt(10^{ripple_dB/10} - 1)
 *  Reference: Chebyshev (1899), later refined by Cauer and Darlington
 *  @param order      Filter order (>= 1)
 *  @param wc         Cutoff frequency (rad/s)
 *  @param ripple_db  Passband ripple in dB
 *  @param poles      Output array (length >= order)
 *  @return           FILTER_OK or error code
 *  Complexity: O(N) */
int chebyshev1_poles(int order, double wc, double ripple_db,
                      double complex *poles);

/** Design a normalized Chebyshev I lowpass prototype.
 *  @param order      Filter order
 *  @param ripple_db  Passband ripple in dB
 *  @return           Normalized prototype (wc=1), or NULL on error */
tf_continuous_t *chebyshev1_prototype(int order, double ripple_db);

/** Estimate required Chebyshev I order given design spec.
 *  L5: N >= acosh(sqrt((10^{As/10}-1)/(10^{Ap/10}-1))) / acosh(ws/wp)
 *  where Ap = passband ripple (dB), As = stopband attenuation (dB),
 *  ws/wp = stopband/passband edge ratio.
 *  Reference: Daniels, Approximation Methods for Electronic Filter Design (1974)
 *  @param passband_ripple_db  Ap
 *  @param stopband_atten_db   As
 *  @param stopband_edge       ws (rad/s)
 *  @param passband_edge       wp (rad/s)
 *  @return                    Minimum required order */
int chebyshev1_estimate_order(double passband_ripple_db,
                               double stopband_atten_db,
                               double stopband_edge,
                               double passband_edge);

/* ==================================================================
 * L5: Chebyshev Type II (Inverse Chebyshev) — Equiripple in Stopband
 * ================================================================== */

/** Compute Chebyshev Type II poles and zeros.
 *  L5: Monotonic in passband, equiripple in stopband.
 *  Has finite transmission zeros on the jw-axis.
 *  Poles are reciprocals of Chebyshev I poles.
 *  Magnitude: |H(jw)| = 1/sqrt(1 + epsilon^2 * T_N^2(wc/w)^{-1})
 *  @param order      Filter order
 *  @param wc         Cutoff frequency (rad/s)
 *  @param atten_db   Minimum stopband attenuation (dB)
 *  @param poles      Output array of pole locations (length >= order)
 *  @param zeros      Output array of zero locations (length >= order)
 *  @return           FILTER_OK or error code */
int chebyshev2_pz(int order, double wc, double atten_db,
                   double complex *poles, double complex *zeros);

/** Design a normalized Chebyshev II lowpass prototype. */
tf_continuous_t *chebyshev2_prototype(int order, double atten_db);

/* ==================================================================
 * L5: Elliptic (Cauer) Filter — Equiripple in Both Bands
 * ================================================================== */

/** Estimate required Elliptic filter order.
 *  L5: Uses the complete elliptic integral ratio.
 *  N >= K(k)*K'(k1) / (K'(k)*K(k1))
 *  where k = wp/ws, k1 = sqrt((10^{Ap/10}-1)/(10^{As/10}-1))
 *  K(k) = complete elliptic integral of the first kind.
 *  Reference: Cauer, Theorie der linearen Wechselstromschaltungen (1941)
 *  @param passband_ripple_db  Ap
 *  @param stopband_atten_db   As
 *  @param wp                  Passband edge (rad/s)
 *  @param ws                  Stopband edge (rad/s)
 *  @return                    Minimum required order */
int elliptic_estimate_order(double passband_ripple_db,
                             double stopband_atten_db,
                             double wp, double ws);

/** Compute the complete elliptic integral K(k) using arithmetic-geometric mean.
 *  L5: Fast, accurate algorithm. K(k) = pi/(2*AGM(1, k'))
 *  where k' = sqrt(1-k^2) is the complementary modulus.
 *  Complexity: O(log(1/eps)) iterations */
double elliptic_k(double k);

/** Design a normalized Elliptic lowpass prototype.
 *  This is the most efficient approximation — steepest transition for given order.
 *  Uses the Jacobi elliptic functions sn, cn, dn for pole/zero placement.
 *  @param order               Filter order (>= 1)
 *  @param passband_ripple_db  Maximum passband ripple
 *  @param stopband_atten_db   Minimum stopband attenuation
 *  @return                    Normalized prototype, or NULL on error */
tf_continuous_t *elliptic_prototype(int order, double passband_ripple_db,
                                     double stopband_atten_db);

/* ==================================================================
 * L5: Bessel Filter — Maximally Flat Group Delay
 * ================================================================== */

/** Compute Bessel prototype poles.
 *  L5: Derived from Bessel polynomials of the first kind.
 *  Bessel filters approximate a constant time delay in the passband.
 *  Group delay is maximally flat at w=0 — optimal phase linearity.
 *  Poles are the roots of the reverse Bessel polynomial.
 *  Reference: Thomson, "Delay Networks Having Maximally Flat Frequency
 *             Characteristics" (1949), also Storch (1954)
 *  @param order  Filter order
 *  @param poles  Output array (length >= order)
 *  @return       FILTER_OK or error code */
int bessel_poles(int order, double complex *poles);

/** Design a normalized Bessel lowpass prototype.
 *  @param order  Filter order (>= 1)
 *  @return       Normalized prototype with unit group delay at DC */
tf_continuous_t *bessel_prototype(int order);

/** Compute Bessel polynomial coefficients up to given order.
 *  L3: B_n(s) = sum_{k=0}^{n} (2n-k)!/(2^{n-k} * k! * (n-k)!) * s^k
 *  These polynomials have all roots in the left half-plane.
 *  @param n        Maximum order
 *  @param coeffs   Output array [n+1][n+1], coeffs[i][k] = coeff of s^k in B_i(s) */
int bessel_polynomials(int n, double *coeffs);

/* ==================================================================
 * L5: Frequency Transformations for Analog Prototypes
 * ================================================================== */

/** LP prototype to LP filter: s -> s/wc
 *  L2: Simple frequency scaling. */
void analog_lp_to_lp(tf_continuous_t *proto, double wc, double gain);

/** LP prototype to HP filter: s -> wc/s
 *  L2: Maps w=0 to w=inf and vice versa.
 *  Each pole at pk becomes a pole at wc/pk.
 *  Also introduces zeros at the origin equal to (order). */
tf_continuous_t *analog_lp_to_hp(const tf_continuous_t *proto,
                                  double wc, double gain);

/** LP prototype to BP filter: s -> (s^2 + w0^2)/(s*B)
 *  L2: w0 = sqrt(wL*wH) = geometric center, B = wH - wL = bandwidth.
 *  Order doubles: each pole becomes a complex pair.
 *  Also introduces zeros at the origin equal to original order. */
tf_continuous_t *analog_lp_to_bp(const tf_continuous_t *proto,
                                  double w0, double bw, double gain);

/** LP prototype to BS filter: s -> (s*B)/(s^2 + w0^2)
 *  L2: Order doubles; introduces zeros on the jw-axis at +-j*w0. */
tf_continuous_t *analog_lp_to_bs(const tf_continuous_t *proto,
                                  double w0, double bw, double gain);

/* ==================================================================
 * L5: Complete Analog Filter Design API
 * ================================================================== */

/** Design a complete analog filter from specifications.
 *  L5/L6: One-stop function that selects the appropriate prototype,
 *  estimates order if needed, computes poles/zeros, and applies
 *  frequency transformation to meet the target spec.
 *
 *  Workflow:
 *  1. Validate spec
 *  2. Estimate order if spec.order == 0
 *  3. Design normalized LP prototype with selected approximation
 *  4. Apply frequency transformation (LP/HP/BP/BS)
 *  5. Denormalize to target frequencies
 *
 *  @param spec  Filter specification
 *  @return      Complete transfer function, or NULL on error */
tf_continuous_t *analog_filter_design(const filter_spec_t *spec);

/** Compute -3dB cutoff frequency from a transfer function.
 *  Uses binary search to find frequency where |H(jw)| = |H(0)|/sqrt(2).
 *  Complexity: O(log(f_range/precision)) evaluations */
double analog_find_cutoff(const tf_continuous_t *tf, double f_start,
                           double f_stop, double precision);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_ANALOG_H */
