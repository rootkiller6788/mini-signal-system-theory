/**
 * filter_analog.c — Analog Filter Prototype Implementation
 *
 * L5: Butterworth, Chebyshev I/II, Elliptic, Bessel pole computation
 * L6: Complete analog filter design from specifications
 *
 * Key algorithms:
 *   Butterworth poles — trigonometric placement on unit circle
 *   Chebyshev I poles — hyperbolic placement on ellipse
 *   Chebyshev II poles — reciprocal Chebyshev I
 *   Elliptic K(k) — arithmetic-geometric mean (AGM) method
 *   Bessel poles — Durand-Kerner root-finding on Bessel polynomials
 *   Bessel polynomials — explicit recursion formula
 *
 * Reference:
 *   Van Valkenburg, Analog Filter Design (1982)
 *   Zverev, Handbook of Filter Synthesis (1967)
 *   Antoniou, Digital Filters: Analysis, Design, Applications (1993)
 *   Thomson, "Delay Networks Having Maximally Flat Frequency
 *            Characteristics", Proc. IEE, 1949
 */

#include "filter_analog.h"
#include "filter_design.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Butterworth Filter — Maximally Flat Magnitude
 * ================================================================== */

/**
 * butterworth_poles — L5 algorithm
 *
 * Theory: Butterworth poles lie on a circle of radius wc in the
 * left half-plane, equally spaced in angle.
 *
 * For order N, the poles are:
 *   p_k = wc * exp(j * (pi/2 + (2k+1)*pi/(2N)))  for k = 0..N-1
 *
 * The magnitude response |H(jw)| = 1/sqrt(1 + (w/wc)^{2N})
 * has its first 2N-1 derivatives zero at w=0 — maximally flat.
 *
 * Reference: S. Butterworth, "On the Theory of Filter Amplifiers",
 *            Experimental Wireless & Wireless Engineer, vol. 7, 1930
 */
int butterworth_poles(int order, double wc, double complex *poles) {
    if (!poles || order < 1 || wc <= 0.0) return FILTER_ERR_NULL_PTR;

    int k;
    for (k = 0; k < order; k++) {
        double angle = M_PI * (0.5 + (2.0 * k + 1.0) / (2.0 * order));
        poles[k] = wc * cos(angle) + wc * sin(angle) * I;
    }
    return FILTER_OK;
}

/**
 * butterworth_prototype — L5 algorithm
 *
 * Builds transfer function H(s) = 1 / prod(s - p_k).
 * For normalized prototype: wc=1, DC gain=1.
 *
 * Implements polynomial multiplication from pole pairs:
 *   (s-p)(s-p*) = s^2 - 2*Re(p)*s + |p|^2
 *
 * Complexity: O(N^2) for polynomial expansion, O(N) for pole computation.
 */
tf_continuous_t *butterworth_prototype(int order) {
    if (order < 1) return NULL;

    double complex *poles =
        (double complex *)malloc(order * sizeof(double complex));
    if (!poles) return NULL;

    int ret = butterworth_poles(order, 1.0, poles);
    if (ret != FILTER_OK) { free(poles); return NULL; }

    int den_len = order + 1;
    tf_continuous_t *tf = tf_continuous_alloc(1, den_len);
    if (!tf) { free(poles); return NULL; }

    /* Numerator: b0 = 1 for normalized LP prototype */
    tf->num[0] = 1.0;
    tf->gain = 1.0;

    /* Denominator: start with D(s) = 1, multiply by (s-p_k) */
    tf->den[0] = 1.0;
    int deg = 0;
    int k;

    for (k = 0; k < order; k++) {
        if (fabs(cimag(poles[k])) < 1e-14) {
            /* Real pole: multiply den by (s + |p|) */
            double a = -creal(poles[k]);
            int i;
            for (i = deg + 1; i >= 1; i--) {
                tf->den[i] = tf->den[i - 1] + a * tf->den[i];
            }
            tf->den[0] = a * tf->den[0];
            deg++;
        } else {
            /* Complex conjugate pair: (s-p)(s-p*) */
            double a1 = -2.0 * creal(poles[k]);
            double a0 = creal(poles[k]) * creal(poles[k])
                      + cimag(poles[k]) * cimag(poles[k]);

            double *temp = (double *)calloc(deg + 3, sizeof(double));
            if (!temp) {
                free(poles); tf_continuous_free(tf); return NULL;
            }
            int i;
            for (i = 0; i <= deg; i++) {
                temp[i]     += tf->den[i] * a0;
                temp[i + 1] += tf->den[i] * a1;
                temp[i + 2] += tf->den[i];
            }
            deg += 2;
            for (i = 0; i <= deg; i++) tf->den[i] = temp[i];
            free(temp);
            k++; /* skip conjugate */
        }
    }

    free(poles);
    return tf;
}

/**
 * butterworth_denormalize — L2 frequency scaling
 *
 * Substitutes s -> s/wc to shift cutoff from 1 to wc rad/s.
 * For D(s/wc): each coefficient d_k is divided by wc^k.
 * Complexity: O(den_len)
 */
void butterworth_denormalize(tf_continuous_t *proto, double wc, double gain) {
    if (!proto || wc <= 0.0) return;

    int i;
    double wc_pow = 1.0;
    for (i = 0; i < proto->den_len; i++) {
        proto->den[i] /= wc_pow;
        wc_pow *= wc;
    }
    if (proto->num_len > 0) {
        proto->num[0] *= gain;
    }
    proto->gain = gain;
}

/* ==================================================================
 * Chebyshev Type I Filter — Equiripple Passband
 * ================================================================== */

/**
 * chebyshev1_poles — L5 algorithm
 *
 * Poles lie on an ellipse in the left half-plane:
 *   beta = asinh(1/epsilon) / N
 *   p_k = -wc*sinh(beta)*sin(theta_k) + j*wc*cosh(beta)*cos(theta_k)
 *   where theta_k = (2k+1)*pi/(2N)
 *
 * epsilon = sqrt(10^{ripple_dB/10} - 1)
 *
 * The ellipse intersection with real axis = wc*sinh(beta)
 * The ellipse intersection with imag axis = wc*cosh(beta)
 *
 * Reference: P.L. Chebyshev, "Sur les questions de minima", 1899
 */
int chebyshev1_poles(int order, double wc, double ripple_db,
                      double complex *poles) {
    if (!poles || order < 1 || wc <= 0.0 || ripple_db < 0.0)
        return FILTER_ERR_NULL_PTR;

    double epsilon = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);
    if (epsilon < 1e-15) epsilon = 1e-15;

    double eta = asinh(1.0 / epsilon) / (double)order;
    double sinh_eta = sinh(eta);
    double cosh_eta = cosh(eta);

    int k;
    for (k = 0; k < order; k++) {
        double theta = M_PI * (2.0 * k + 1.0) / (2.0 * order);
        poles[k] = (-wc * sinh_eta * sin(theta))
                 + ( wc * cosh_eta * cos(theta)) * I;
    }
    return FILTER_OK;
}

/**
 * chebyshev1_prototype — L5 algorithm
 *
 * Builds normalized Chebyshev I transfer function.
 * DC gain: |H(0)| = 1 for odd N, |H(0)| = 1/sqrt(1+epsilon^2) for even N.
 * The gain adjustment ensures 0 dB maximum in the passband.
 */
tf_continuous_t *chebyshev1_prototype(int order, double ripple_db) {
    if (order < 1 || ripple_db < 0.0) return NULL;

    double complex *poles =
        (double complex *)malloc(order * sizeof(double complex));
    if (!poles) return NULL;

    int ret = chebyshev1_poles(order, 1.0, ripple_db, poles);
    if (ret != FILTER_OK) { free(poles); return NULL; }

    int den_len = order + 1;
    tf_continuous_t *tf = tf_continuous_alloc(1, den_len);
    if (!tf) { free(poles); return NULL; }

    tf->den[0] = 1.0;
    int deg = 0;
    int k;

    for (k = 0; k < order; k++) {
        if (fabs(cimag(poles[k])) < 1e-14) {
            double a = -creal(poles[k]);
            int i;
            for (i = deg + 1; i >= 1; i--)
                tf->den[i] = tf->den[i - 1] + a * tf->den[i];
            tf->den[0] *= a;
            deg++;
        } else {
            double a1 = -2.0 * creal(poles[k]);
            double a0 = creal(poles[k]) * creal(poles[k])
                      + cimag(poles[k]) * cimag(poles[k]);
            double *temp = (double *)calloc(deg + 3, sizeof(double));
            if (!temp) { free(poles); tf_continuous_free(tf); return NULL; }
            int i;
            for (i = 0; i <= deg; i++) {
                temp[i]     += tf->den[i] * a0;
                temp[i + 1] += tf->den[i] * a1;
                temp[i + 2] += tf->den[i];
            }
            deg += 2;
            for (i = 0; i <= deg; i++) tf->den[i] = temp[i];
            free(temp);
            k++;
        }
    }

    double epsilon = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);
    double dc_gain = (order % 2 == 0) ? 1.0 / sqrt(1.0 + epsilon * epsilon) : 1.0;
    tf->num[0] = dc_gain * tf->den[0];
    tf->gain = 1.0;

    free(poles);
    return tf;
}

/**
 * chebyshev1_estimate_order — L5 algorithm
 *
 * N >= acosh(sqrt((10^{As/10} - 1) / (10^{Ap/10} - 1))) / acosh(ws/wp)
 *
 * This is the exact formula using inverse hyperbolic cosine:
 *   acosh(x) = ln(x + sqrt(x^2 - 1)) for x >= 1.
 *
 * Complexity: O(1)
 */
int chebyshev1_estimate_order(double passband_ripple_db,
                               double stopband_atten_db,
                               double stopband_edge,
                               double passband_edge) {
    if (passband_edge <= 0.0 || stopband_edge <= passband_edge) return -1;

    double Ap = passband_ripple_db;
    double As = stopband_atten_db;
    double ratio = stopband_edge / passband_edge;

    double num = sqrt((pow(10.0, As / 10.0) - 1.0)
                    / (pow(10.0, Ap / 10.0) - 1.0));
    double N = acosh(num) / acosh(ratio);
    return (int)ceil(N);
}

/* ==================================================================
 * Chebyshev Type II (Inverse Chebyshev) — Equiripple Stopband
 * ================================================================== */

/**
 * chebyshev2_pz — L5 algorithm
 *
 * Type II poles are reciprocals of Type I poles (scaled):
 *   p_k^(II) = wc^2 / conj(p_k^(I))
 *
 * Type II zeros lie on the jw-axis (finite transmission zeros):
 *   z_k = j * wc / cos((2k+1)*pi/(2N))
 *
 * These zeros create the equiripple characteristic in the stopband
 * while the passband remains monotonic.
 */
int chebyshev2_pz(int order, double wc, double atten_db,
                   double complex *poles, double complex *zeros) {
    if (!poles || !zeros || order < 1 || wc <= 0.0 || atten_db < 0.0)
        return FILTER_ERR_NULL_PTR;

    double complex *cheb1_poles =
        (double complex *)malloc(order * sizeof(double complex));
    if (!cheb1_poles) return FILTER_ERR_MEMORY;

    int ret = chebyshev1_poles(order, 1.0, atten_db, cheb1_poles);
    if (ret != FILTER_OK) { free(cheb1_poles); return ret; }

    int k;
    for (k = 0; k < order; k++) {
        double complex p_conj = conj(cheb1_poles[k]);
        if (cabs(p_conj) < 1e-15) {
            poles[k] = -1e10 + 0 * I;
        } else {
            poles[k] = (wc * wc) / p_conj;
        }
    }

    for (k = 0; k < order; k++) {
        double cos_val = cos(M_PI * (2.0 * k + 1.0) / (2.0 * order));
        if (fabs(cos_val) < 1e-15) {
            zeros[k] = 1e10 * I;
        } else {
            zeros[k] = (wc / cos_val) * I;
        }
    }

    free(cheb1_poles);
    return FILTER_OK;
}

/**
 * chebyshev2_prototype — L5 algorithm
 *
 * Builds normalized (wc=1) Chebyshev II prototype.
 * The numerator contains finite zeros creating stopband notches.
 */
tf_continuous_t *chebyshev2_prototype(int order, double atten_db) {
    if (order < 1 || atten_db < 0.0) return NULL;

    double complex *poles =
        (double complex *)malloc(order * sizeof(double complex));
    double complex *zeros =
        (double complex *)malloc(order * sizeof(double complex));
    if (!poles || !zeros) {
        free(poles); free(zeros); return NULL;
    }

    int ret = chebyshev2_pz(order, 1.0, atten_db, poles, zeros);
    if (ret != FILTER_OK) { free(poles); free(zeros); return NULL; }

    int num_len = order + 1;
    int den_len = order + 1;
    tf_continuous_t *tf = tf_continuous_alloc(num_len, den_len);
    if (!tf) { free(poles); free(zeros); return NULL; }

    /* Build numerator from zeros */
    tf->num[0] = 1.0;
    int num_deg = 0;
    int k;
    for (k = 0; k < order; k++) {
        /* Type II zeros are purely imaginary, in conjugate pairs */
        if (cimag(zeros[k]) > 1e-9 || cimag(zeros[k]) < -1e-9) {
            /* Zero at +-j*wz: (s - j*wz)(s + j*wz) = s^2 + wz^2 */
            double wz2 = cimag(zeros[k]) * cimag(zeros[k]);
            double *temp = (double *)calloc(num_deg + 3, sizeof(double));
            if (!temp) { free(poles); free(zeros); tf_continuous_free(tf); return NULL; }
            int i;
            for (i = 0; i <= num_deg; i++) {
                temp[i]     += tf->num[i] * wz2;
                temp[i + 2] += tf->num[i];
            }
            num_deg += 2;
            for (i = 0; i <= num_deg; i++) tf->num[i] = temp[i];
            free(temp);
            k++;
        }
    }

    /* Build denominator from poles */
    tf->den[0] = 1.0;
    int den_deg = 0;
    for (k = 0; k < order; k++) {
        if (fabs(cimag(poles[k])) < 1e-14) {
            double a = -creal(poles[k]);
            int i;
            for (i = den_deg + 1; i >= 1; i--)
                tf->den[i] = tf->den[i - 1] + a * tf->den[i];
            tf->den[0] *= a;
            den_deg++;
        } else {
            double a1 = -2.0 * creal(poles[k]);
            double a0 = creal(poles[k]) * creal(poles[k])
                      + cimag(poles[k]) * cimag(poles[k]);
            double *temp = (double *)calloc(den_deg + 3, sizeof(double));
            if (!temp) { free(poles); free(zeros); tf_continuous_free(tf); return NULL; }
            int i;
            for (i = 0; i <= den_deg; i++) {
                temp[i]     += tf->den[i] * a0;
                temp[i + 1] += tf->den[i] * a1;
                temp[i + 2] += tf->den[i];
            }
            den_deg += 2;
            for (i = 0; i <= den_deg; i++) tf->den[i] = temp[i];
            free(temp);
            k++;
        }
    }

    /* Normalize: |H(0)| = 1 → num[0] = den[0] */
    if (fabs(tf->den[0]) > 1e-15) {
        tf->num[0] = tf->den[0];
        tf->gain = 1.0;
    }

    free(poles);
    free(zeros);
    return tf;
}

/* ==================================================================
 * Elliptic / Cauer Filter — Equiripple Both Bands
 * ================================================================== */

/**
 * elliptic_k — L5 algorithm
 *
 * Complete elliptic integral of the first kind K(k):
 *   K(k) = pi / (2 * AGM(1, sqrt(1-k^2)))
 *
 * The Arithmetic-Geometric Mean (AGM) converges quadratically:
 *   a_{n+1} = (a_n + g_n) / 2
 *   g_{n+1} = sqrt(a_n * g_n)
 *
 * Typically 5-8 iterations suffice for double precision.
 *
 * Reference: M. Abramowitz & I.A. Stegun, Handbook of Mathematical
 *            Functions, Sec. 17.6, Dover, 1965.
 * Complexity: O(log(1/eps)) per evaluation
 */
double elliptic_k(double k) {
    if (k < 0.0 || k > 1.0) return -1.0;
    if (k == 1.0) return INFINITY;
    if (k == 0.0) return M_PI / 2.0;

    double a = 1.0;
    double g = sqrt(1.0 - k * k);
    int i;
    for (i = 0; i < 50; i++) {
        double a_next = (a + g) / 2.0;
        double g_next = sqrt(a * g);
        if (fabs(a - g) < 1e-15) break;
        a = a_next;
        g = g_next;
    }
    return M_PI / (2.0 * a);
}

/**
 * elliptic_estimate_order — L5 algorithm
 *
 * Uses the ratio of complete elliptic integrals:
 *   N >= K(k) * K'(k1) / (K'(k) * K(k1))
 * where:
 *   k  = wp/ws  (discrimination factor)
 *   k1 = sqrt((10^{Ap/10} - 1) / (10^{As/10} - 1))  (selectivity factor)
 *   K'(x) = K(sqrt(1 - x^2))
 *
 * This gives the theoretically minimum order — no filter of lower
 * order can meet the given specifications with the same ripple constraints.
 */
int elliptic_estimate_order(double passband_ripple_db,
                             double stopband_atten_db,
                             double wp, double ws) {
    if (wp <= 0.0 || ws <= wp) return -1;

    double k  = wp / ws;
    double k1_num = pow(10.0, passband_ripple_db / 10.0) - 1.0;
    double k1_den = pow(10.0, stopband_atten_db / 10.0) - 1.0;
    if (k1_den < 1e-15) return -1;
    double k1 = sqrt(k1_num / k1_den);

    double Kk   = elliptic_k(k);
    double Kkp  = elliptic_k(sqrt(1.0 - k * k));
    double Kk1  = elliptic_k(k1);
    double Kk1p = elliptic_k(sqrt(1.0 - k1 * k1));

    double N = (Kk * Kk1p) / (Kkp * Kk1);
    return (int)ceil(N);
}

/**
 * elliptic_prototype — L5 algorithm
 *
 * Full elliptic design requires Jacobi elliptic functions (sn, cn, dn)
 * for precise pole/zero placement. This implementation delegates to
 * Chebyshev I as a practical approximation for low orders.
 *
 * For production-quality elliptic filters, use established libraries
 * (e.g., scipy.signal.ellip, MATLAB ellip).
 */
tf_continuous_t *elliptic_prototype(int order, double passband_ripple_db,
                                     double stopband_atten_db) {
    if (order < 1 || order > 10) return NULL;
    if (passband_ripple_db < 0.0 || stopband_atten_db < 0.0) return NULL;
    return chebyshev1_prototype(order, passband_ripple_db);
}

/* ==================================================================
 * Bessel Filter — Maximally Flat Group Delay
 * ================================================================== */

/**
 * bessel_polynomials — L3/L5 algorithm
 *
 * Bessel polynomials B_n(s) via recurrence:
 *   B_0(s) = 1
 *   B_1(s) = s + 1
 *   B_n(s) = (2n - 1) * B_{n-1}(s) + s^2 * B_{n-2}(s)
 *
 * Explicit form:
 *   B_n(s) = sum_{k=0}^{n} (2n-k)! / (2^{n-k} * k! * (n-k)!) * s^k
 *
 * The Bessel filter denominator is the reverse Bessel polynomial:
 *   theta_n(s) = B_n(1/s) * s^n = sum_{k=0}^{n} B_n[k] * s^{n-k}
 *
 * Reference: W.E. Thomson, "Delay Networks Having Maximally Flat
 *            Frequency Characteristics", Proc. IEE, Part C, 1949.
 */
int bessel_polynomials(int n, double *coeffs) {
    if (!coeffs || n < 0 || n > 20) {
        return (n > 20) ? FILTER_ERR_INVALID_ORDER : FILTER_ERR_NULL_PTR;
    }

    int stride = n + 1;
    int i, k;

    /* Zero-initialize */
    for (i = 0; i <= n; i++) {
        for (k = 0; k <= n; k++) {
            coeffs[i * stride + k] = 0.0;
        }
    }

    /* B_0(s) = 1 */
    coeffs[0 * stride + 0] = 1.0;

    if (n >= 1) {
        /* B_1(s) = s + 1 */
        coeffs[1 * stride + 0] = 1.0;
        coeffs[1 * stride + 1] = 1.0;
    }

    /* Recurrence: B_n(s) = (2n-1)*B_{n-1}(s) + s^2*B_{n-2}(s) */
    for (i = 2; i <= n; i++) {
        double factor = 2.0 * i - 1.0;
        for (k = 0; k <= i; k++) {
            double val = 0.0;
            if (k <= i - 1)
                val += factor * coeffs[(i - 1) * stride + k];
            if (k >= 2 && (k - 2) <= (i - 2))
                val += coeffs[(i - 2) * stride + (k - 2)];
            coeffs[i * stride + k] = val;
        }
    }
    return FILTER_OK;
}

/**
 * bessel_poles — L5 algorithm
 *
 * Finds roots of the reverse Bessel polynomial using the
 * Durand-Kerner simultaneous root-finding method (Weierstrass method).
 *
 * Initial guesses: equally spaced on a circle.
 * Iteration: p_i^{new} = p_i - P(p_i) / prod_{j!=i} (p_i - p_j)
 *
 * Reference: E. Durand, "Solutions Numeriques des Equations
 *            Algebriques", 1960; I.O. Kerner, Numer. Math., 1966.
 */
int bessel_poles(int order, double complex *poles) {
    if (!poles || order < 1 || order > 20) {
        return (order > 20) ? FILTER_ERR_INVALID_ORDER : FILTER_ERR_NULL_PTR;
    }

    int stride = order + 1;
    double *poly_coeffs =
        (double *)malloc((order + 1) * stride * sizeof(double));
    if (!poly_coeffs) return FILTER_ERR_MEMORY;

    int ret = bessel_polynomials(order, poly_coeffs);
    if (ret != FILTER_OK) { free(poly_coeffs); return ret; }

    /* Reverse the polynomial: theta_n(s) = sum B_n[k] * s^{n-k} */
    double *coeff = (double *)malloc((order + 1) * sizeof(double));
    if (!coeff) { free(poly_coeffs); return FILTER_ERR_MEMORY; }
    int k;
    for (k = 0; k <= order; k++) {
        coeff[order - k] = poly_coeffs[order * stride + k];
    }
    free(poly_coeffs);

    /* Normalize leading coefficient to 1 */
    if (fabs(coeff[0]) < 1e-15) { free(coeff); return FILTER_ERR_NUMERICAL; }
    for (k = 1; k <= order; k++) coeff[k] /= coeff[0];
    coeff[0] = 1.0;

    /* Initial guesses on a circle */
    int i;
    for (i = 0; i < order; i++) {
        double angle = M_PI * (0.5 + (2.0 * i + 1.0) / (2.0 * order));
        poles[i] = cos(angle) + sin(angle) * I;
    }

    /* Durand-Kerner iteration */
    int iter, max_iter = 100;
    for (iter = 0; iter < max_iter; iter++) {
        int converged = 1;
        for (i = 0; i < order; i++) {
            /* Evaluate polynomial at p_i */
            double complex p_val = 0;
            int j;
            double complex p_pow = 1;
            for (j = 0; j <= order; j++) {
                p_val += coeff[j] * p_pow;
                p_pow *= poles[i];
            }

            /* Compute product (p_i - p_j) for j != i */
            double complex prod = 1;
            for (j = 0; j < order; j++) {
                if (j != i) prod *= (poles[i] - poles[j]);
            }
            if (cabs(prod) < 1e-15) prod = 1e-10;

            double complex delta = p_val / prod;
            poles[i] -= delta;
            if (cabs(delta) > 1e-12) converged = 0;
        }
        if (converged) break;
    }

    free(coeff);
    return FILTER_OK;
}

/**
 * bessel_prototype — L5 algorithm
 *
 * Builds normalized Bessel lowpass transfer function.
 * The Bessel filter is unique in having maximally flat group delay
 * at omega=0, making it optimal for preserving pulse/waveform shape
 * in the time domain (minimal overshoot, no ringing).
 *
 * Trade-off: requires higher order than Butterworth for equivalent
 * stopband attenuation.
 */
tf_continuous_t *bessel_prototype(int order) {
    if (order < 1 || order > 20) return NULL;

    double complex *poles =
        (double complex *)malloc(order * sizeof(double complex));
    if (!poles) return NULL;

    int ret = bessel_poles(order, poles);
    if (ret != FILTER_OK) { free(poles); return NULL; }

    int den_len = order + 1;
    tf_continuous_t *tf = tf_continuous_alloc(1, den_len);
    if (!tf) { free(poles); return NULL; }

    tf->num[0] = 1.0;
    tf->gain = 1.0;
    tf->den[0] = 1.0;

    int deg = 0;
    int k;
    for (k = 0; k < order; k++) {
        if (fabs(cimag(poles[k])) < 1e-12) {
            double a = -creal(poles[k]);
            int i;
            for (i = deg + 1; i >= 1; i--)
                tf->den[i] = tf->den[i - 1] + a * tf->den[i];
            tf->den[0] *= a;
            deg++;
        } else {
            double a1 = -2.0 * creal(poles[k]);
            double a0 = creal(poles[k]) * creal(poles[k])
                      + cimag(poles[k]) * cimag(poles[k]);
            double *temp = (double *)calloc(deg + 3, sizeof(double));
            if (!temp) { free(poles); tf_continuous_free(tf); return NULL; }
            int i;
            for (i = 0; i <= deg; i++) {
                temp[i]     += tf->den[i] * a0;
                temp[i + 1] += tf->den[i] * a1;
                temp[i + 2] += tf->den[i];
            }
            deg += 2;
            for (i = 0; i <= deg; i++) tf->den[i] = temp[i];
            free(temp);
            k++;
        }
    }

    tf->num[0] = tf->den[0]; /* DC gain = 1 */
    free(poles);
    return tf;
}

/* ==================================================================
 * Frequency Transformations
 * ================================================================== */

void analog_lp_to_lp(tf_continuous_t *proto, double wc, double gain) {
    if (!proto || wc <= 0.0) return;

    double wc_pow = 1.0;
    int i;
    for (i = 0; i < proto->den_len; i++) {
        proto->den[i] /= wc_pow;
        wc_pow *= wc;
    }
    wc_pow = 1.0;
    for (i = 0; i < proto->num_len; i++) {
        proto->num[i] /= wc_pow;
        wc_pow *= wc;
    }
    proto->gain *= gain;
}

tf_continuous_t *analog_lp_to_hp(const tf_continuous_t *proto,
                                  double wc, double gain) {
    if (!proto || wc <= 0.0) return NULL;

    int N = proto->den_len - 1;
    int new_len = N + 1;
    tf_continuous_t *result = tf_continuous_alloc(new_len, new_len);
    if (!result) return NULL;

    /* Numerator of HP = denominator of LP, reversed and scaled */
    int k;
    for (k = 0; k <= N; k++) {
        result->num[N - k] = proto->den[k] * pow(wc, k);
    }
    /* Denominator of HP = numerator of LP, reversed and scaled */
    int M = proto->num_len - 1;
    for (k = 0; k <= M; k++) {
        result->den[N - k] = proto->num[k] * pow(wc, k);
    }
    for (k = M + 1; k <= N; k++) {
        result->den[N - k] = 0.0;
    }
    result->gain = proto->gain * gain;
    return result;
}

/**
 * analog_lp_to_bp — L2/L5 LP prototype to bandpass transformation
 *
 * Substitution: s -> (s^2 + w0^2) / (s * B)
 *
 * This doubles the filter order: an N-th order LP becomes a 2N-th order BP.
 * w0 = sqrt(wL * wH) is the geometric center frequency.
 * B = wH - wL is the bandwidth.
 *
 * Method: factor the LP prototype into analog biquad sections,
 * apply the s-domain substitution to each section, then multiply.
 *
 * The substitution maps each LP pole/zero as follows:
 *   LP pole at s = -sigma + j*w_k  →  BP poles at roots of s^2 - p*s + w0^2 = 0
 *   where p = -sigma + j*w_k.
 *
 * Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 10
 * Complexity: O(N^2)
 */
tf_continuous_t *analog_lp_to_bp(const tf_continuous_t *proto,
                                  double w0, double bw, double gain) {
    if (!proto || w0 <= 0.0 || bw <= 0.0) return NULL;

    int N = proto->den_len - 1;
    if (N < 1 || N > 12) return NULL;

    /* Factor into biquads */
    int max_bq = (N + 1) / 2;
    analog_biquad_t *bqs =
        (analog_biquad_t *)calloc(max_bq, sizeof(analog_biquad_t));
    if (!bqs) return NULL;

    int num_bq = tf_continuous_to_biquads(proto, bqs, max_bq);
    if (num_bq < 1) { free(bqs); return NULL; }

    /* Build result from cascaded sections */
    tf_continuous_t *result = tf_continuous_alloc(1, 1);
    if (!result) { free(bqs); return NULL; }
    result->num[0] = 1.0;
    result->den[0] = 1.0;
    result->gain = 1.0;

    int k;
    for (k = 0; k < num_bq; k++) {
        double B2 = bw * bw;
        double w02 = w0 * w0;
        double w04 = w02 * w02;

        /* Each LP biquad: H_LP(s) = (b2 s^2 + b1 s + b0) / (s^2 + a1 s + a0)
         * After s -> (s^2 + w0^2)/(s*B):
         *   Numerator of BP:  b2*(s^2+w0^2)^2 + b1*B*s*(s^2+w0^2) + b0*B^2*s^2  (divided by B^2*s^2)
         *   Denominator of BP: (s^2+w0^2)^2 + a1*B*s*(s^2+w0^2) + a0*B^2*s^2  (divided by B^2*s^2)
         *
         * Multiply through by B^2 to clear denominator:
         */
        tf_continuous_t *section = tf_continuous_alloc(5, 5);
        if (!section) { tf_continuous_free(result); free(bqs); return NULL; }

        /* Numerator (degree 4) */
        section->num[4] = bqs[k].b2;
        section->num[3] = bqs[k].b1 * bw;
        section->num[2] = 2.0 * bqs[k].b2 * w02 + bqs[k].b0 * B2;
        section->num[1] = bqs[k].b1 * bw * w02;
        section->num[0] = bqs[k].b2 * w04;

        /* Denominator (degree 4, a2 = 1 for normalized biquad) */
        section->den[4] = 1.0;
        section->den[3] = bqs[k].a1 * bw;
        section->den[2] = 2.0 * w02 + bqs[k].a0 * B2;
        section->den[1] = bqs[k].a1 * bw * w02;
        section->den[0] = w04;

        section->gain = 1.0;

        tf_continuous_t *new_result = tf_continuous_multiply(result, section);
        tf_continuous_free(result);
        tf_continuous_free(section);
        result = new_result;
        if (!result) { free(bqs); return NULL; }
    }

    /* Apply overall gain */
    result->gain = proto->gain * gain;

    /* Normalize numerator: for BP, DC gain is 0, so we match peak gain */
    int i;
    for (i = 0; i < result->num_len; i++) result->num[i] *= proto->gain * gain;

    free(bqs);
    return result;
}

/**
 * analog_lp_to_bs — L2/L5 LP prototype to bandstop transformation
 *
 * Substitution: s -> (s * B) / (s^2 + w0^2)
 *
 * Dual of LP→BP: maps the LP passband to two frequency bands
 * (below wL and above wH), with the stopband centered at w0.
 *
 * Like LP→BP, this doubles the order: N → 2N.
 * DC and infinite frequency map to the passbands.
 *
 * Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 10
 *             Zverev, Handbook of Filter Synthesis (1967)
 * Complexity: O(N^2)
 */
tf_continuous_t *analog_lp_to_bs(const tf_continuous_t *proto,
                                  double w0, double bw, double gain) {
    if (!proto || w0 <= 0.0 || bw <= 0.0) return NULL;

    int N = proto->den_len - 1;
    if (N < 1 || N > 12) return NULL;

    int max_bq = (N + 1) / 2;
    analog_biquad_t *bqs =
        (analog_biquad_t *)calloc(max_bq, sizeof(analog_biquad_t));
    if (!bqs) return NULL;

    int num_bq = tf_continuous_to_biquads(proto, bqs, max_bq);
    if (num_bq < 1) { free(bqs); return NULL; }

    tf_continuous_t *result = tf_continuous_alloc(1, 1);
    if (!result) { free(bqs); return NULL; }
    result->num[0] = 1.0;
    result->den[0] = 1.0;
    result->gain = 1.0;

    int k;
    for (k = 0; k < num_bq; k++) {
        double B = bw;
        double B2 = B * B;
        double w02 = w0 * w0;

        /* For each LP biquad: H_LP(s) = (b2 s^2 + b1 s + b0)/(s^2 + a1 s + a0)
         *
         * Apply s -> (s*B)/(s^2 + w0^2):
         *   Numerator: b2*B^2*s^2 + b1*B*s*(s^2+w0^2) + b0*(s^2+w0^2)^2
         *   Denominator: B^2*s^2 + a1*B*s*(s^2+w0^2) + a0*(s^2+w0^2)^2
         *
         * After expansion (both degree 4):
         */
        tf_continuous_t *section = tf_continuous_alloc(5, 5);
        if (!section) { tf_continuous_free(result); free(bqs); return NULL; }

        /* Numerator */
        section->num[4] = bqs[k].b0;
        section->num[3] = bqs[k].b1 * B;
        section->num[2] = 2.0 * bqs[k].b0 * w02 + bqs[k].b2 * B2;
        section->num[1] = bqs[k].b1 * B * w02;
        section->num[0] = bqs[k].b0 * w02 * w02;

        /* Denominator */
        section->den[4] = bqs[k].a0;
        section->den[3] = bqs[k].a1 * B;
        section->den[2] = 2.0 * bqs[k].a0 * w02 + B2;
        section->den[1] = bqs[k].a1 * B * w02;
        section->den[0] = bqs[k].a0 * w02 * w02;

        section->gain = 1.0;

        tf_continuous_t *new_result = tf_continuous_multiply(result, section);
        tf_continuous_free(result);
        tf_continuous_free(section);
        result = new_result;
        if (!result) { free(bqs); return NULL; }
    }

    result->gain = proto->gain * gain;
    int i;
    for (i = 0; i < result->num_len; i++) result->num[i] *= proto->gain * gain;

    free(bqs);
    return result;
}

/* ==================================================================
 * Complete Analog Filter Design
 * ================================================================== */

/**
 * analog_filter_design — L6 canonical problem
 *
 * One-stop analog filter design pipeline:
 * 1. Validate specification
 * 2. Estimate order (if spec->order == 0)
 * 3. Design normalized LP prototype
 * 4. Apply frequency transformation
 * 5. Adjust gain
 *
 * Works for all five approximation families (Butterworth, Chebyshev I/II,
 * Elliptic, Bessel) and four filter types (LP, HP, BP, BS).
 */
tf_continuous_t *analog_filter_design(const filter_spec_t *spec) {
    if (!spec || !filter_spec_is_valid(spec)) return NULL;

    int order = spec->order;
    double wc = 2.0 * M_PI * spec->fc1;

    if (order <= 0) {
        double wp = 2.0 * M_PI * spec->passband_edge;
        double ws = 2.0 * M_PI * spec->stopband_edge;
        if (wp <= 0.0) wp = wc * 0.9;
        if (ws <= 0.0) ws = wc * 1.5;

        order = butterworth_estimate_order(
            spec->passband_ripple_db, spec->stopband_atten_db, wp, ws);
        if (order < 1) order = 1;
        if (order > 20) return NULL;
    }

    /* Design prototype */
    tf_continuous_t *proto = NULL;
    switch (spec->approx) {
        case APPROX_BUTTERWORTH:
            proto = butterworth_prototype(order); break;
        case APPROX_CHEBYSHEV_I:
            proto = chebyshev1_prototype(order, spec->passband_ripple_db); break;
        case APPROX_CHEBYSHEV_II:
            proto = chebyshev2_prototype(order, spec->stopband_atten_db); break;
        case APPROX_ELLIPTIC:
            proto = elliptic_prototype(order, spec->passband_ripple_db,
                                        spec->stopband_atten_db); break;
        case APPROX_BESSEL:
            proto = bessel_prototype(order); break;
        default:
            proto = butterworth_prototype(order); break;
    }

    if (!proto) return NULL;

    double linear_gain = pow(10.0, spec->gain_db / 20.0);
    tf_continuous_t *result = NULL;

    switch (spec->type) {
        case FILTER_LOWPASS:
            analog_lp_to_lp(proto, wc, linear_gain);
            result = proto;
            break;
        case FILTER_HIGHPASS:
            result = analog_lp_to_hp(proto, wc, linear_gain);
            tf_continuous_free(proto);
            break;
        case FILTER_BANDPASS:
            result = analog_lp_to_bp(proto,
                2.0 * M_PI * sqrt(spec->fc1 * spec->fc2),
                2.0 * M_PI * (spec->fc2 - spec->fc1), linear_gain);
            tf_continuous_free(proto);
            break;
        case FILTER_BANDSTOP:
            result = analog_lp_to_bs(proto,
                2.0 * M_PI * sqrt(spec->fc1 * spec->fc2),
                2.0 * M_PI * (spec->fc2 - spec->fc1), linear_gain);
            tf_continuous_free(proto);
            break;
        default:
            analog_lp_to_lp(proto, wc, linear_gain);
            result = proto;
            break;
    }

    return result;
}

double analog_find_cutoff(const tf_continuous_t *tf, double f_start,
                           double f_stop, double precision) {
    if (!tf || f_start < 0.0 || f_stop <= f_start) return -1.0;

    double mag_dc = cabs(tf_continuous_eval(tf, 0));
    if (mag_dc < 1e-15) return -1.0;
    double target = mag_dc / sqrt(2.0);

    double lo = f_start, hi = f_stop;
    int iter;
    for (iter = 0; iter < 100; iter++) {
        double mid = (lo + hi) / 2.0;
        double w = 2.0 * M_PI * mid;
        double mag = cabs(tf_continuous_eval(tf, w * I));
        if (fabs(mag - target) < precision) return mid;
        if (mag > target) lo = mid;
        else hi = mid;
        if ((hi - lo) < precision) return (lo + hi) / 2.0;
    }
    return (lo + hi) / 2.0;
}

int butterworth_estimate_order(double passband_ripple_db,
                                double stopband_atten_db,
                                double wp, double ws) {
    if (wp <= 0.0 || ws <= wp) return -1;
    double Ap = passband_ripple_db, As = stopband_atten_db;
    double ratio = ws / wp;
    double num = (pow(10.0, As / 10.0) - 1.0)
               / (pow(10.0, Ap / 10.0) - 1.0);
    double N = log10(num) / (2.0 * log10(ratio));
    return (int)ceil(N);
}

int bessel_estimate_order(double passband_ripple_db,
                           double stopband_atten_db,
                           double wp, double ws) {
    int bw_order = butterworth_estimate_order(
        passband_ripple_db, stopband_atten_db, wp, ws);
    if (bw_order < 1) return 1;
    return (int)ceil(bw_order * 1.3);
}
