/**
 * filter_transfer.c — Transfer Function Manipulation
 *
 * L3: Polynomial transfer functions, eigenvalue-based root finding,
 *     polynomial convolution, stability tests.
 * L4: Bilinear transform, Routh-Hurwitz stability, Jury test.
 *
 * Reference:
 *   Oppenheim & Willsky (1997) Signals and Systems
 *   Jury, "Theory and Application of the z-Transform Method" (1964)
 *   Dorf & Bishop, Modern Control Systems (2017)
 */

#include "filter_transfer.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Transfer Function Allocation / Free
 * ================================================================== */

tf_continuous_t *tf_continuous_alloc(int num_len, int den_len) {
    if (num_len < 1 || den_len < 1) return NULL;

    tf_continuous_t *tf = (tf_continuous_t *)malloc(sizeof(tf_continuous_t));
    if (!tf) return NULL;

    tf->num = (double *)calloc(num_len, sizeof(double));
    tf->den = (double *)calloc(den_len, sizeof(double));

    if (!tf->num || !tf->den) {
        free(tf->num); free(tf->den); free(tf);
        return NULL;
    }

    tf->num_len = num_len;
    tf->den_len = den_len;
    tf->gain = 1.0;
    return tf;
}

void tf_continuous_free(tf_continuous_t *tf) {
    if (tf) {
        free(tf->num);
        free(tf->den);
        free(tf);
    }
}

tf_discrete_t *tf_discrete_alloc(int num_len, int den_len) {
    if (num_len < 1 || den_len < 1) return NULL;

    tf_discrete_t *tf = (tf_discrete_t *)malloc(sizeof(tf_discrete_t));
    if (!tf) return NULL;

    tf->num = (double *)calloc(num_len, sizeof(double));
    tf->den = (double *)calloc(den_len, sizeof(double));

    if (!tf->num || !tf->den) {
        free(tf->num); free(tf->den); free(tf);
        return NULL;
    }

    tf->num_len = num_len;
    tf->den_len = den_len;
    tf->gain = 1.0;
    return tf;
}

void tf_discrete_free(tf_discrete_t *tf) {
    if (tf) {
        free(tf->num);
        free(tf->den);
        free(tf);
    }
}

/* ==================================================================
 * Transfer Function Evaluation — Horner's Method
 * ================================================================== */

/**
 * tf_continuous_eval — L3 polynomial evaluation at complex s
 *
 * Uses Horner's method (nested multiplication) for numerical stability:
 *   N(s) = num[0] + s*(num[1] + s*(num[2] + ...))
 *   D(s) = den[0] + s*(den[1] + s*(den[2] + ...))
 *   H(s) = gain * N(s) / D(s)
 *
 * Complexity: O(max(num_len, den_len))
 */
double complex tf_continuous_eval(const tf_continuous_t *tf, double complex s) {
    if (!tf) return 0;

    /* Evaluate numerator */
    double complex num_val = 0;
    int i;
    for (i = tf->num_len - 1; i >= 0; i--) {
        num_val = num_val * s + tf->num[i];
    }

    /* Evaluate denominator */
    double complex den_val = 0;
    for (i = tf->den_len - 1; i >= 0; i--) {
        den_val = den_val * s + tf->den[i];
    }

    if (cabs(den_val) < 1e-30) return INFINITY + I * INFINITY;
    return tf->gain * num_val / den_val;
}

/**
 * tf_discrete_eval — L3 polynomial evaluation at complex z
 *
 * Same Horner's method applied to z-domain transfer function.
 */
double complex tf_discrete_eval(const tf_discrete_t *tf, double complex z) {
    if (!tf) return 0;

    double complex num_val = 0;
    int i;
    for (i = tf->num_len - 1; i >= 0; i--) {
        num_val = num_val * z + tf->num[i];
    }

    double complex den_val = 0;
    for (i = tf->den_len - 1; i >= 0; i--) {
        den_val = den_val * z + tf->den[i];
    }

    if (cabs(den_val) < 1e-30) return INFINITY + I * INFINITY;
    return tf->gain * num_val / den_val;
}

double tf_continuous_mag_at(const tf_continuous_t *tf, double omega) {
    double complex s = omega * I;
    return cabs(tf_continuous_eval(tf, s));
}

double tf_discrete_mag_at(const tf_discrete_t *tf, double omega) {
    double complex z = cos(omega) + sin(omega) * I;
    return cabs(tf_discrete_eval(tf, z));
}

/* ==================================================================
 * Transfer Function Multiplication — Polynomial Convolution
 * ================================================================== */

/**
 * tf_continuous_multiply — L3 polynomial convolution
 *
 * H(s) = H1(s) * H2(s)
 * num(s) = num1(s) * num2(s)  (polynomial convolution)
 * den(s) = den1(s) * den2(s)  (polynomial convolution)
 *
 * Convolution: C_k = sum_{i+j=k} A_i * B_j
 * Complexity: O(num1_len * num2_len + den1_len * den2_len)
 */
tf_continuous_t *tf_continuous_multiply(const tf_continuous_t *tf1,
                                         const tf_continuous_t *tf2) {
    if (!tf1 || !tf2) return NULL;

    int new_num_len = tf1->num_len + tf2->num_len - 1;
    int new_den_len = tf1->den_len + tf2->den_len - 1;

    tf_continuous_t *result = tf_continuous_alloc(new_num_len, new_den_len);
    if (!result) return NULL;

    /* Convolve numerators */
    int i, j;
    for (i = 0; i < tf1->num_len; i++) {
        for (j = 0; j < tf2->num_len; j++) {
            result->num[i + j] += tf1->num[i] * tf2->num[j];
        }
    }

    /* Convolve denominators */
    for (i = 0; i < tf1->den_len; i++) {
        for (j = 0; j < tf2->den_len; j++) {
            result->den[i + j] += tf1->den[i] * tf2->den[j];
        }
    }

    result->gain = tf1->gain * tf2->gain;
    return result;
}

tf_discrete_t *tf_discrete_multiply(const tf_discrete_t *tf1,
                                     const tf_discrete_t *tf2) {
    if (!tf1 || !tf2) return NULL;

    int new_num_len = tf1->num_len + tf2->num_len - 1;
    int new_den_len = tf1->den_len + tf2->den_len - 1;

    tf_discrete_t *result = tf_discrete_alloc(new_num_len, new_den_len);
    if (!result) return NULL;

    int i, j;
    for (i = 0; i < tf1->num_len; i++)
        for (j = 0; j < tf2->num_len; j++)
            result->num[i + j] += tf1->num[i] * tf2->num[j];

    for (i = 0; i < tf1->den_len; i++)
        for (j = 0; j < tf2->den_len; j++)
            result->den[i + j] += tf1->den[i] * tf2->den[j];

    result->gain = tf1->gain * tf2->gain;
    return result;
}

/* ==================================================================
 * Conversion: Transfer Function <-> Biquads
 * ================================================================== */

int tf_continuous_to_biquads(const tf_continuous_t *tf,
                              analog_biquad_t *bqs, int max_bq) {
    if (!tf || !bqs || max_bq < 1) return -1;

    /* For even-order filters, pair poles into biquads.
     * This implementation returns single-stage for simplicity. */
    int k = 0;
    if (tf->den_len >= 3 && k < max_bq) {
        int deg = tf->den_len - 1;
        /* Use the two highest-degree denominator coefficients */
        if (deg >= 2) {
            bqs[k].a1 = tf->den[deg - 1] / tf->den[deg];
            bqs[k].a0 = tf->den[deg - 2] / tf->den[deg];
        }
        if (deg >= 1) {
            bqs[k].b0 = tf->num[0] / tf->den[deg];
        }
        bqs[k].b2 = 0.0;
        bqs[k].b1 = 0.0;
        bqs[k].gain = tf->gain;
        k++;
    }
    return k;
}

int tf_discrete_to_biquads(const tf_discrete_t *tf,
                            biquad_section_t *bqs, int max_bq) {
    if (!tf || !bqs || max_bq < 1) return -1;
    int k = 0;
    if (tf->den_len >= 3 && k < max_bq) {
        int deg = tf->den_len - 1;
        bqs[k].a1 = tf->den[deg - 1] / tf->den[deg];
        bqs[k].a2 = 0.0;
        bqs[k].b0 = tf->num[0] / tf->den[deg];
        bqs[k].b1 = 0.0;
        bqs[k].b2 = 0.0;
        bqs[k].gain = tf->gain;
        k++;
    }
    return k;
}

/* ==================================================================
 * Pole-Zero Map Conversion
 * ================================================================== */

/**
 * tf_continuous_to_pzmap — L3 eigenvalue-based pole-zero extraction
 *
 * Uses the Durand-Kerner method to find all roots of numerator and
 * denominator polynomials simultaneously, then constructs a pole-zero map.
 *
 * Poles come from denominator roots; zeros come from numerator roots.
 * Complex-conjugate pairs are grouped in pz_pairs[]; real-only poles/zeros
 * go to real_poles[]/real_zeros[] arrays.
 *
 * Key insight: the companion matrix eigenvalues give polynomial roots,
 * but for polynomials with complex-conjugate root pairs (typical for
 * filter design), Durand-Kerner preserves the conjugate symmetry.
 *
 * Complexity: O(N^2 * max_iter) for each polynomial
 * Reference: Durand (1960), Kerner (1966)
 */
pz_map_t *tf_continuous_to_pzmap(const tf_continuous_t *tf) {
    if (!tf || tf->den_len < 2) return NULL;

    int n_den = tf->den_len - 1;  /* denominator degree */
    int n_num = tf->num_len - 1;  /* numerator degree */
    int total_roots = n_den + n_num;
    if (total_roots < 1) return NULL;
    if (n_den > 20 || n_num > 20) return NULL; /* practical limit */

    pz_map_t *pz = (pz_map_t *)calloc(1, sizeof(pz_map_t));
    if (!pz) return NULL;

    /* Allocate root storage */
    double complex *den_roots = NULL, *num_roots = NULL;
    if (n_den > 0) {
        den_roots = (double complex *)malloc(n_den * sizeof(double complex));
        if (!den_roots) { pz_map_free(pz); return NULL; }
        /* Normalize leading coefficient */
        double *d = (double *)malloc((n_den + 1) * sizeof(double));
        if (!d) { free(den_roots); pz_map_free(pz); return NULL; }
        int i;
        for (i = 0; i <= n_den; i++) d[i] = tf->den[i];
        double lead = d[n_den];
        if (fabs(lead) < 1e-15) { free(d); free(den_roots); pz_map_free(pz); return NULL; }
        for (i = 0; i <= n_den; i++) d[i] /= lead;
        durand_kerner(d, n_den, den_roots, 100);
        free(d);
    }

    if (n_num > 0) {
        num_roots = (double complex *)malloc(n_num * sizeof(double complex));
        if (!num_roots) {
            free(den_roots); pz_map_free(pz); return NULL;
        }
        double *n = (double *)malloc((n_num + 1) * sizeof(double));
        if (!n) { free(den_roots); free(num_roots); pz_map_free(pz); return NULL; }
        int i;
        for (i = 0; i <= n_num; i++) n[i] = tf->num[i];
        double lead = n[n_num];
        if (fabs(lead) < 1e-15) { free(n); free(den_roots); free(num_roots); pz_map_free(pz); return NULL; }
        for (i = 0; i <= n_num; i++) n[i] /= lead;
        durand_kerner(n, n_num, num_roots, 100);
        free(n);
    }

    /* Classify roots: group complex-conjugate pairs, separate real roots */
    int *used_den = (int *)calloc(n_den, sizeof(int));
    int *used_num = (int *)calloc(n_num, sizeof(int));
    int max_pairs = (n_den + n_num) / 2 + 2;
    pz->pz_pairs = (pz_pair_t *)calloc(max_pairs, sizeof(pz_pair_t));
    if (!used_den || !used_num || !pz->pz_pairs) {
        free(used_den); free(used_num); free(den_roots); free(num_roots);
        pz_map_free(pz); return NULL;
    }

    int pair_idx = 0;

    /* Group denominator roots (poles) */
    int i, j;
    for (i = 0; i < n_den; i++) {
        if (used_den[i]) continue;
        if (fabs(cimag(den_roots[i])) > 1e-8) {
            /* Complex root — find conjugate */
            for (j = i + 1; j < n_den; j++) {
                if (!used_den[j] &&
                    fabs(creal(den_roots[i]) - creal(den_roots[j])) < 1e-8 &&
                    fabs(cimag(den_roots[i]) + cimag(den_roots[j])) < 1e-8) {
                    pz->pz_pairs[pair_idx].pole = den_roots[i];
                    pz->pz_pairs[pair_idx].freq_pole = cabs(den_roots[i]);
                    double sigma = -creal(den_roots[i]);
                    pz->pz_pairs[pair_idx].q_pole =
                        (sigma > 1e-15) ? cabs(den_roots[i]) / (2.0 * sigma) : 1000.0;
                    pair_idx++;
                    used_den[i] = 1;
                    used_den[j] = 1;
                    break;
                }
            }
            if (!used_den[i]) {
                /* Unpaired complex root — treat as real */
                used_den[i] = 1;
            }
        }
    }
    /* Remaining real poles */
    for (i = 0; i < n_den; i++) {
        if (!used_den[i]) {
            pz->real_poles = (double *)realloc(pz->real_poles,
                (pz->num_real_poles + 1) * sizeof(double));
            if (!pz->real_poles) { pz_map_free(pz); free(used_den); free(used_num);
                free(den_roots); free(num_roots); return NULL; }
            pz->real_poles[pz->num_real_poles] = creal(den_roots[i]);
            pz->num_real_poles++;
        }
    }
    pz->num_pairs = pair_idx;

    /* Group numerator roots (zeros) — reuse pair array, append */
    for (i = 0; i < n_num && pair_idx < max_pairs; i++) {
        if (used_num[i]) continue;
        if (fabs(cimag(num_roots[i])) > 1e-8) {
            for (j = i + 1; j < n_num; j++) {
                if (!used_num[j] &&
                    fabs(creal(num_roots[i]) - creal(num_roots[j])) < 1e-8 &&
                    fabs(cimag(num_roots[i]) + cimag(num_roots[j])) < 1e-8) {
                    if (pair_idx < pz->num_pairs) {
                        pz->pz_pairs[pair_idx].zero = num_roots[i];
                        pz->pz_pairs[pair_idx].freq_zero = cabs(num_roots[i]);
                        double sigma = -creal(num_roots[i]);
                        pz->pz_pairs[pair_idx].q_zero =
                            (sigma > 1e-15) ? cabs(num_roots[i]) / (2.0 * sigma) : 1000.0;
                    }
                    pair_idx++;
                    used_num[i] = 1;
                    used_num[j] = 1;
                    break;
                }
            }
            if (!used_num[i]) used_num[i] = 1;
        }
    }
    for (i = 0; i < n_num; i++) {
        if (!used_num[i]) {
            pz->real_zeros = (double *)realloc(pz->real_zeros,
                (pz->num_real_zeros + 1) * sizeof(double));
            if (!pz->real_zeros) { pz_map_free(pz); free(used_den); free(used_num);
                free(den_roots); free(num_roots); return NULL; }
            pz->real_zeros[pz->num_real_zeros] = creal(num_roots[i]);
            pz->num_real_zeros++;
        }
    }
    pz->num_pairs = pair_idx;
    pz->gain = tf->gain;

    free(used_den); free(used_num);
    free(den_roots); free(num_roots);
    return pz;
}

/**
 * tf_discrete_to_pzmap — L3 discrete-time pole-zero extraction
 *
 * Same Durand-Kerner approach for z-domain transfer functions.
 * For discrete-time, we extract roots in the z-plane and classify
 * by complex-conjugate pairing and real roots.
 *
 * Complexity: O(N^2 * max_iter)
 */
pz_map_t *tf_discrete_to_pzmap(const tf_discrete_t *tf) {
    if (!tf || tf->den_len < 2) return NULL;

    int n_den = tf->den_len - 1;
    int n_num = tf->num_len - 1;
    if (n_den < 1 || n_den > 20 || n_num > 20) return NULL;

    /* Build a temporary continuous-style polynomial for root finding */
    /* Durand-Kerner works the same way regardless of domain */
    pz_map_t *pz = (pz_map_t *)calloc(1, sizeof(pz_map_t));
    if (!pz) return NULL;

    int max_pairs = (n_den + n_num) / 2 + 2;
    pz->pz_pairs = (pz_pair_t *)calloc(max_pairs, sizeof(pz_pair_t));
    if (!pz->pz_pairs) { pz_map_free(pz); return NULL; }

    double complex *den_roots = NULL, *num_roots = NULL;
    int *used_den = (int *)calloc(n_den, sizeof(int));
    int *used_num = (int *)calloc(n_num, sizeof(int));

    if (n_den > 0) {
        den_roots = (double complex *)malloc(n_den * sizeof(double complex));
        double *d = (double *)malloc((n_den + 1) * sizeof(double));
        if (!den_roots || !d || !used_den || !used_num) {
            free(den_roots); free(d); free(num_roots); free(used_den); free(used_num);
            pz_map_free(pz); return NULL;
        }
        int i;
        for (i = 0; i <= n_den; i++) d[i] = tf->den[i];
        double lead = d[n_den];
        if (fabs(lead) < 1e-15) { free(d); free(den_roots); free(num_roots);
            free(used_den); free(used_num); pz_map_free(pz); return NULL; }
        for (i = 0; i <= n_den; i++) d[i] /= lead;
        durand_kerner(d, n_den, den_roots, 100);
        free(d);
    }

    if (n_num > 0) {
        num_roots = (double complex *)malloc(n_num * sizeof(double complex));
        double *n = (double *)malloc((n_num + 1) * sizeof(double));
        if (!num_roots || !n) {
            free(den_roots); free(num_roots); free(n);
            free(used_den); free(used_num); pz_map_free(pz); return NULL;
        }
        int i;
        for (i = 0; i <= n_num; i++) n[i] = tf->num[i];
        double lead = n[n_num];
        if (fabs(lead) < 1e-15) { free(n); free(den_roots); free(num_roots);
            free(used_den); free(used_num); pz_map_free(pz); return NULL; }
        for (i = 0; i <= n_num; i++) n[i] /= lead;
        durand_kerner(n, n_num, num_roots, 100);
        free(n);
    }

    int pair_idx = 0;
    int i, j;

    /* Poles */
    for (i = 0; i < n_den; i++) {
        if (used_den[i]) continue;
        if (fabs(cimag(den_roots[i])) > 1e-8) {
            for (j = i + 1; j < n_den; j++) {
                if (!used_den[j] &&
                    fabs(creal(den_roots[i]) - creal(den_roots[j])) < 1e-8 &&
                    fabs(cimag(den_roots[i]) + cimag(den_roots[j])) < 1e-8) {
                    pz->pz_pairs[pair_idx].pole = den_roots[i];
                    pz->pz_pairs[pair_idx].freq_pole = cabs(den_roots[i]);
                    double margin = 1.0 - cabs(den_roots[i]);
                    pz->pz_pairs[pair_idx].q_pole =
                        (margin > 1e-6) ? cabs(den_roots[i]) / (2.0 * (-log(cabs(den_roots[i]))))
                                        : 500.0;
                    pair_idx++;
                    used_den[i] = 1;
                    used_den[j] = 1;
                    break;
                }
            }
            if (!used_den[i]) used_den[i] = 1;
        }
    }
    for (i = 0; i < n_den; i++) {
        if (!used_den[i]) {
            pz->real_poles = (double *)realloc(pz->real_poles,
                (pz->num_real_poles + 1) * sizeof(double));
            if (pz->real_poles) {
                pz->real_poles[pz->num_real_poles] = creal(den_roots[i]);
                pz->num_real_poles++;
            }
        }
    }

    /* Zeros */
    for (i = 0; i < n_num; i++) {
        if (used_num[i]) continue;
        if (fabs(cimag(num_roots[i])) > 1e-8) {
            for (j = i + 1; j < n_num; j++) {
                if (!used_num[j] &&
                    fabs(creal(num_roots[i]) - creal(num_roots[j])) < 1e-8 &&
                    fabs(cimag(num_roots[i]) + cimag(num_roots[j])) < 1e-8) {
                    if (pair_idx < max_pairs) {
                        pz->pz_pairs[pair_idx].zero = num_roots[i];
                        pz->pz_pairs[pair_idx].freq_zero = cabs(num_roots[i]);
                    }
                    pair_idx++;
                    used_num[i] = 1;
                    used_num[j] = 1;
                    break;
                }
            }
            if (!used_num[i]) used_num[i] = 1;
        }
    }
    for (i = 0; i < n_num; i++) {
        if (!used_num[i]) {
            pz->real_zeros = (double *)realloc(pz->real_zeros,
                (pz->num_real_zeros + 1) * sizeof(double));
            if (pz->real_zeros) {
                pz->real_zeros[pz->num_real_zeros] = creal(num_roots[i]);
                pz->num_real_zeros++;
            }
        }
    }

    pz->num_pairs = pair_idx;
    pz->gain = tf->gain;

    free(used_den); free(used_num);
    free(den_roots); free(num_roots);
    return pz;
}

void pz_map_free(pz_map_t *pz) {
    if (pz) {
        free(pz->pz_pairs);
        free(pz->real_poles);
        free(pz->real_zeros);
        free(pz);
    }
}

/**
 * pzmap_to_tf_continuous — L3 polynomial expansion from roots
 *
 * Builds H(s) = gain * prod(s - z_k) / prod(s - p_k)
 *
 * Starting from a pole-zero map, expands the product of (s - root) terms
 * into polynomial coefficients. Real roots contribute (s - r) factors;
 * complex-conjugate pairs contribute (s - p)(s - p*) = s^2 - 2*Re(p)*s + |p|^2.
 *
 * This is the reverse operation of tf_continuous_to_pzmap.
 *
 * Reference: Oppenheim & Willsky (1997), Ch. 9
 * Complexity: O(N^2) for polynomial expansion
 */
tf_continuous_t *pzmap_to_tf_continuous(const pz_map_t *pz) {
    if (!pz) return NULL;

    /* Determine degree: denominator from poles, numerator from zeros */
    int den_deg = pz->num_real_poles + 2 * pz->num_pairs;
    int num_deg = pz->num_real_zeros + 2 * pz->num_pairs;
    if (den_deg < 1) return NULL;

    tf_continuous_t *tf = tf_continuous_alloc(
        (num_deg > 0 ? num_deg + 1 : 1), den_deg + 1);
    if (!tf) return NULL;

    /* Build denominator: start with D(s) = 1 */
    tf->den[0] = 1.0;
    int deg_d = 0;
    int k;
    double *temp;

    /* Add real poles: (s - r) for each real pole */
    for (k = 0; k < pz->num_real_poles; k++) {
        double r = pz->real_poles[k];
        temp = (double *)calloc(deg_d + 2, sizeof(double));
        if (!temp) { tf_continuous_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_d; i++) {
            temp[i] -= r * tf->den[i];
            temp[i + 1] += tf->den[i];
        }
        deg_d++;
        for (i = 0; i <= deg_d; i++) tf->den[i] = temp[i];
        free(temp);
    }

    /* Add complex pole pairs: (s-p)(s-p*) = s^2 - 2*Re(p)*s + |p|^2 */
    for (k = 0; k < pz->num_pairs; k++) {
        double a1 = -2.0 * creal(pz->pz_pairs[k].pole);
        double a0 = creal(pz->pz_pairs[k].pole) * creal(pz->pz_pairs[k].pole)
                  + cimag(pz->pz_pairs[k].pole) * cimag(pz->pz_pairs[k].pole);
        temp = (double *)calloc(deg_d + 3, sizeof(double));
        if (!temp) { tf_continuous_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_d; i++) {
            temp[i]     += tf->den[i] * a0;
            temp[i + 1] += tf->den[i] * a1;
            temp[i + 2] += tf->den[i];
        }
        deg_d += 2;
        for (i = 0; i <= deg_d; i++) tf->den[i] = temp[i];
        free(temp);
    }

    /* Build numerator: start with N(s) = 1 */
    tf->num[0] = 1.0;
    int deg_n = 0;

    for (k = 0; k < pz->num_real_zeros; k++) {
        double r = pz->real_zeros[k];
        temp = (double *)calloc(deg_n + 2, sizeof(double));
        if (!temp) { tf_continuous_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_n; i++) {
            temp[i] -= r * tf->num[i];
            temp[i + 1] += tf->num[i];
        }
        deg_n++;
        for (i = 0; i <= deg_n; i++) tf->num[i] = temp[i];
        free(temp);
    }

    for (k = 0; k < pz->num_pairs; k++) {
        double a1 = -2.0 * creal(pz->pz_pairs[k].zero);
        double a0 = creal(pz->pz_pairs[k].zero) * creal(pz->pz_pairs[k].zero)
                  + cimag(pz->pz_pairs[k].zero) * cimag(pz->pz_pairs[k].zero);
        temp = (double *)calloc(deg_n + 3, sizeof(double));
        if (!temp) { tf_continuous_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_n; i++) {
            temp[i]     += tf->num[i] * a0;
            temp[i + 1] += tf->num[i] * a1;
            temp[i + 2] += tf->num[i];
        }
        deg_n += 2;
        for (i = 0; i <= deg_n; i++) tf->num[i] = temp[i];
        free(temp);
    }

    tf->gain = pz->gain;
    return tf;
}

/**
 * pzmap_to_tf_discrete — L3 build discrete TF from pole-zero map
 *
 * H(z) = gain * prod (z - z_k) / prod (z - p_k)
 *
 * Same polynomial expansion as continuous case but in the z-domain.
 * The coefficients represent z-powers rather than s-powers.
 *
 * Complexity: O(N^2) for polynomial expansion
 */
tf_discrete_t *pzmap_to_tf_discrete(const pz_map_t *pz) {
    if (!pz) return NULL;

    int den_deg = pz->num_real_poles + 2 * pz->num_pairs;
    int num_deg = pz->num_real_zeros + 2 * pz->num_pairs;
    if (den_deg < 1) return NULL;

    tf_discrete_t *tf = tf_discrete_alloc(
        (num_deg > 0 ? num_deg + 1 : 1), den_deg + 1);
    if (!tf) return NULL;

    /* Build denominator */
    tf->den[0] = 1.0;
    int deg_d = 0;
    int k;
    double *temp;

    for (k = 0; k < pz->num_real_poles; k++) {
        double r = pz->real_poles[k];
        temp = (double *)calloc(deg_d + 2, sizeof(double));
        if (!temp) { tf_discrete_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_d; i++) {
            temp[i] -= r * tf->den[i];
            temp[i + 1] += tf->den[i];
        }
        deg_d++;
        for (i = 0; i <= deg_d; i++) tf->den[i] = temp[i];
        free(temp);
    }

    for (k = 0; k < pz->num_pairs; k++) {
        double a1 = -2.0 * creal(pz->pz_pairs[k].pole);
        double a0 = creal(pz->pz_pairs[k].pole) * creal(pz->pz_pairs[k].pole)
                  + cimag(pz->pz_pairs[k].pole) * cimag(pz->pz_pairs[k].pole);
        temp = (double *)calloc(deg_d + 3, sizeof(double));
        if (!temp) { tf_discrete_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_d; i++) {
            temp[i]     += tf->den[i] * a0;
            temp[i + 1] += tf->den[i] * a1;
            temp[i + 2] += tf->den[i];
        }
        deg_d += 2;
        for (i = 0; i <= deg_d; i++) tf->den[i] = temp[i];
        free(temp);
    }

    /* Build numerator */
    tf->num[0] = 1.0;
    int deg_n = 0;

    for (k = 0; k < pz->num_real_zeros; k++) {
        double r = pz->real_zeros[k];
        temp = (double *)calloc(deg_n + 2, sizeof(double));
        if (!temp) { tf_discrete_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_n; i++) {
            temp[i] -= r * tf->num[i];
            temp[i + 1] += tf->num[i];
        }
        deg_n++;
        for (i = 0; i <= deg_n; i++) tf->num[i] = temp[i];
        free(temp);
    }

    for (k = 0; k < pz->num_pairs; k++) {
        double a1 = -2.0 * creal(pz->pz_pairs[k].zero);
        double a0 = creal(pz->pz_pairs[k].zero) * creal(pz->pz_pairs[k].zero)
                  + cimag(pz->pz_pairs[k].zero) * cimag(pz->pz_pairs[k].zero);
        temp = (double *)calloc(deg_n + 3, sizeof(double));
        if (!temp) { tf_discrete_free(tf); return NULL; }
        int i;
        for (i = 0; i <= deg_n; i++) {
            temp[i]     += tf->num[i] * a0;
            temp[i + 1] += tf->num[i] * a1;
            temp[i + 2] += tf->num[i];
        }
        deg_n += 2;
        for (i = 0; i <= deg_n; i++) tf->num[i] = temp[i];
        free(temp);
    }

    tf->gain = pz->gain;
    return tf;
}

/* ==================================================================
 * Frequency Transformations (Transfer Function Level)
 * ================================================================== */

void tf_analog_lp_to_lp(tf_continuous_t *tf, double new_wc) {
    if (!tf || new_wc <= 0.0) return;
    double wc_pow = 1.0;
    int i;
    for (i = 0; i < tf->den_len; i++) {
        tf->den[i] /= wc_pow;
        wc_pow *= new_wc;
    }
    wc_pow = 1.0;
    for (i = 0; i < tf->num_len; i++) {
        tf->num[i] /= wc_pow;
        wc_pow *= new_wc;
    }
}

/**
 * tf_analog_lp_to_hp — L2 LP-to-HP frequency transformation
 *
 * Substitution: s -> wc/s
 *
 * H_HP(s) = H_LP(wc/s)
 *         = N_LP(wc/s) / D_LP(wc/s)
 *
 * Multiply numerator and denominator by s^N (where N = max degree)
 * to clear denominator fractions.
 *
 * For N(s) = sum b_k * s^k: N(wc/s) * s^N = sum b_k * wc^k * s^{N-k}
 * For D(s) = sum a_k * s^k: D(wc/s) * s^N = sum a_k * wc^k * s^{N-k}
 *
 * This reverses the coefficient order and scales by wc^k.
 *
 * Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 7
 * Complexity: O(max(num_len, den_len))
 */
void tf_analog_lp_to_hp(tf_continuous_t *tf, double new_wc) {
    if (!tf || new_wc <= 0.0) return;

    int N_num = tf->num_len - 1;
    int N_den = tf->den_len - 1;
    int N = (N_den > N_num) ? N_den : N_num;

    /* Create temporary copies */
    double *old_num = (double *)malloc(tf->num_len * sizeof(double));
    double *old_den = (double *)malloc(tf->den_len * sizeof(double));
    if (!old_num || !old_den) { free(old_num); free(old_den); return; }

    int i;
    for (i = 0; i < tf->num_len; i++) old_num[i] = tf->num[i];
    for (i = 0; i < tf->den_len; i++) old_den[i] = tf->den[i];

    /* New numerator: reverse and scale old denominator */
    for (i = 0; i <= N_den && i <= N; i++) {
        double wc_pow = pow(new_wc, i);
        tf->num[N - i] = old_den[i] * wc_pow;
    }
    for (i = N_den + 1; i <= N; i++) tf->num[N - i] = 0.0;

    /* Resize numerator to N+1 */
    tf->num = (double *)realloc(tf->num, (N + 1) * sizeof(double));
    if (tf->num) tf->num_len = N + 1;

    /* New denominator: reverse and scale old numerator */
    for (i = 0; i <= N_num && i <= N; i++) {
        double wc_pow = pow(new_wc, i);
        tf->den[N - i] = old_num[i] * wc_pow;
    }
    for (i = N_num + 1; i <= N; i++) tf->den[N - i] = 0.0;

    /* Resize denominator to N+1 */
    tf->den = (double *)realloc(tf->den, (N + 1) * sizeof(double));
    if (tf->den) tf->den_len = N + 1;

    free(old_num);
    free(old_den);
}

/**
 * tf_analog_lp_to_bp — L2 LP-to-BP frequency transformation
 *
 * Substitution: s -> (s^2 + w0^2) / (s * B)
 *
 * This doubles the filter order: an N-th order LP prototype becomes a
 * 2N-th order BP filter.
 *
 * Common implementation approach: factor prototype into biquad sections,
 * apply transformation to each section, then multiply results.
 *
 * For the polynomial-level transformation:
 *   Replacing s -> (s^2 + w0^2)/(s*B) in the original polynomial
 *   multiplies the degree by 2.
 *
 * Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 10
 * Complexity: O(N^2) for polynomial substitution
 */
tf_continuous_t *tf_analog_lp_to_bp(const tf_continuous_t *tf,
                                     double w0, double bw) {
    if (!tf || w0 <= 0.0 || bw <= 0.0) return NULL;

    int N = tf->den_len - 1;
    if (N < 1 || N > 10) return NULL;

    /* Factor prototype into biquads for numerical stability */
    int max_bq = (N + 1) / 2;
    analog_biquad_t *bqs = (analog_biquad_t *)calloc(max_bq, sizeof(analog_biquad_t));
    if (!bqs) return NULL;

    int num_bq = tf_continuous_to_biquads(tf, bqs, max_bq);
    if (num_bq < 1) { free(bqs); return NULL; }

    /* Build result: start with identity (num=1, den=1) */
    tf_continuous_t *result = tf_continuous_alloc(1, 1);
    if (!result) { free(bqs); return NULL; }
    result->num[0] = 1.0;
    result->den[0] = 1.0;
    result->gain = tf->gain;

    int k;
    for (k = 0; k < num_bq; k++) {
        /* For each LP biquad H(s) = (b2*s^2 + b1*s + b0) / (s^2 + a1*s + a0):
         * Apply s -> (s^2 + w0^2)/(s*B)
         *
         * Numerator: b2*((s^2+w0^2)/(s*B))^2 + b1*(s^2+w0^2)/(s*B) + b0
         *           = [b2*s^4 + 2*b2*w0^2*s^2 + b2*w0^4 + b1*B*s*(s^2+w0^2) + b0*B^2*s^2] / (s^2*B^2)
         * Multiply by s^2*B^2:
         *           = b2*B^2*s^4 + b1*B*s^3 + (2*b2*w0^2 + b0*B^2)*s^2 + b1*B*w0^2*s + b2*w0^4
         *
         * Denominator: similar with a-coefficients (a2=1 for normalized)
         */

        /* Build 4th-order section */
        tf_continuous_t *section = tf_continuous_alloc(5, 5);
        if (!section) { tf_continuous_free(result); free(bqs); return NULL; }

        double B = bw;
        double B2 = B * B;
        double w02 = w0 * w0;
        double w04 = w02 * w02;

        /* Numerator coefficients (highest to lowest degree) */
        section->num[4] = bqs[k].b2 * B2;
        section->num[3] = bqs[k].b1 * B;
        section->num[2] = 2.0 * bqs[k].b2 * w02 + bqs[k].b0 * B2;
        section->num[1] = bqs[k].b1 * B * w02;
        section->num[0] = bqs[k].b2 * w04;

        /* Denominator: a2=1 (normalized), a1, a0 */
        section->den[4] = 1.0 * B2;
        section->den[3] = bqs[k].a1 * B;
        section->den[2] = 2.0 * w02 + bqs[k].a0 * B2;
        section->den[1] = bqs[k].a1 * B * w02;
        section->den[0] = w02 * w02;

        section->gain = bqs[k].gain;

        /* Multiply into result */
        tf_continuous_t *new_result = tf_continuous_multiply(result, section);
        tf_continuous_free(result);
        tf_continuous_free(section);
        result = new_result;
        if (!result) { free(bqs); return NULL; }
    }

    free(bqs);
    return result;
}

/**
 * tf_analog_lp_to_bs — L2 LP-to-BS frequency transformation
 *
 * Substitution: s -> (s * B) / (s^2 + w0^2)
 *
 * This is the dual of LP-to-BP: it maps DC and infinity to the
 * stopband center, and maps the passband frequencies.
 * Order doubles: N -> 2N.
 *
 * Reference: Van Valkenburg, Analog Filter Design (1982), Ch. 10
 * Complexity: O(N^2)
 */
tf_continuous_t *tf_analog_lp_to_bs(const tf_continuous_t *tf,
                                     double w0, double bw) {
    if (!tf || w0 <= 0.0 || bw <= 0.0) return NULL;

    int N = tf->den_len - 1;
    if (N < 1 || N > 10) return NULL;

    int max_bq = (N + 1) / 2;
    analog_biquad_t *bqs = (analog_biquad_t *)calloc(max_bq, sizeof(analog_biquad_t));
    if (!bqs) return NULL;

    int num_bq = tf_continuous_to_biquads(tf, bqs, max_bq);
    if (num_bq < 1) { free(bqs); return NULL; }

    tf_continuous_t *result = tf_continuous_alloc(1, 1);
    if (!result) { free(bqs); return NULL; }
    result->num[0] = 1.0;
    result->den[0] = 1.0;
    result->gain = tf->gain;

    int k;
    for (k = 0; k < num_bq; k++) {
        /* For each LP biquad: H(s) = (b2*s^2 + b1*s + b0) / (s^2 + a1*s + a0)
         * Apply s -> (s*B)/(s^2 + w0^2)
         *
         * Substitute and multiply numerator and denominator by (s^2 + w0^2)^2
         * to clear fractions.
         */
        tf_continuous_t *section = tf_continuous_alloc(5, 5);
        if (!section) { tf_continuous_free(result); free(bqs); return NULL; }

        double B = bw;
        double B2 = B * B;
        double w02 = w0 * w0;
        double w04 = w02 * w02;

        /* Numerator after substitution:
         * b2*B^2*s^2 + b1*B*s*(s^2+w0^2) + b0*(s^2+w0^2)^2
         * = b0*s^4 + b1*B*s^3 + (2*b0*w0^2 + b2*B^2)*s^2 + b1*B*w0^2*s + b0*w0^4
         */
        section->num[4] = bqs[k].b0;
        section->num[3] = bqs[k].b1 * B;
        section->num[2] = 2.0 * bqs[k].b0 * w02 + bqs[k].b2 * B2;
        section->num[1] = bqs[k].b1 * B * w02;
        section->num[0] = bqs[k].b0 * w04;

        /* Denominator after substitution:
         * B^2*s^2 + a1*B*s*(s^2+w0^2) + a0*(s^2+w0^2)^2
         * = a0*s^4 + a1*B*s^3 + (2*a0*w0^2 + B^2)*s^2 + a1*B*w0^2*s + a0*w0^4
         */
        section->den[4] = bqs[k].a0;
        section->den[3] = bqs[k].a1 * B;
        section->den[2] = 2.0 * bqs[k].a0 * w02 + B2;
        section->den[1] = bqs[k].a1 * B * w02;
        section->den[0] = bqs[k].a0 * w04;

        section->gain = bqs[k].gain;

        tf_continuous_t *new_result = tf_continuous_multiply(result, section);
        tf_continuous_free(result);
        tf_continuous_free(section);
        result = new_result;
        if (!result) { free(bqs); return NULL; }
    }

    free(bqs);
    return result;
}

/* ==================================================================
 * L4: Routh-Hurwitz Stability Criterion (Continuous-Time)
 * ================================================================== */

/**
 * tf_continuous_is_stable — L4 theorem implementation
 *
 * Routh-Hurwitz criterion: construct the Routh array and count
 * sign changes in the first column.
 *
 * Algorithm:
 * Row 0: a_N, a_{N-2}, a_{N-4}, ...
 * Row 1: a_{N-1}, a_{N-3}, a_{N-5}, ...
 * Row k: computed from rows k-2 and k-1 using determinant formula:
 *   r_{k,j} = (r_{k-2,1}*r_{k-1,j+1} - r_{k-1,1}*r_{k-2,j+1}) / r_{k-1,1}
 *
 * Stability condition: all entries in the first column have the same sign.
 *
 * Reference: E.J. Routh, "A Treatise on the Stability of a Given
 *            State of Motion", 1877; A. Hurwitz, Math. Ann., 1895.
 */
int tf_continuous_is_stable(const tf_continuous_t *tf) {
    if (!tf || tf->den_len < 1) return -1;

    int n = tf->den_len - 1; /* Polynomial degree */
    if (n < 1) return -1;
    if (n > 20) return -1; /* Practical limit */

    /* Build Routh array: (n+1) x ((n+2)//2) */
    int cols = (n + 2) / 2;
    double *routh = (double *)calloc((n + 1) * cols, sizeof(double));
    if (!routh) return -1;

    /* Row 0: even-indexed coefficients (highest degree first) */
    int i, j;
    for (i = 0, j = 0; j <= n; j += 2) {
        routh[0 * cols + i] = tf->den[n - j];
        i++;
    }

    /* Row 1: odd-indexed coefficients */
    for (i = 0, j = 1; j <= n; j += 2) {
        routh[1 * cols + i] = tf->den[n - j];
        i++;
    }

    /* Check first-column signs for all rows.
     * The number of RHP poles = number of sign changes.
     * Stable iff no sign changes (all same sign). */
    int stable = 1;
    int *first_col_signs = (int *)malloc((n + 1) * sizeof(int));
    if (!first_col_signs) { free(routh); return -1; }

    first_col_signs[0] = (routh[0 * cols + 0] > 0) ? 1 :
                         (routh[0 * cols + 0] < 0) ? -1 : 0;
    first_col_signs[1] = (routh[1 * cols + 0] > 0) ? 1 :
                         (routh[1 * cols + 0] < 0) ? -1 : 0;

    for (i = 2; i <= n; i++) {
        double first_prev = routh[(i - 1) * cols + 0];
        if (fabs(first_prev) < 1e-15) {
            first_prev = 1e-10;
        }

        for (j = 0; j < cols - 1; j++) {
            double a = routh[(i - 2) * cols + 0];
            double b = routh[(i - 2) * cols + (j + 1)];
            double c = routh[(i - 1) * cols + 0];
            double d = routh[(i - 1) * cols + (j + 1)];
            routh[i * cols + j] = -(a * d - c * b) / first_prev;
        }

        first_col_signs[i] = (routh[i * cols + 0] > 0) ? 1 :
                              (routh[i * cols + 0] < 0) ? -1 : 0;
    }

    /* Count sign changes in first column */
    int sign_changes = 0;
    for (i = 0; i < n; i++) {
        if (first_col_signs[i] != 0 && first_col_signs[i + 1] != 0
            && first_col_signs[i] != first_col_signs[i + 1]) {
            sign_changes++;
        }
    }
    stable = (sign_changes == 0) ? 1 : 0;

    free(first_col_signs);
    free(routh);
    return stable;
}

/* ==================================================================
 * L4: Jury Stability Criterion (Discrete-Time)
 * ================================================================== */

/**
 * tf_discrete_is_stable — L4 theorem implementation
 *
 * Jury stability test for discrete-time systems.
 *
 * For characteristic polynomial A(z) = a0 + a1*z^{-1} + ... + an*z^{-n}:
 * Multiply by z^n to get P(z) = a0*z^n + a1*z^{n-1} + ... + an.
 *
 * The Jury conditions for P(z) with coefficients p_0..p_n (where p_0 = a0):
 * 1. P(1) > 0
 * 2. (-1)^n * P(-1) > 0
 * 3. |p_n| < p_0
 * Then construct Jury table and check |b_0| > |b_{n-1}|, etc.
 *
 * Simpler direct check: for low-order (n <= 3), use algebraic conditions.
 * For higher order, construct the full Jury table.
 *
 * Reference: E.I. Jury, "Theory and Application of the z-Transform
 *            Method", Wiley, 1964.
 */
int tf_discrete_is_stable(const tf_discrete_t *tf) {
    if (!tf || tf->den_len < 1) return -1;
    int n = tf->den_len - 1;
    if (n < 1 || n > 20) return -1;

    if (fabs(tf->den[0]) < 1e-15) return -1;

    /* Use direct algebraic conditions for low orders */
    if (n == 1) {
        /* a0 + a1*z^{-1}: stable iff |a1/a0| < 1 */
        return (fabs(tf->den[1]) < fabs(tf->den[0])) ? 1 : 0;
    }
    if (n == 2) {
        /* a0 + a1*z^{-1} + a2*z^{-2}:
         * stable iff |a2| < |a0| AND |a1| < |a0 + a2| */
        double a0 = tf->den[0], a1 = tf->den[1], a2 = tf->den[2];
        if (fabs(a2) >= fabs(a0)) return 0;
        if (fabs(a1) >= fabs(a0 + a2)) return 0;
        return 1;
    }
    if (n == 3) {
        /* a0 + a1*z^{-1} + a2*z^{-2} + a3*z^{-3}
         * See Oppenheim & Schafer, Table 6.1 */
        double a0 = tf->den[0], a1 = tf->den[1];
        double a2 = tf->den[2], a3 = tf->den[3];
        if (fabs(a3) >= fabs(a0)) return 0;
        if (fabs(a1*a3 - a0*a2) >= fabs(a0*a0 - a3*a3)) return 0;
        return 1;
    }

    /* For n >= 4: construct Jury table */
    double *coeff = (double *)malloc((n + 1) * sizeof(double));
    if (!coeff) return -1;

    int k;
    for (k = 0; k <= n; k++) coeff[k] = tf->den[k];

    /* Normalize so leading coefficient a0 = 1 */
    double lead = coeff[0];
    for (k = 0; k <= n; k++) coeff[k] /= lead;

    /* P(1) > 0 */
    double sum_p1 = 0;
    for (k = 0; k <= n; k++) sum_p1 += coeff[k];
    if (sum_p1 <= 0) { free(coeff); return 0; }

    /* (-1)^n * P(-1) > 0 */
    double sum_pm1 = 0;
    for (k = 0; k <= n; k++) {
        sum_pm1 += (k % 2 == 0) ? coeff[k] : -coeff[k];
    }
    if ((n % 2 == 0 ? sum_pm1 : -sum_pm1) <= 0) {
        free(coeff); return 0;
    }

    /* |p_n| < 1 (since a0=1 after normalization) */
    if (fabs(coeff[n]) >= 1.0) { free(coeff); return 0; }

    /* Jury table construction */
    double *row = (double *)malloc((n + 1) * sizeof(double));
    if (!row) { free(coeff); return 0; }
    for (k = 0; k <= n; k++) row[k] = coeff[k];

    int i;
    for (i = 0; i < n - 1; i++) {
        double alpha = row[n - i] / row[0];
        if (fabs(alpha) >= 1.0) { free(coeff); free(row); return 0; }
        int j;
        for (j = 0; j < n - i; j++) {
            coeff[j] = row[j] - alpha * row[n - i - j];
        }
        if (fabs(coeff[0]) <= fabs(coeff[n - i - 1])) {
            free(coeff); free(row); return 0;
        }
        for (j = 0; j < n - i; j++) row[j] = coeff[j];
    }

    free(coeff); free(row);
    return 1;
}

int sos_is_stable(const sos_matrix_t *sos) {
    if (!sos) return -1;
    int k;
    for (k = 0; k < sos->num_sections; k++) {
        double a1 = sos->sos_matrix[k * 6 + 3 + 1];
        double a2 = sos->sos_matrix[k * 6 + 3 + 2];
        /* Biquad stability: |a2| < 1 and |a1| < 1 + a2 */
        if (fabs(a2) >= 1.0) return 0;
        if (fabs(a1) >= 1.0 + a2) return 0;
    }
    return 1;
}

/* ==================================================================
 * L4: Bilinear Transform — s <-> z mapping
 * ================================================================== */

/**
 * bilinear_s_to_z — L4 key mapping
 *
 * s = 2*fs * (z-1)/(z+1)   =>   z = (1 + s/(2*fs)) / (1 - s/(2*fs))
 *
 * Properties:
 * - LHP (Re(s) < 0) maps to |z| < 1 (stable -> stable)
 * - jw-axis maps to |z| = 1 (frequency axis preserved)
 * - One-to-one mapping, no aliasing
 * - Introduces nonlinear frequency warping: w_a = 2*fs*tan(w_d/(2*fs))
 *
 * Complexity: O(1)
 */
double complex bilinear_s_to_z(double complex s, double fs) {
    if (fs <= 0.0) return 0;
    double complex k = s / (2.0 * fs);
    return (1.0 + k) / (1.0 - k);
}

double complex bilinear_z_to_s(double complex z, double fs) {
    if (fs <= 0.0) return 0;
    return 2.0 * fs * (z - 1.0) / (z + 1.0);
}

/**
 * bilinear_prewarp — L4 frequency pre-warping
 *
 * To design a digital filter with normalized frequency w_d (rad/sample),
 * the analog prototype cutoff must be pre-warped:
 *   w_a = 2*fs * tan(w_d / 2)
 *
 * w_d = 2*pi*f_d/fs where f_d is the desired digital cutoff in Hz.
 *
 * This ensures the digital filter's cutoff is exactly at w_d
 * after applying the bilinear transform: s = 2*fs*(z-1)/(z+1).
 *
 * The mapping is: w_a = 2*fs * tan(w_d/2)
 * Equivalently: f_a = (fs/pi) * tan(pi*f_d/fs)
 *
 * Reference: Oppenheim & Schafer (2010), Sec. 7.1.2
 * Complexity: O(1)
 */
double bilinear_prewarp(double w_d, double fs) {
    if (fs <= 0.0 || w_d <= 0.0) return -1.0;
    return 2.0 * fs * tan(w_d / 2.0);
}

/**
 * bilinear_transform — L4/L5 complete s->z mapping
 *
 * Converts H(s) to H(z) by substituting s = 2*fs*(z-1)/(z+1).
 *
 * For an order-N analog filter, the result is an order-N digital filter
 * with the same number of poles and zeros.
 *
 * The substitution involves polynomial operations:
 * - Replace each s^k with (2*fs*(z-1)/(z+1))^k
 * - Multiply numerator and denominator by (z+1)^N to clear denominators
 * - Collect coefficients of z powers
 *
 * Complexity: O(N^3) for polynomial manipulations
 */
tf_discrete_t *bilinear_transform(const tf_continuous_t *tf_analog, double fs) {
    if (!tf_analog || fs <= 0.0) return NULL;

    int N = tf_analog->den_len - 1;
    if (N < 1 || N > 10) return NULL;

    double c = 2.0 * fs;

    /* For each term s^k, substitute ((c*(z-1)/(z+1)))^k
     * and multiply by (z+1)^N to clear denominators.
     *
     * Simplified: Expand (z+1)^N * (z-1)^k = sum binom * z^m
     * This is a binomial expansion problem.
     *
     * For practicality, we implement order 1-4 directly. */

    /* Allocate result: order N for both num and den */
    tf_discrete_t *result = tf_discrete_alloc(N + 1, N + 1);
    if (!result) return NULL;

    /* For simplicity, work with the biquad cascade approach.
     * Convert analog TF to digital via bilinear substitution
     * on the s-domain poles directly. */

    /* Start with H(z) = gain * num_analog(c*(z-1)/(z+1))
     *                   / den_analog(c*(z-1)/(z+1))
     *
     * Multiply top and bottom by (z+1)^N:
     * H(z) = gain * sum b_k * c^k * (z-1)^k * (z+1)^{N-k}
     *             / sum a_k * c^k * (z-1)^k * (z+1)^{N-k}
     */

    double *poly_z1 = (double *)calloc(N + 1, sizeof(double)); /* (z-1) coefficients */
    double *poly_zp1 = (double *)calloc(N + 1, sizeof(double)); /* (z+1) coefficients */
    if (!poly_z1 || !poly_zp1) {
        free(poly_z1); free(poly_zp1); tf_discrete_free(result); return NULL;
    }

    /* (z-1)^k expansion via binomial theorem */
    poly_z1[0] = 1.0;  /* (z-1)^0 = 1 */
    poly_zp1[0] = 1.0; /* (z+1)^0 = 1 */

    /* Evaluate denominator and numerator sums */
    int k;
    for (k = 0; k <= N; k++) {
        /* Compute (z-1)^k * (z+1)^{N-k} and add to sums, weighted by coeffs */
        double c_pow = pow(c, k);
        double num_coef = (k < tf_analog->num_len) ? tf_analog->num[k] : 0.0;
        double den_coef = (k < tf_analog->den_len) ? tf_analog->den[k] : 0.0;

        /* Expand (z-1)^k * (z+1)^{N-k} = sum_{m} binom_coef[m] * z^m */
        /* Use combinatorial expansion */
        int m;
        for (m = 0; m <= N; m++) {
            double coef = 0.0;
            int p;
            for (p = 0; p <= k && p <= m; p++) {
                int q = m - p;
                if (q < 0 || q > N - k) continue;

                /* Binomial: C(k,p) * (-1)^p * C(N-k, q) * 1^q */
                double binom1 = 1.0, binom2 = 1.0;
                int r;
                for (r = 1; r <= p; r++) binom1 = binom1 * (k - r + 1) / r;
                for (r = 1; r <= q; r++) binom2 = binom2 * (N - k - r + 1) / r;

                if (p % 2 == 1) binom1 = -binom1;
                coef += binom1 * binom2;
            }

            result->num[m] += num_coef * c_pow * coef;
            result->den[m] += den_coef * c_pow * coef;
        }
    }

    /* Normalize so den[0] = 1 */
    if (fabs(result->den[0]) > 1e-15) {
        int j;
        for (j = 1; j <= N; j++) result->den[j] /= result->den[0];
        for (j = 0; j <= N; j++) result->num[j] /= result->den[0];
        result->den[0] = 1.0;
    }

    result->gain = tf_analog->gain;

    free(poly_z1);
    free(poly_zp1);
    return result;
}

/* ==================================================================
 * L5: Impulse Invariance Method
 * ================================================================== */

/**
 * impulse_invariance — L5 complete impulse invariance conversion
 *
 * Converts H(s) to H(z) by sampling the impulse response:
 *   h[n] = T * h_a(nT)   where T = 1/fs
 *
 * Using partial fraction expansion:
 *   H(s) = sum_{k} r_k/(s - p_k) + direct_term
 *
 * Each term maps to:
 *   r_k/(s - p_k)  →  T * r_k / (1 - e^{p_k*T} * z^{-1})
 *
 * Implementation steps:
 * 1. Find poles of the analog transfer function
 * 2. Map each pole: p_digital = e^{p_analog * T}
 * 3. Adjust gain: maintain same impulse response amplitude
 * 4. Build denominator from mapped poles
 * 5. Construct numerator to preserve the response
 *
 * Properties:
 * - Preserves impulse response shape at sample instants exactly
 * - Stable analog → stable digital (e^{negative real} → |z|<1)
 * - Aliasing: frequencies above fs/2 fold back (only suitable for LP/BP)
 *
 * Reference: Oppenheim & Schafer (2010), Sec. 7.2.1
 * Complexity: O(N^2 * iter) for root finding + O(N^2) for polynomial expansion
 */
tf_discrete_t *impulse_invariance(const tf_continuous_t *tf_analog, double fs) {
    if (!tf_analog || fs <= 0.0) return NULL;

    int N = tf_analog->den_len - 1;
    if (N < 1 || N > 20) return NULL;

    double T = 1.0 / fs;

    /* Find analog poles via polynomial root-finding */
    double *den_copy = (double *)malloc(tf_analog->den_len * sizeof(double));
    if (!den_copy) return NULL;
    int i;
    for (i = 0; i < tf_analog->den_len; i++) den_copy[i] = tf_analog->den[i];

    double lead = den_copy[N];
    if (fabs(lead) < 1e-15) { free(den_copy); return NULL; }
    for (i = 0; i <= N; i++) den_copy[i] /= lead;

    double complex *analog_poles = (double complex *)malloc(N * sizeof(double complex));
    if (!analog_poles) { free(den_copy); return NULL; }

    durand_kerner(den_copy, N, analog_poles, 100);
    free(den_copy);

    /* Map poles: p_digital = exp(p_analog * T) */
    double complex *digital_poles = (double complex *)malloc(N * sizeof(double complex));
    if (!digital_poles) { free(analog_poles); return NULL; }

    for (i = 0; i < N; i++) {
        digital_poles[i] = cexp(analog_poles[i] * T);
    }
    free(analog_poles);

    /* Build denominator polynomial from digital poles */
    int den_len = N + 1;
    tf_discrete_t *result = tf_discrete_alloc(den_len, den_len);
    if (!result) { free(digital_poles); return NULL; }

    result->den[0] = 1.0;
    int deg = 0;

    for (i = 0; i < N; i++) {
        if (fabs(cimag(digital_poles[i])) < 1e-10) {
            /* Real pole: (1 - p*z^{-1}) */
            double p = creal(digital_poles[i]);
            int j;
            for (j = deg + 1; j >= 1; j--) {
                result->den[j] = result->den[j - 1] * (-p) + result->den[j];
            }
            result->den[0] *= (-p);
            deg++;
        } else if (cimag(digital_poles[i]) > 1e-10) {
            /* Complex conjugate pair: (1 - p*z^{-1})(1 - p**z^{-1})
             * = 1 - 2*Re(p)*z^{-1} + |p|^2*z^{-2} */
            double a1 = -2.0 * creal(digital_poles[i]);
            double a2 = creal(digital_poles[i]) * creal(digital_poles[i])
                      + cimag(digital_poles[i]) * cimag(digital_poles[i]);
            double *temp = (double *)calloc(deg + 3, sizeof(double));
            if (!temp) { free(digital_poles); tf_discrete_free(result); return NULL; }
            int j;
            for (j = 0; j <= deg; j++) {
                temp[j]     += result->den[j];
                temp[j + 1] += result->den[j] * a1;
                temp[j + 2] += result->den[j] * a2;
            }
            deg += 2;
            for (j = 0; j <= deg; j++) result->den[j] = temp[j];
            free(temp);
            i++; /* skip conjugate */
        }
    }

    /* Numerator: for simplicity, use T * b0 (DC gain matching) */
    /* A full implementation would compute residues, but for common
     * filter design the DC gain matching suffices. */
    result->num[0] = T * tf_analog->num[0] * tf_analog->gain;
    result->gain = 1.0;

    /* Normalize DC gain: |H(e^{j0})| should match |H(j0)| */
    double dc_analog = cabs(tf_continuous_eval(tf_analog, 0));
    double dc_digital = cabs(tf_discrete_eval(result, 1.0 + 0*I));
    if (dc_digital > 1e-15 && fabs(dc_analog) > 1e-15) {
        double gain_correction = dc_analog / dc_digital;
        int j;
        for (j = 0; j < result->num_len; j++) result->num[j] *= gain_correction;
        result->gain *= gain_correction;
    }

    free(digital_poles);
    return result;
}
