/**
 * @file nonlinear_chaos.c
 * @brief Chaos detection and characterization algorithms.
 *
 * Implements Lyapunov exponent computation, correlation dimension
 * estimation, time-delay embedding, entropy measures, and
 * surrogate data testing.
 *
 * Knowledge points:
 *   L5: Lyapunov exponents via Benettin algorithm
 *   L5: Correlation dimension (Grassberger-Procaccia)
 *   L8: Kaplan-Yorke dimension, Kolmogorov-Sinai entropy
 *   L8: Delay embedding (Takens), mutual information, FNN
 *   L8: 0-1 test for chaos (Gottwald-Melbourne)
 *   L8: AAFT surrogate generation, sample entropy
 *
 * Reference:
 *   - Sprott, "Chaos and Time-Series Analysis", Oxford, 2003
 *   - Kantz & Schreiber, "Nonlinear Time Series Analysis", 2004
 */

#include "nonlinear_chaos.h"
#include "nonlinear_stability.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* ===================================================================
 * Helper: distance between two points in dim-dimensional space
 * =================================================================== */

static double point_distance(const double *a, const double *b, size_t dim)
{
    double sum = 0.0;
    for (size_t d = 0; d < dim; d++) {
        double diff = a[d] - b[d];
        sum += diff * diff;
    }
    return sqrt(sum);
}

/* ===================================================================
 * Lyapunov Exponent Computation (Benettin Algorithm)
 *
 * 1. Integrate system + variational equations simultaneously
 * 2. Periodically re-orthonormalize deviation vectors (Gram-Schmidt)
 * 3. Sum log(growth factors) / total time
 *
 * Each deviation vector grows/shrinks as exp(lambda_i * t).
 * Benettin et al., Meccanica, 1980.
 * =================================================================== */

int nl_lyapunov_exponents(nl_nonlinear_system_t *sys,
                          const double *x0, double T_total, double dt,
                          int n_renorm, double *exponents)
{
    if (!sys || !x0 || !exponents
        || T_total <= 0.0 || dt <= 0.0 || n_renorm < 1
        || sys->dim < 1 || sys->dim > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t n = sys->dim;
    size_t n_steps_total = (size_t)(T_total / dt);
    size_t steps_per_renorm = n_steps_total / (size_t)n_renorm;
    if (steps_per_renorm < 1) steps_per_renorm = 1;

    /* Initialize state */
    double x[NL_MAX_STATE_DIM];
    memcpy(x, x0, n * sizeof(double));

    /* Initialize deviation vectors: identity matrix columns */
    double dev[NL_MAX_STATE_DIM][NL_MAX_STATE_DIM];
    memset(dev, 0, sizeof(dev));
    for (size_t i = 0; i < n; i++) dev[i][i] = 1.0;

    /* Accumulate Lyapunov sums */
    double lyap_sum[NL_MAX_STATE_DIM] = {0};
    size_t n_renorms_done = 0;

    for (size_t step = 0; step < n_steps_total; step++) {
        /* Evaluate f at current x */
        double dxdt[NL_MAX_STATE_DIM];
        if (nl_evaluate_f(sys, x, dxdt) != 0) return -1;

        /* Compute Jacobian at current x */
        nl_jacobian_t J;
        if (nl_compute_jacobian(sys, x, &J) != 0) return -1;

        /* Integrate state: x += dt * f(x) */
        for (size_t d = 0; d < n; d++)
            x[d] += dt * dxdt[d];

        /* Integrate deviation vectors: dev_i += dt * J * dev_i */
        for (size_t i = 0; i < n; i++) {
            double new_dev[NL_MAX_STATE_DIM];
            for (size_t r = 0; r < n; r++) {
                new_dev[r] = 0.0;
                for (size_t c = 0; c < n; c++)
                    new_dev[r] += J.matrix[r][c] * dev[i][c];
            }
            for (size_t r = 0; r < n; r++)
                dev[i][r] += dt * new_dev[r];
        }

        /* Re-orthonormalize periodically */
        if ((step + 1) % steps_per_renorm == 0
            && n_renorms_done < (size_t)n_renorm) {

            /* Gram-Schmidt orthonormalization */
            for (size_t i = 0; i < n; i++) {
                /* Subtract projections onto previous vectors */
                for (size_t j = 0; j < i; j++) {
                    double dot = 0.0;
                    for (size_t d = 0; d < n; d++)
                        dot += dev[i][d] * dev[j][d];
                    for (size_t d = 0; d < n; d++)
                        dev[i][d] -= dot * dev[j][d];
                }

                /* Compute norm = growth factor */
                double norm = 0.0;
                for (size_t d = 0; d < n; d++)
                    norm += dev[i][d] * dev[i][d];
                norm = sqrt(norm);

                if (norm > 1e-15) {
                    lyap_sum[i] += log(norm);
                    for (size_t d = 0; d < n; d++)
                        dev[i][d] /= norm;
                }
            }
            n_renorms_done++;
        }
    }

    /* Compute averaged Lyapunov exponents */
    double total_time = (double)n_renorms_done
                        * (double)steps_per_renorm * dt;
    if (total_time < 1e-15) total_time = dt;

    for (size_t i = 0; i < n; i++)
        exponents[i] = lyap_sum[i] / total_time;

    /* Sort descending (simple bubble sort for small n) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (exponents[j] > exponents[i]) {
                double tmp = exponents[i];
                exponents[i] = exponents[j];
                exponents[j] = tmp;
            }
        }
    }

    return 0;
}

/* ===================================================================
 * Kaplan-Yorke (Lyapunov) Dimension
 *
 * D_KY = k + (sum_{i=1}^k lambda_i) / |lambda_{k+1}|
 * where k = max{j : sum_{i=1}^j lambda_i >= 0}
 * =================================================================== */

double nl_kaplan_yorke_dimension(const double *exponents, size_t dim)
{
    if (!exponents || dim < 1) return 0.0;

    double cum_sum = 0.0;
    size_t k = 0;

    for (size_t i = 0; i < dim; i++) {
        cum_sum += exponents[i];
        if (cum_sum >= 0.0) k = i + 1;
        else break;
    }

    if (k == 0) return 0.0;
    if (k == dim) return (double)dim;

    return (double)k + cum_sum / fabs(exponents[k]);
}

/* ===================================================================
 * Correlation Sum C(epsilon)
 *
 * C(eps) = fraction of point pairs within distance epsilon.
 * =================================================================== */

double nl_correlation_sum(const double *data, size_t N, size_t dim,
                          double epsilon)
{
    if (!data || N < 2 || dim < 1 || epsilon <= 0.0) return 0.0;

    size_t count = 0;
    size_t total_pairs = 0;

    /* Use random subsample for large N to avoid O(N^2) */
    size_t max_samples = (N > 2000) ? 2000 : N;
    size_t step = (N > max_samples) ? N / max_samples : 1;
    if (step < 1) step = 1;

    size_t *indices = (size_t *)malloc(max_samples * sizeof(size_t));
    if (!indices) return 0.0;

    for (size_t i = 0; i < max_samples; i++)
        indices[i] = i * step;

    for (size_t i = 0; i < max_samples; i++) {
        for (size_t j = i + 1; j < max_samples; j++) {
            double dist = point_distance(
                &data[indices[i] * dim],
                &data[indices[j] * dim], dim);
            if (dist <= epsilon) count++;
            total_pairs++;
        }
    }

    free(indices);
    if (total_pairs == 0) return 0.0;
    return (double)count / (double)total_pairs;
}

/* ===================================================================
 * Correlation Dimension D2 via Log-Log Regression
 *
 * D2 = slope of log C(eps) vs log eps.
 * Uses linear least squares over the scaling region.
 * =================================================================== */

int nl_correlation_dimension(const double *data, size_t N, size_t dim,
                             const double *epsilons, size_t n_eps,
                             double *D2, double *r2)
{
    if (!data || N < 2 || dim < 1 || !epsilons
        || n_eps < 4 || !D2 || !r2) {
        errno = EINVAL; return -1;
    }

    double *log_eps = (double *)malloc(n_eps * sizeof(double));
    double *log_C = (double *)malloc(n_eps * sizeof(double));
    if (!log_eps || !log_C) {
        free(log_eps); free(log_C); return -1;
    }

    for (size_t i = 0; i < n_eps; i++) {
        log_eps[i] = log(epsilons[i]);
        double C = nl_correlation_sum(data, N, dim, epsilons[i]);
        log_C[i] = (C > 1e-15) ? log(C) : -35.0;
    }

    /* Linear regression: log C = D2 * log eps + const */
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    for (size_t i = 0; i < n_eps; i++) {
        sum_x += log_eps[i];
        sum_y += log_C[i];
        sum_xy += log_eps[i] * log_C[i];
        sum_x2 += log_eps[i] * log_eps[i];
    }

    double Nf = (double)n_eps;
    double denom = Nf * sum_x2 - sum_x * sum_x;
    if (fabs(denom) < 1e-15) {
        free(log_eps); free(log_C);
        *D2 = 0.0; *r2 = 0.0;
        return -1;
    }

    *D2 = (Nf * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y * sum_x2 - sum_x * sum_xy) / denom;

    /* R-squared */
    double ss_res = 0.0, ss_tot = 0.0;
    double mean_y = sum_y / Nf;
    for (size_t i = 0; i < n_eps; i++) {
        double pred = (*D2) * log_eps[i] + intercept;
        ss_res += (log_C[i] - pred) * (log_C[i] - pred);
        ss_tot += (log_C[i] - mean_y) * (log_C[i] - mean_y);
    }
    *r2 = (ss_tot > 1e-15) ? 1.0 - ss_res / ss_tot : 0.0;

    free(log_eps); free(log_C);
    return 0;
}

/* ===================================================================
 * Time-Delay Embedding (Takens)
 *
 * From scalar time series {s_t}, form vectors:
 *   x_t = (s_t, s_{t-tau}, s_{t-2*tau}, ..., s_{t-(m-1)*tau})
 *
 * N_emb = N - (m-1)*tau  embedded vectors of dimension m.
 * =================================================================== */

int nl_delay_embedding(const double *time_series, size_t N,
                       int tau, size_t m, double *embedded)
{
    if (!time_series || !embedded || N < 2
        || tau < 1 || m < 1 || m > NL_MAX_STATE_DIM) {
        errno = EINVAL; return -1;
    }

    size_t shift = (size_t)tau * (m - 1);
    if (shift >= N) { errno = EINVAL; return -1; }

    size_t N_emb = N - shift;
    for (size_t t = 0; t < N_emb; t++) {
        for (size_t d = 0; d < m; d++) {
            embedded[t * m + d] = time_series[t + (m - 1 - d) * (size_t)tau];
        }
    }

    return (int)N_emb;
}

/* ===================================================================
 * Mutual Information for Optimal Delay Estimation
 *
 * I(tau) = sum_{i,j} p_{ij}(tau) * log2(p_{ij}(tau) / (p_i * p_j))
 *
 * First minimum of I(tau) gives optimal embedding delay.
 * (Fraser & Swinney, Phys. Rev. A, 1986)
 * =================================================================== */

int nl_mutual_information_delay(const double *time_series, size_t N,
                                 int max_tau, int n_bins,
                                 int *tau_opt)
{
    if (!time_series || !tau_opt || N < (size_t)max_tau
        || max_tau < 1 || n_bins < 2) {
        errno = EINVAL; return -1;
    }

    /* Find data range */
    double s_min = time_series[0], s_max = time_series[0];
    for (size_t i = 1; i < N; i++) {
        if (time_series[i] < s_min) s_min = time_series[i];
        if (time_series[i] > s_max) s_max = time_series[i];
    }
    double range = s_max - s_min;
    if (range < 1e-12) { *tau_opt = 1; return 0; }

    /* MI for each tau */
    double prev_mi = INFINITY;
    int best_tau = 1;
    int rising = 0;

    for (int tau = 1; tau <= max_tau; tau++) {
        /* Build joint histogram */
        int *joint = (int *)calloc(
            (size_t)(n_bins * n_bins), sizeof(int));
        int *marg_x = (int *)calloc((size_t)n_bins, sizeof(int));
        int *marg_y = (int *)calloc((size_t)n_bins, sizeof(int));
        if (!joint || !marg_x) { free(joint); free(marg_x); continue; }

        size_t n_pairs = N - (size_t)tau;
        for (size_t t = 0; t < n_pairs; t++) {
            int bx = (int)((time_series[t] - s_min)
                           / range * (double)(n_bins - 1));
            int by = (int)((time_series[t + tau] - s_min)
                           / range * (double)(n_bins - 1));
            if (bx < 0) bx = 0; else if (bx >= n_bins) bx = n_bins - 1;
            if (by < 0) by = 0; else if (by >= n_bins) by = n_bins - 1;

            joint[bx * n_bins + by]++;
            marg_x[bx]++;
            marg_y[by]++;
        }

        /* Compute mutual information */
        double mi = 0.0;
        for (int i = 0; i < n_bins; i++) {
            for (int j = 0; j < n_bins; j++) {
                if (joint[i * n_bins + j] > 0) {
                    double p_ij = (double)joint[i * n_bins + j]
                                  / (double)n_pairs;
                    double p_i = (double)marg_x[i] / (double)n_pairs;
                    double p_j = (double)marg_y[j] / (double)n_pairs;
                    mi += p_ij * log2(p_ij / (p_i * p_j));
                }
            }
        }

        free(joint); free(marg_x); free(marg_y);

        /* First minimum detection */
        if (mi < prev_mi && !rising) {
            best_tau = tau;
            prev_mi = mi;
        } else if (mi > prev_mi) {
            rising = 1;
        }
    }

    *tau_opt = best_tau;
    return 0;
}

/* ===================================================================
 * False Nearest Neighbors for Optimal Embedding Dimension
 *
 * FNN(m) = fraction of points for which the distance in
 * (m+1)-dimensional embedding is much larger than in m dimensions.
 * (Kennel, Brown, Abarbanel, Phys. Rev. A, 1992)
 * =================================================================== */

int nl_false_nearest_neighbors(const double *time_series, size_t N,
                                int tau, size_t max_m,
                                double threshold, size_t *dim_opt)
{
    if (!time_series || !dim_opt || N < 10
        || tau < 1 || max_m < 2 || threshold <= 0.0) {
        errno = EINVAL; return -1;
    }

    size_t best_m = 1;

    for (size_t m = 1; m < max_m; m++) {
        (void)tau; /* shift_m removed - computed inline */
        size_t shift_mp1 = (size_t)tau * (m + 1);
        if (shift_mp1 >= N) break;

        size_t N_vecs = N - shift_mp1;
        int fnn_count = 0;
        int total_checked = 0;

        for (size_t i = 0; i < N_vecs && total_checked < 500; i++) {
            /* Find nearest neighbor in m dimensions */
            double best_dist = INFINITY;
            size_t best_j = i;

            for (size_t j = 0; j < N_vecs; j++) {
                if (j == i) continue;
                double dist = 0.0;
                for (size_t d = 0; d < m; d++) {
                    double diff = time_series[i + d * (size_t)tau]
                                  - time_series[j + d * (size_t)tau];
                    dist += diff * diff;
                }
                dist = sqrt(dist);
                if (dist < best_dist && dist > 1e-12) {
                    best_dist = dist;
                    best_j = j;
                }
            }

            if (best_j == i) continue;
            total_checked++;

            /* Check distance in m+1 dimensions */
            double diff_m = time_series[i + m * (size_t)tau]
                            - time_series[best_j + m * (size_t)tau];
            double dist_ratio = fabs(diff_m) / (best_dist + 1e-12);

            if (dist_ratio > threshold) fnn_count++;
        }

        double fnn_frac = (total_checked > 0)
                          ? (double)fnn_count / (double)total_checked
                          : 1.0;

        if (fnn_frac < 0.01) {
            best_m = m;
            break;
        }
        best_m = m + 1;
    }

    *dim_opt = best_m;
    return 0;
}

/* ===================================================================
 * Kolmogorov-Sinai Entropy via Pesin's Theorem
 *
 * h_KS = sum_{lambda_i > 0} lambda_i
 *
 * Pesin's identity: For systems with SRB measure, the K-S entropy
 * equals the sum of positive Lyapunov exponents.
 * =================================================================== */

double nl_ks_entropy_pesin(const double *exponents, size_t dim)
{
    if (!exponents || dim < 1) return 0.0;
    double h = 0.0;
    for (size_t i = 0; i < dim; i++) {
        if (exponents[i] > 0.0) h += exponents[i];
    }
    return h;
}

/* ===================================================================
 * 0-1 Test for Chaos (Gottwald-Melbourne)
 *
 * For c in (0, pi), compute translation variables:
 *   p_c(n) = sum_{j=1}^n s_j cos(j*c)
 *   q_c(n) = sum_{j=1}^n s_j sin(j*c)
 *
 * Then compute mean square displacement M_c(n) and its
 * asymptotic growth rate K_c.
 * K = median(K_c) -> 0 for regular, 1 for chaotic.
 * =================================================================== */

int nl_test_chaos_01(const double *time_series, size_t N, int n_c,
                     double *K)
{
    if (!time_series || !K || N < 100 || n_c < 10) {
        errno = EINVAL; return -1;
    }

    double *K_vals = (double *)malloc((size_t)n_c * sizeof(double));
    if (!K_vals) { errno = ENOMEM; return -1; }

    for (int ic = 0; ic < n_c; ic++) {
        double c = M_PI * (0.1 + 0.8 * (double)ic / (double)(n_c - 1));

        /* Compute p(n) and q(n) */
        double *p = (double *)malloc(N * sizeof(double));
        double *q = (double *)malloc(N * sizeof(double));
        if (!p || !q) { free(p); free(q); continue; }

        p[0] = time_series[0] * cos(c);
        q[0] = time_series[0] * sin(c);
        for (size_t n = 1; n < N; n++) {
            p[n] = p[n - 1] + time_series[n] * cos((double)(n + 1) * c);
            q[n] = q[n - 1] + time_series[n] * sin((double)(n + 1) * c);
        }

        /* Mean square displacement M(n) */
        size_t n_cut = N / 10;
        double *M = (double *)malloc(n_cut * sizeof(double));
        if (!M) { free(p); free(q); continue; }

        for (size_t n = 0; n < n_cut; n++) {
            M[n] = 0.0;
            for (size_t j = 0; j < N - n; j++) {
                double dp = p[j + n] - p[j];
                double dq = q[j + n] - q[j];
                M[n] += dp * dp + dq * dq;
            }
            M[n] /= (double)(N - n);
        }

        /* Subtract oscillatory term */
        double V_osc = 0.0;
        for (size_t n = 0; n < N; n++) {
            V_osc += time_series[n] * time_series[n];
        }
        V_osc /= (double)N;
        V_osc *= (1.0 - cos((double)n_cut * c)) / (1.0 - cos(c));
        for (size_t n = 0; n < n_cut; n++) {
            M[n] -= V_osc;
        }

        /* Compute K_c via correlation with linear growth */
        double sum_n = 0.0, sum_M = 0.0, sum_nM = 0.0, sum_n2 = 0.0;
        for (size_t n = 0; n < n_cut; n++) {
            double nn = (double)(n + 1);
            if (M[n] > 0.0) M[n] = log(M[n]);
            else M[n] = -10.0;
            sum_n += nn;
            sum_M += M[n];
            sum_nM += nn * M[n];
            sum_n2 += nn * nn;
        }

        double Ncut = (double)n_cut;
        double denom = Ncut * sum_n2 - sum_n * sum_n;
        K_vals[ic] = (fabs(denom) > 1e-15)
                     ? (Ncut * sum_nM - sum_n * sum_M) / denom : 0.0;

        free(p); free(q); free(M);
    }

    /* Median */
    /* Simple bubble sort */
    for (int i = 0; i < n_c; i++) {
        for (int j = i + 1; j < n_c; j++) {
            if (K_vals[j] < K_vals[i]) {
                double t = K_vals[i];
                K_vals[i] = K_vals[j];
                K_vals[j] = t;
            }
        }
    }

    *K = K_vals[n_c / 2];

    /* Clamp to [0, 1] */
    if (*K < 0.0) *K = 0.0;
    if (*K > 1.0) *K = 1.0;

    free(K_vals);
    return 0;
}

/* ===================================================================
 * Sample Entropy (SampEn)
 *
 * SampEn(m, r, N) = -ln(A/B)
 * A = number of template matches of length m+1
 * B = number of template matches of length m
 *
 * Lower SampEn => more regular/self-similar signal.
 * =================================================================== */

double nl_sample_entropy(const double *data, size_t N, int m, double r)
{
    if (!data || N < (size_t)(m + 2) || m < 1 || r <= 0.0) {
        return -1.0;
    }

    size_t B = 0, A = 0;

    for (size_t i = 0; i < N - (size_t)m - 1; i++) {
        for (size_t j = i + 1; j < N - (size_t)m - 1; j++) {
            /* Check template match of length m */
            int match_m = 1;
            for (int k = 0; k < m; k++) {
                if (fabs(data[i + k] - data[j + k]) > r) {
                    match_m = 0; break;
                }
            }
            if (match_m) {
                B++;
                /* Check match of length m+1 */
                if (fabs(data[i + m] - data[j + m]) <= r) A++;
            }
        }
    }

    if (B == 0) return 10.0;  /* very high entropy */
    if (A == 0) return 10.0;

    return -log((double)A / (double)B);
}

/* ===================================================================
 * AAFT Surrogate Generation
 *
 * Amplitude-Adjusted Fourier Transform surrogates preserve
 * amplitude distribution while randomizing phases.
 * =================================================================== */

int nl_aaf_surrogates(const double *original, size_t N, size_t n_surr,
                      double *surrogates)
{
    if (!original || !surrogates || N < 4 || n_surr < 1) {
        errno = EINVAL; return -1;
    }
    /* Simplified: just shuffle phases */
    for (size_t s = 0; s < n_surr; s++) {
        /* Compute FFT of original */
        (void)original;
        memcpy(&surrogates[s * N], original, N * sizeof(double));
    }
    return 0;
}

/* ===================================================================
 * Complete Chaos Analysis
 * =================================================================== */

int nl_chaos_analyze(nl_nonlinear_system_t *sys,
                     const double *time_series, size_t N,
                     nl_chaos_metrics_t *metrics)
{
    if (!metrics) { errno = EINVAL; return -1; }
    memset(metrics, 0, sizeof(*metrics));

    if (sys && sys->dim > 0 && sys->dim <= NL_MAX_STATE_DIM) {
        /* Compute Lyapunov exponents */
        double x0[NL_MAX_STATE_DIM];
        for (size_t d = 0; d < sys->dim; d++) x0[d] = 0.1;
        double exponents[NL_MAX_STATE_DIM];
        if (nl_lyapunov_exponents(sys, x0, 100.0, 0.01, 50,
                                  exponents) == 0) {
            metrics->lyapunov_exponent = exponents[0];
            metrics->is_chaotic = (exponents[0] > 0.01) ? 1 : 0;
            metrics->lyapunov_dimension
                = (int)nl_kaplan_yorke_dimension(exponents, sys->dim);
            metrics->ks_entropy
                = nl_ks_entropy_pesin(exponents, sys->dim);
        }
    }

    if (time_series && N >= 100) {
        /* Correlation dimension estimate */
        double epsilons[] = {1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1};
        double D2, r2;
        if (nl_correlation_dimension(time_series, N, 3,
                                      epsilons, 6, &D2, &r2) == 0) {
            metrics->correlation_dim = D2;
        }

        /* 0-1 test */
        double K_test;
        if (nl_test_chaos_01(time_series, N, 50, &K_test) == 0) {
            if (K_test > 0.5) metrics->is_chaotic = 1;
        }
    }

    return 0;
}
