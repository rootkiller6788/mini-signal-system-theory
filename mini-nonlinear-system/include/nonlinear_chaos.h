/**
 * @file nonlinear_chaos.h
 * @brief Chaos detection, characterization, and analysis.
 *
 * Implements methods for detecting and quantifying chaotic behavior
 * in nonlinear dynamical systems, including Lyapunov exponent
 * computation, fractal dimension estimation, and entropy measures.
 *
 * L2: Core Concepts — deterministic chaos, strange attractors,
 *     sensitive dependence on initial conditions
 * L5: Algorithms — Lyapunov exponent computation (Benettin),
 *     correlation dimension (Grassberger-Procaccia),
 *     delay embedding (Takens)
 * L6: Canonical Problems — Lorenz chaos, Rossler chaos,
 *     period-doubling route, intermittency
 * L8: Advanced Topics — fractal dimension, Kolmogorov-Sinai entropy,
 *     surrogate data testing, recurrence quantification
 *
 * Course Mapping:
 *   - MIT 18.385: Nonlinear Dynamics and Chaos
 *   - Cornell MAE 6780: Nonlinear Dynamics
 *   - Tsinghua: 非线性动力学
 *
 * Reference:
 *   - Strogatz, "Nonlinear Dynamics and Chaos", Perseus, 1994
 *   - Kantz & Schreiber, "Nonlinear Time Series Analysis", 2nd ed.,
 *     Cambridge, 2004
 *   - Sprott, "Chaos and Time-Series Analysis", Oxford, 2003
 */

#ifndef NONLINEAR_CHAOS_H
#define NONLINEAR_CHAOS_H

#include "nonlinear_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * L5: Lyapunov Exponents
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute all Lyapunov exponents using Benettin's algorithm (L5).
 *
 * The Lyapunov exponents λ₁ ≥ λ₂ ≥ ... ≥ λ_n characterize the average
 * exponential rates of divergence (or convergence) of nearby trajectories.
 *
 * If λ₁ > 0, the system exhibits sensitive dependence on initial
 * conditions, a hallmark of chaos.
 *
 * The algorithm (Benettin et al., 1980):
 * 1. Integrate the system and the variational equations simultaneously
 * 2. Periodically re-orthonormalize the deviation vectors via Gram-Schmidt
 * 3. The Lyapunov exponents are the time-averaged logarithms of the
 *    growth factors of each orthonormal vector
 *
 * Complexity: O(T/dt · dim² · cost(f)) where T is total integration time.
 *
 * Reference: Benettin, Galgani, Giorgilli, Strelcyn,
 * "Lyapunov Characteristic Exponents for Smooth Dynamical Systems...",
 * Meccanica, 1980.
 *
 * @param sys         Nonlinear system.
 * @param x0          Initial condition (length = dim).
 * @param T_total     Total integration time.
 * @param dt          Integration step.
 * @param n_renorm    Number of re-orthonormalization steps.
 * @param exponents   Output: Lyapunov exponents λ₁, ..., λ_dim.
 * @return 0 on success, -1 on error.
 */
int nl_lyapunov_exponents(nl_nonlinear_system_t *sys,
                          const double *x0, double T_total, double dt,
                          int n_renorm, double *exponents);

/**
 * @brief Compute the Kaplan-Yorke (Lyapunov) dimension (L5, L8).
 *
 * The Kaplan-Yorke dimension D_KY is defined as:
 *   D_KY = k + (Σ_{i=1}^k λ_i) / |λ_{k+1}|
 *
 * where k is the largest integer such that Σ_{i=1}^k λ_i ≥ 0,
 * and λ₁ ≥ λ₂ ≥ ... ≥ λ_n are the Lyapunov exponents.
 *
 * This provides an estimate of the fractal dimension of a strange
 * attractor and relates the dynamical (Lyapunov) properties to
 * geometric (dimension) properties.
 *
 * Reference: Kaplan & Yorke, "Chaotic behavior of multidimensional
 * difference equations", in Peitgen & Walther (eds.),
 * Functional Differential Equations..., Springer, 1979.
 *
 * @param exponents Lyapunov exponents (sorted descending).
 * @param dim       Number of exponents.
 * @return Kaplan-Yorke dimension (≥ 0).
 */
double nl_kaplan_yorke_dimension(const double *exponents, size_t dim);

/* ---------------------------------------------------------------------------
 * L5, L8: Correlation Dimension (Grassberger-Procaccia)
 * -------------------------------------------------------------------------*/

/**
 * @brief Compute correlation sum C(ε) for a time series (L5, L8).
 *
 * The correlation sum is the fraction of pairs of points whose
 * distance is less than ε:
 *
 *   C(ε) = 2/(N(N-1)) · Σ_{i=1}^{N} Σ_{j=i+1}^{N} Θ(ε - ||x_i - x_j||)
 *
 * where Θ is the Heaviside step function.
 *
 * For small ε, C(ε) ∼ ε^{D₂}, where D₂ is the correlation dimension.
 *
 * Reference: Grassberger & Procaccia, "Characterization of Strange
 * Attractors", Phys. Rev. Lett., 1983.
 *
 * @param data     Time series data (trajectory points).
 * @param N        Number of points.
 * @param dim      Embedding dimension.
 * @param epsilon  Distance threshold.
 * @return C(ε) ∈ [0, 1].
 */
double nl_correlation_sum(const double *data, size_t N, size_t dim,
                          double epsilon);

/**
 * @brief Estimate correlation dimension D₂ via log-log regression (L5, L8).
 *
 * Computes C(ε) for a range of ε values and performs linear regression
 * on log C(ε) vs log ε. The slope gives D₂.
 *
 * @param data      Time series (trajectory) data.
 * @param N         Number of points.
 * @param dim       Embedding dimension.
 * @param epsilons  Array of ε values.
 * @param n_eps     Number of ε values.
 * @param D2        Output: estimated correlation dimension.
 * @param r2        Output: R² of the log-log fit.
 * @return 0 on success, -1 on error.
 */
int nl_correlation_dimension(const double *data, size_t N, size_t dim,
                             const double *epsilons, size_t n_eps,
                             double *D2, double *r2);

/* ---------------------------------------------------------------------------
 * L5, L8: Time-Delay Embedding (Takens)
 * -------------------------------------------------------------------------*/

/**
 * @brief Create delay-coordinate embedding from scalar time series (L8).
 *
 * Takens' embedding theorem: a scalar time series {s_t} from a
 * d-dimensional dynamical system can be reconstructed in an
 * m-dimensional embedding space using delay coordinates:
 *
 *   x_t = (s_t, s_{t-τ}, s_{t-2τ}, ..., s_{t-(m-1)τ})
 *
 * provided m > 2d (for generic observables and almost every τ).
 *
 * The embedding preserves the topological invariants of the
 * original attractor (diffeomorphic reconstruction).
 *
 * Reference: Takens, "Detecting strange attractors in turbulence",
 * in Rand & Young (eds.), Dynamical Systems and Turbulence,
 * Springer Lecture Notes in Math, 1981.
 *
 * @param time_series Scalar time series.
 * @param N           Length of time series.
 * @param tau         Delay (in samples).
 * @param m           Embedding dimension.
 * @param embedded    Output: embedded vectors (N - (m-1)τ by m, row-major).
 * @return Number of embedded vectors, or -1 on error.
 */
int nl_delay_embedding(const double *time_series, size_t N,
                       int tau, size_t m, double *embedded);

/**
 * @brief Estimate optimal delay τ using mutual information (L8).
 *
 * The first minimum of the mutual information between s_t and s_{t-τ}
 * provides a heuristic for the optimal embedding delay τ.
 *
 * I(τ) = Σ_{i,j} p_{ij}(τ) log₂(p_{ij}(τ) / (p_i p_j))
 *
 * Using equiprobable partition of the amplitude range.
 *
 * Reference: Fraser & Swinney, "Independent coordinates for strange
 * attractors from mutual information", Phys. Rev. A, 1986.
 *
 * @param time_series Scalar time series.
 * @param N           Length.
 * @param max_tau     Maximum delay to search.
 * @param n_bins      Number of bins for mutual information.
 * @param tau_opt     Output: optimal delay.
 * @return 0 on success.
 */
int nl_mutual_information_delay(const double *time_series, size_t N,
                                 int max_tau, int n_bins,
                                 int *tau_opt);

/**
 * @brief Estimate optimal embedding dimension m using false nearest neighbors (L8).
 *
 * The percentage of false nearest neighbors drops to (near) zero when
 * the embedding dimension is sufficient to unfold the attractor.
 *
 * FNN(m) = fraction of points for which:
 *   |s_{t+mτ} - s_{r+mτ}| / R_m(t, r) > threshold
 *
 * where R_m is the distance in m-dimensional embedding space.
 *
 * Reference: Kennel, Brown, Abarbanel, "Determining embedding dimension
 * for phase-space reconstruction...", Phys. Rev. A, 1992.
 *
 * @param time_series Scalar time series.
 * @param N           Length.
 * @param tau         Delay (from mutual information).
 * @param max_m       Maximum embedding dimension to try.
 * @param threshold   FNN threshold (typically 10-15).
 * @param dim_opt     Output: optimal embedding dimension.
 * @return 0 on success.
 */
int nl_false_nearest_neighbors(const double *time_series, size_t N,
                                int tau, size_t max_m,
                                double threshold, size_t *dim_opt);

/* ---------------------------------------------------------------------------
 * L5, L8: Entropy and Chaos Quantification
 * -------------------------------------------------------------------------*/

/**
 * @brief Approximate Kolmogorov-Sinai (K-S) entropy via Pesin's theorem (L8).
 *
 * Pesin's theorem relates the K-S entropy h_KS to the Lyapunov exponents:
 *   h_KS = Σ_{λ_i > 0} λ_i
 *
 * This is an upper bound in general, but equality holds for systems
 * with an SRB (Sinai-Ruelle-Bowen) measure.
 *
 * Reference: Pesin, "Characteristic Lyapunov exponents and smooth
 * ergodic theory", Russian Math. Surveys, 1977.
 *
 * @param exponents Lyapunov exponents (sorted descending).
 * @param dim       Number of exponents.
 * @return K-S entropy estimate (≥ 0).
 */
double nl_ks_entropy_pesin(const double *exponents, size_t dim);

/**
 * @brief 0-1 test for chaos (Gottwald-Melbourne) (L8).
 *
 * The 0-1 test distinguishes regular (K ≈ 0) from chaotic (K ≈ 1)
 * dynamics directly from a time series, without requiring phase space
 * reconstruction or Lyapunov exponent computation.
 *
 * Method:
 * 1. Compute translation variables p_c(n) and q_c(n) for random c ∈ (0,π).
 * 2. Compute mean square displacement M_c(n).
 * 3. The growth rate K_c = lim_{n→∞} log M_c(n) / log n.
 * 4. K = median(K_c) over random c.
 *
 * Reference: Gottwald & Melbourne, "A new test for chaos in
 * deterministic systems", Proc. R. Soc. Lond. A, 2004.
 *
 * @param time_series Scalar time series.
 * @param N           Length of time series.
 * @param n_c         Number of random c values.
 * @param K           Output: test statistic (≈0 = regular, ≈1 = chaotic).
 * @return 0 on success.
 */
int nl_test_chaos_01(const double *time_series, size_t N, int n_c,
                     double *K);

/* ---------------------------------------------------------------------------
 * L8: Surrogate Data Testing
 * -------------------------------------------------------------------------*/

/**
 * @brief Generate amplitude-adjusted Fourier transform (AAFT) surrogates (L8).
 *
 * AAFT surrogates preserve the amplitude distribution and approximate
 * the power spectrum of the original data. They serve as null
 * hypothesis for testing nonlinearity/chaos: if a nonlinear statistic
 * computed on the original data is significantly different from that
 * on surrogates, nonlinearity is present.
 *
 * This is a two-step iterative procedure:
 * 1. Rescale Gaussian data to match the original's rank order
 * 2. Adjust power spectrum to match original
 * 3. Iterate to convergence
 *
 * Reference: Theiler, Eubank, Longtin, Galdrikian, Farmer,
 * "Testing for nonlinearity in time series...", Physica D, 1992.
 *
 * @param original   Original time series.
 * @param N          Length.
 * @param n_surr     Number of surrogates to generate.
 * @param surrogates Output: surrogates (n_surr × N, row-major).
 * @return 0 on success.
 */
int nl_aaf_surrogates(const double *original, size_t N, size_t n_surr,
                      double *surrogates);

/**
 * @brief Compute sample entropy (SampEn) for regularity quantification (L8).
 *
 * SampEn(m, r, N) = -ln(A/B) where A is the number of template matches
 * of length m+1 and B is the number of matches of length m.
 *
 * Lower values indicate more regular (predictable) time series.
 *
 * Reference: Richman & Moorman, "Physiological time-series analysis
 * using approximate entropy and sample entropy",
 * Am. J. Physiol. Heart Circ. Physiol., 2000.
 *
 * @param data  Time series.
 * @param N     Length.
 * @param m     Template length (typically 2).
 * @param r     Tolerance (typically 0.2 × std(data)).
 * @return Sample entropy value, or -1 on error.
 */
double nl_sample_entropy(const double *data, size_t N, int m, double r);

/**
 * @brief Perform complete chaos analysis on a time series (L8).
 *
 * Convenience function that combines Lyapunov exponent, correlation
 * dimension, K-S entropy, and the 0-1 test.
 *
 * @param sys          System (used for Lyapunov computation if available).
 * @param time_series  Time series data (can be NULL if using sys).
 * @param N            Length.
 * @param metrics      Output: chaos metrics struct.
 * @return 0 on success.
 */
int nl_chaos_analyze(nl_nonlinear_system_t *sys,
                     const double *time_series, size_t N,
                     nl_chaos_metrics_t *metrics);

#ifdef __cplusplus
}
#endif

#endif /* NONLINEAR_CHAOS_H */
