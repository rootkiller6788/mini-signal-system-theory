/**
 * @file sysid_nonlinear.c
 * @brief Nonlinear system identification methods (L8 Advanced Topics)
 *
 * Implements NARX (Nonlinear ARX), Hammerstein, and Wiener model structures.
 *
 * Knowledge points:
 *   L8 - NARX: y(k) = f(y(k-1),...,y(k-na), u(k-nk),...,u(k-nk-nb+1))
 *   L8 - Hammerstein: static nonlinearity followed by linear dynamics
 *   L8 - Wiener: linear dynamics followed by static nonlinearity
 *   L8 - Polynomial basis expansion for nonlinear regression
 *   L8 - Sigmoid network (simple neural network) for NARX
 *
 * References:
 *   Billings, S.A. (2013) "Nonlinear System Identification: NARMAX Methods"
 *   Sjöberg et al. (1995) "Nonlinear black-box modeling in system identification"
 *   Haber & Unbehauen (1990) "Structure identification of nonlinear dynamic systems"
 */

#include "sysid_types.h"
#include "sysid_models.h"
#include "sysid_estimation.h"
#include "sysid_nonlinear.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/*---------------------------------------------------------------------------
 * NARX: Nonlinear ARX model
 * y(k) = f(φ(k); θ)
 * where φ(k) = [y(k-1),...,y(k-na), u(k-nk),...,u(k-nk-nb+1)]
 *
 * f is parameterized as a polynomial expansion or a sigmoid network.
 *---------------------------------------------------------------------------*/

/** Polynomial NARX model:
 *  f(φ) = θ₀ + Σ θᵢ φᵢ + Σ θᵢⱼ φᵢ φⱼ + ... (up to given polynomial degree)
 */
typedef struct {
    int    na;           /**< AR order */
    int    nb;           /**< Input order */
    int    nk;           /**< Input delay */
    int    degree;       /**< Polynomial degree (1=linear, 2=quadratic, 3=cubic) */
    int    n_regressors; /**< Number of basis regressors after expansion */
    int    n_linear;     /**< Number of linear terms (na + nb) */
    double *theta;       /**< Polynomial coefficients [n_regressors] */
    double Ts;
} sysid_narx_poly_t;

/**
 * Compute the number of polynomial basis terms for given number of linear
 * regressors d and polynomial degree p:
 *   n_terms = C(d+p, p) = (d+p)! / (d! * p!)
 * Includes constant term.
 */
int sysid_narx_poly_nterms(int d, int p) {
    /* Compute binomial coefficient C(d+p, p) using multiplicative formula */
    if (p == 0) return 1;
    int n = d + p;
    int k = (p < d) ? p : d;
    long long result = 1;
    for (int i = 1; i <= k; i++) {
        result = result * (n - k + i) / i;
    }
    return (int)result;
}

/**
 * Generate polynomial basis expansion of regressor vector φ.
 * For d=2, p=2: φ = [x₁, x₂] → basis = [1, x₁, x₂, x₁², x₁x₂, x₂²]
 *
 * @param phi       Linear regressors [d]
 * @param d         Number of linear regressors
 * @param degree    Polynomial degree
 * @param basis     Output basis vector [sysid_narx_poly_nterms(d, degree)]
 * @param n_basis   Output: actual number of basis terms
 * @return          0 on success
 */
int sysid_narx_poly_expand(const double *phi, int d, int degree,
                            double *basis, int *n_basis) {
    if (!phi || !basis || !n_basis || d <= 0 || degree < 1) return -1;
    int nb_terms = sysid_narx_poly_nterms(d, degree);
    *n_basis = nb_terms;

    /* For degree 1: just linear terms + constant */
    basis[0] = 1.0;
    for (int i = 0; i < d; i++) basis[1 + i] = phi[i];

    if (degree >= 2) {
        /* Quadratic cross terms */
        int idx = 1 + d; /* after linear terms */
        for (int i = 0; i < d; i++) {
            for (int j = i; j < d; j++) {
                basis[idx++] = phi[i] * phi[j];
            }
        }
    }

    if (degree >= 3) {
        /* Cubic terms (partial implementation: include up to cubic for small d) */
        int idx = 1 + d + d * (d + 1) / 2;
        for (int i = 0; i < d && idx < nb_terms; i++) {
            for (int j = i; j < d && idx < nb_terms; j++) {
                for (int k = j; k < d && idx < nb_terms; k++) {
                    basis[idx++] = phi[i] * phi[j] * phi[k];
                }
            }
        }
    }

    return 0;
}

/**
 * Estimate NARX polynomial model from I/O data using LS on expanded regressors.
 *
 * @param u         Input data [N]
 * @param y         Output data [N]
 * @param N         Number of samples
 * @param na        AR order
 * @param nb        Input order
 * @param nk        Input delay
 * @param degree    Polynomial degree (1, 2, or 3)
 * @param theta     Output: polynomial coefficients
 * @param n_terms   Output: number of terms
 * @param sigma2    Output: noise variance estimate
 * @return          0 on success
 */
int sysid_narx_poly_estimate(const double *u, const double *y, int N,
                              int na, int nb, int nk, int degree,
                              double *theta, int *n_terms, double *sigma2) {
    if (!u || !y || !theta || !n_terms || N <= 0) return -1;
    int d = na + nb; /* number of linear regressors */
    int nb_terms = sysid_narx_poly_nterms(d, degree);
    *n_terms = nb_terms;

    int start_row = nk + nb;
    if (na > start_row) start_row = na;
    int nrows = N - start_row;
    if (nrows <= nb_terms) return -1;

    /* Build expanded regression matrix */
    double *Phi_exp = (double *)calloc((size_t)(nrows * nb_terms), sizeof(double));
    double *Y_vec = (double *)calloc((size_t)nrows, sizeof(double));
    double *phi_lin = (double *)calloc((size_t)d, sizeof(double));
    double *basis_row = (double *)calloc((size_t)nb_terms, sizeof(double));
    if (!Phi_exp || !Y_vec || !phi_lin || !basis_row) {
        free(Phi_exp); free(Y_vec); free(phi_lin); free(basis_row);
        return -1;
    }

    for (int k = start_row; k < N; k++) {
        /* Build linear regressor vector */
        int col = 0;
        for (int i = 1; i <= na; i++) phi_lin[col++] = -y[k - i];
        for (int i = 0; i < nb; i++) phi_lin[col++] = u[k - nk - i];

        int num_basis;
        sysid_narx_poly_expand(phi_lin, d, degree, basis_row, &num_basis);

        for (int j = 0; j < num_basis; j++) {
            Phi_exp[(k - start_row) * nb_terms + j] = basis_row[j];
        }
        Y_vec[k - start_row] = y[k];
    }

    /* LS solve */
    extern int sysid_ls_solve(const double *Phi, const double *Y, int nrows,
                               int ncols, double *theta, double *cov, double *sigma2);
    int ret = sysid_ls_solve(Phi_exp, Y_vec, nrows, nb_terms, theta, NULL, sigma2);

    free(Phi_exp); free(Y_vec); free(phi_lin); free(basis_row);
    return ret;
}

/**
 * Predict NARX polynomial model output.
 */
int sysid_narx_poly_predict(const double *u, const double *y, int N,
                             int na, int nb, int nk, int degree,
                             const double *theta, double *ypred) {
    if (!u || !y || !theta || !ypred || N <= 0) return -1;
    int d = na + nb;
    int nb_terms = sysid_narx_poly_nterms(d, degree);

    int start_row = nk + nb;
    if (na > start_row) start_row = na;

    double *phi_lin = (double *)calloc((size_t)d, sizeof(double));
    double *basis_row = (double *)calloc((size_t)nb_terms, sizeof(double));
    if (!phi_lin || !basis_row) { free(phi_lin); free(basis_row); return -1; }

    for (int k = 0; k < N; k++) {
        if (k < start_row) {
            ypred[k] = 0.0;
            continue;
        }
        int col = 0;
        for (int i = 1; i <= na; i++) phi_lin[col++] = -y[k - i];
        for (int i = 0; i < nb; i++) phi_lin[col++] = u[k - nk - i];

        int n_basis;
        sysid_narx_poly_expand(phi_lin, d, degree, basis_row, &n_basis);

        double pred = 0.0;
        for (int j = 0; j < n_basis; j++) pred += theta[j] * basis_row[j];
        ypred[k] = pred;
    }

    free(phi_lin); free(basis_row);
    return 0;
}

/*---------------------------------------------------------------------------
 * Sigmoid Network NARX (simple 1-hidden-layer neural network)
 * f(φ) = Σ wⱼ σ(vⱼᵀ φ + bⱼ) + w₀
 * where σ(x) = 1/(1+e^{-x}) is the sigmoid function.
 *---------------------------------------------------------------------------*/

/* sysid_sig_net_t defined in sysid_nonlinear.h */

/** Sigmoid activation */
static inline double sigmoid(double x) {
    if (x > 50.0) return 1.0;
    if (x < -50.0) return 0.0;
    return 1.0 / (1.0 + exp(-x));
}

/**
 * Initialize sigmoid network for NARX with random small weights.
 */
int sysid_sig_net_alloc(sysid_sig_net_t *net, int na, int nb, int nk,
                         int n_hidden, double Ts, unsigned int seed) {
    if (!net || na < 0 || nb < 0 || n_hidden <= 0) return -1;
    memset(net, 0, sizeof(sysid_sig_net_t));
    net->na = na; net->nb = nb; net->nk = nk;
    net->n_hidden = n_hidden;
    net->n_input = na + nb;
    net->Ts = Ts;

    int n_in = net->n_input;
    net->W1 = (double *)calloc((size_t)(n_hidden * n_in), sizeof(double));
    net->b1 = (double *)calloc((size_t)n_hidden, sizeof(double));
    net->W2 = (double *)calloc((size_t)n_hidden, sizeof(double));
    if (!net->W1 || !net->b1 || !net->W2) {
        free(net->W1); free(net->b1); free(net->W2); return -1;
    }

    /* Xavier/He initialization with small weights */
    unsigned int state = seed;
    double scale = sqrt(2.0 / (double)n_in);
    for (int i = 0; i < n_hidden * n_in; i++) {
        /* Simple LCG for [-scale, scale] */
        state = (1103515245U * state + 12345U) & 0x7fffffffU;
        net->W1[i] = scale * (2.0 * (double)state / (double)0x80000000U - 1.0);
    }
    net->b2 = 0.0;
    for (int i = 0; i < n_hidden; i++) {
        state = (1103515245U * state + 12345U) & 0x7fffffffU;
        net->W2[i] = scale * (2.0 * (double)state / (double)0x80000000U - 1.0);
    }

    return 0;
}

void sysid_sig_net_free(sysid_sig_net_t *net) {
    if (!net) return;
    free(net->W1); free(net->b1); free(net->W2);
    net->W1 = NULL; net->b1 = NULL; net->W2 = NULL;
}

/**
 * Forward pass of sigmoid network: f(φ) = W2 ᵀ σ(W1 φ + b1) + b2
 */
double sysid_sig_net_forward(const sysid_sig_net_t *net, const double *phi) {
    if (!net || !phi) return 0.0;
    int n_h = net->n_hidden, n_in = net->n_input;
    double output = net->b2;
    for (int j = 0; j < n_h; j++) {
        double z = net->b1[j];
        for (int i = 0; i < n_in; i++) z += net->W1[j * n_in + i] * phi[i];
        output += net->W2[j] * sigmoid(z);
    }
    return output;
}

/**
 * Train sigmoid network using stochastic gradient descent with backpropagation.
 */
int sysid_sig_net_train(sysid_sig_net_t *net, const double *u, const double *y,
                         int N, int epochs, double lr) {
    if (!net || !u || !y || N <= 0) return -1;
    int na = net->na, nb = net->nb, nk = net->nk;
    int n_in = net->n_input, n_h = net->n_hidden;

    int start_row = nk + nb;
    if (na > start_row) start_row = na;

    double *phi = (double *)calloc((size_t)n_in, sizeof(double));
    double *hidden_z = (double *)calloc((size_t)n_h, sizeof(double));
    double *hidden_a = (double *)calloc((size_t)n_h, sizeof(double));
    double *delta_W1 = (double *)calloc((size_t)(n_h * n_in), sizeof(double));
    double *delta_b1 = (double *)calloc((size_t)n_h, sizeof(double));
    double *delta_W2 = (double *)calloc((size_t)n_h, sizeof(double));

    for (int epoch = 0; epoch < epochs; epoch++) {
        for (int k = start_row; k < N; k++) {
            /* Build regressor */
            int col = 0;
            for (int i = 1; i <= na; i++) phi[col++] = -y[k - i];
            for (int i = 0; i < nb; i++) phi[col++] = u[k - nk - i];

            /* Forward pass */
            for (int j = 0; j < n_h; j++) {
                hidden_z[j] = net->b1[j];
                for (int i = 0; i < n_in; i++) hidden_z[j] += net->W1[j * n_in + i] * phi[i];
                hidden_a[j] = sigmoid(hidden_z[j]);
            }
            double ypred = net->b2;
            for (int j = 0; j < n_h; j++) ypred += net->W2[j] * hidden_a[j];

            double err = y[k] - ypred;

            /* Backpropagation */
            /* dL/dW2 = -err * hidden_a */
            /* dL/db2 = -err */
            /* dL/dW1 = -err * W2_j * σ'(z) * φ_i */
            /* dL/db1 = -err * W2_j * σ'(z) */

            double db2_delta = lr * err;

            for (int j = 0; j < n_h; j++) {
                double dW2 = lr * err * hidden_a[j];
                net->W2[j] += dW2;
                delta_W2[j] = dW2;

                double sigma_prime = hidden_a[j] * (1.0 - hidden_a[j]);
                double d_hidden = err * net->W2[j] * sigma_prime;

                for (int i = 0; i < n_in; i++) {
                    net->W1[j * n_in + i] += lr * d_hidden * phi[i];
                }
                net->b1[j] += lr * d_hidden;
            }
            net->b2 += db2_delta;
        }
    }

    free(phi); free(hidden_z); free(hidden_a);
    free(delta_W1); free(delta_b1); free(delta_W2);
    return 0;
}

/**
 * Predict using trained sigmoid network.
 */
int sysid_sig_net_predict(const sysid_sig_net_t *net, const double *u,
                           const double *y, int N, double *ypred) {
    if (!net || !u || !y || !ypred || N <= 0) return -1;
    int na = net->na, nb = net->nb, nk = net->nk;
    int n_in = net->n_input;

    int start_row = nk + nb;
    if (na > start_row) start_row = na;

    double *phi = (double *)calloc((size_t)n_in, sizeof(double));
    if (!phi) return -1;

    for (int k = 0; k < N; k++) {
        if (k < start_row) { ypred[k] = 0.0; continue; }
        int col = 0;
        for (int i = 1; i <= na; i++) phi[col++] = -y[k - i];
        for (int i = 0; i < nb; i++) phi[col++] = u[k - nk - i];
        ypred[k] = sysid_sig_net_forward(net, phi);
    }
    free(phi);
    return 0;
}

/*---------------------------------------------------------------------------
 * Hammerstein model: static nonlinearity → linear dynamics
 *
 * Structure:
 *   v(k) = g(u(k); α)       [static nonlinearity, e.g., polynomial]
 *   y(k) = [B(q)/F(q)] v(k) + e(k)  [linear OE dynamics]
 *
 * Identification approach (iterative / Narendra-Gallman, 1966):
 *   1. Fix linear part, estimate nonlinearity via LS
 *   2. Fix nonlinearity, estimate linear part via PEM
 *   3. Iterate until convergence
 *---------------------------------------------------------------------------*/

/**
 * Estimate Hammerstein model (simplified: polynomial NL + FIR linear part)
 *
 * @param u        Input data [N]
 * @param y        Output data [N]
 * @param N        Samples
 * @param nb       FIR length
 * @param nk       Delay
 * @param poly_order  Polynomial order for static nonlinearity
 * @param alpha    Output: polynomial coefficients [poly_order]
 * @param b        Output: FIR coefficients [nb]
 * @param max_iter Maximum iterations
 * @return         0 on success
 */
int sysid_hammerstein_estimate(const double *u, const double *y, int N,
                                int nb, int nk, int poly_order,
                                double *alpha, double *b,
                                double *sigma2, int max_iter) {
    if (!u || !y || !alpha || !b || N <= 0 || nb <= 0 || poly_order <= 0) return -1;

    /* Initialize: b = simple impulse estimate */
    double *v = (double *)calloc((size_t)N, sizeof(double));
    if (!v) return -1;

    /* Copy u as initial v (identity nonlinearity) */
    memcpy(v, u, (size_t)N * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        /* Step 1: Given v, estimate linear part (FIR) */
        int nrows_fir = N - nk - nb + 1;
        if (nrows_fir <= 0) continue;

        double *Phi_fir = (double *)calloc((size_t)(nrows_fir * nb), sizeof(double));
        double *Y_fir = (double *)calloc((size_t)nrows_fir, sizeof(double));
        if (Phi_fir && Y_fir) {
            for (int k = nk + nb - 1; k < N; k++) {
                int row = k - (nk + nb - 1);
                for (int i = 0; i < nb; i++) {
                    Phi_fir[row * nb + i] = v[k - nk - i];
                }
                Y_fir[row] = y[k];
            }
            double *theta_b = (double *)calloc((size_t)nb, sizeof(double));
            extern int sysid_ls_solve(const double *Phi, const double *Y, int nrows,
                                       int ncols, double *theta, double *cov, double *sigma2);
            sysid_ls_solve(Phi_fir, Y_fir, nrows_fir, nb, theta_b, NULL, NULL);
            memcpy(b, theta_b, (size_t)nb * sizeof(double));
            free(theta_b);
            free(Phi_fir); free(Y_fir);
        }

        /* Step 2: Given linear part, estimate nonlinearity.
         * We have: y(k) ≈ Σ b_i * g(u(k-nk-i))
         * This is nonlinear in g; use polynomial parameterization:
         * g(u) = α₁ u + α₂ u² + ... + α_p u^p
         * Then y(k) ≈ Σ b_i * Σ α_j * u^j(k-nk-i) = Σ α_j * Σ b_i u^j(k-nk-i)
         * This is LINEAR in α!
         */

        int nrows_nl = N - nk - nb + 1;
        double *Phi_nl = (double *)calloc((size_t)(nrows_nl * poly_order), sizeof(double));
        if (Phi_nl) {
            for (int k = nk + nb - 1; k < N; k++) {
                int row = k - (nk + nb - 1);
                for (int p = 1; p <= poly_order; p++) {
                    double sum = 0.0;
                    for (int i = 0; i < nb; i++) {
                        int idx = k - nk - i;
                        if (idx >= 0) {
                            sum += b[i] * pow(u[idx], (double)p);
                        }
                    }
                    Phi_nl[row * poly_order + (p - 1)] = sum;
                }
            }
            double *theta_alpha = (double *)calloc((size_t)poly_order, sizeof(double));
            sysid_ls_solve(Phi_nl, Y_fir, nrows_nl, poly_order,
                          theta_alpha, NULL, NULL);
            memcpy(alpha, theta_alpha, (size_t)poly_order * sizeof(double));
            free(theta_alpha);
        }
        free(Phi_nl);

        /* Recompute v from updated nonlinearity */
        for (int k = 0; k < N; k++) {
            double vu = 0.0;
            for (int p = 1; p <= poly_order; p++) {
                vu += alpha[p - 1] * pow(u[k], (double)p);
            }
            v[k] = vu;
        }
    }

    if (sigma2) {
        /* Compute final loss */
        double ss = 0.0;
        for (int k = nk + nb - 1; k < N; k++) {
            double yh = 0.0;
            for (int i = 0; i < nb; i++) {
                yh += b[i] * v[k - nk - i];
            }
            ss += (y[k] - yh) * (y[k] - yh);
        }
        *sigma2 = ss / (double)(N - nk - nb);
    }

    free(v);
    return 0;
}

/*---------------------------------------------------------------------------
 * Wiener model: linear dynamics → static nonlinearity
 *
 * Structure:
 *   w(k) = [B(q)/F(q)] u(k)        [linear OE dynamics]
 *   y(k) = h(w(k); β) + e(k)       [static nonlinearity]
 *
 * Identification (simplified: FIR + polynomial output NL)
 *---------------------------------------------------------------------------*/

/**
 * Estimate Wiener model (simplified).
 */
int sysid_wiener_estimate(const double *u, const double *y, int N,
                           int nb, int nk, int poly_order,
                           double *b, double *beta, double *sigma2) {
    if (!u || !y || !b || !beta || N <= 0 || nb <= 0) return -1;

    /* Wiener: estimate intermediate signal w from inverse NL, then LS for FIR.
     * Simplified iterative approach:
     * 1. Assume identity NL, estimate FIR
     * 2. Given FIR, estimate output NL
     * 3. Given NL, estimate FIR, etc.
     */

    /* Step 1: FIR from u→y (assuming identity NL) */
    int nrows = N - nk - nb + 1;
    if (nrows <= 0) return -1;
    double *Phi = (double *)calloc((size_t)(nrows * nb), sizeof(double));
    double *Y_vec = (double *)calloc((size_t)nrows, sizeof(double));
    double *w_est = (double *)calloc((size_t)N, sizeof(double));
    if (!Phi || !Y_vec || !w_est) {
        free(Phi); free(Y_vec); free(w_est); return -1;
    }

    for (int k = nk + nb - 1; k < N; k++) {
        int row = k - (nk + nb - 1);
        for (int i = 0; i < nb; i++) Phi[row * nb + i] = u[k - nk - i];
        Y_vec[row] = y[k];
    }

    double *theta_b = (double *)calloc((size_t)nb, sizeof(double));
    extern int sysid_ls_solve(const double *Phi, const double *Y, int nrows,
                               int ncols, double *theta, double *cov, double *sigma2);
    sysid_ls_solve(Phi, Y_vec, nrows, nb, theta_b, NULL, NULL);
    memcpy(b, theta_b, (size_t)nb * sizeof(double));

    /* Step 2: Estimate w = Σ b_i u (linear prediction), then fit NL from w to y */
    for (int k = 0; k < N; k++) {
        w_est[k] = 0.0;
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) w_est[k] += b[i] * u[idx];
        }
    }

    /* Fit polynomial: y ≈ β₀ + β₁ w + β₂ w² + ... + β_p w^p */
    int nrows_nl = N;
    double *Phi_nl = (double *)calloc((size_t)(nrows_nl * (poly_order + 1)), sizeof(double));
    if (Phi_nl) {
        for (int k = 0; k < nrows_nl; k++) {
            Phi_nl[k * (poly_order + 1)] = 1.0; /* constant */
            for (int p = 1; p <= poly_order; p++) {
                Phi_nl[k * (poly_order + 1) + p] = pow(w_est[k], (double)p);
            }
        }
        double *theta_beta = (double *)calloc((size_t)(poly_order + 1), sizeof(double));
        sysid_ls_solve(Phi_nl, y, N, poly_order + 1, theta_beta, NULL, sigma2);
        memcpy(beta, theta_beta, (size_t)(poly_order + 1) * sizeof(double));
        free(theta_beta);
    }
    free(Phi_nl);

    free(Phi); free(Y_vec); free(w_est); free(theta_b);
    return 0;
}
