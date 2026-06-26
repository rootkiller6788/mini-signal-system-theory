/**
 * @file sysid_nonlinear.h
 * @brief Nonlinear system identification methods (L8 Advanced Topics)
 *
 * Knowledge points:
 *   L8 - NARX (Nonlinear ARX): y(k) = f(phi(k); theta) with polynomial basis
 *   L8 - Hammerstein model: static nonlinearity + linear dynamics
 *   L8 - Wiener model: linear dynamics + static nonlinearity
 *   L8 - Polynomial basis expansion for nonlinear regression
 *   L8 - Sigmoid network (1-hidden-layer NN) for NARX modeling
 *   L8 - Block-structured nonlinear model identification
 *
 * References:
 *   Billings, S.A. (2013) "Nonlinear System Identification: NARMAX Methods", Wiley
 *   Sjoberg et al. (1995) "Nonlinear black-box modeling in system identification"
 *   Haber & Unbehauen (1990) "Structure identification of nonlinear dynamic systems"
 *   Nelles, O. (2001) "Nonlinear System Identification", Springer
 *   Narendra & Gallman (1966) "An iterative method for the identification of
 *     nonlinear systems using a Hammerstein model"
 */

#ifndef SYSID_NONLINEAR_H
#define SYSID_NONLINEAR_H

#include "sysid_types.h"
#include "sysid_models.h"

/*===========================================================================
 * NARX — Nonlinear AutoRegressive with eXogenous inputs
 *
 * y(k) = f(phi(k); theta)
 * where phi(k) = [y(k-1), ..., y(k-na), u(k-nk), ..., u(k-nk-nb+1)]
 *
 * f is parameterized as either:
 *   1. Polynomial expansion: f(phi) = sum_k theta_k * basis_k(phi)
 *   2. Sigmoid network: f(phi) = sum_j w_j * sigma(v_j^T phi + b_j) + w0
 *===========================================================================*/

/*---------------------------------------------------------------------------
 * L8: Polynomial NARX model
 *
 * The polynomial basis includes all monomials up to a given degree:
 *   basis = {1, phi_1, ..., phi_d, phi_1*phi_1, phi_1*phi_2, ..., phi_d*phi_d,
 *            cubic terms, ...}
 *
 * Number of terms for d linear regressors and degree p:
 *   n_terms = C(d+p, p) = (d+p)! / (d! * p!)
 *
 * For p=2: n_terms = (d+2)(d+1)/2 = 1 + d + d(d+1)/2
 * For p=3: n_terms = (d+3)(d+2)(d+1)/6
 *
 * This is a powerful but parsimonious way to capture nonlinear dynamics
 * without requiring a priori structural knowledge (Billings, 2013).
 *---------------------------------------------------------------------------*/

/**
 * Compute the number of basis terms for polynomial expansion.
 *
 * @param d      Number of linear regressors (na + nb)
 * @param degree Polynomial degree (1 = linear, 2 = quadratic, 3 = cubic)
 * @return       Number of basis terms (includes constant)
 */
int sysid_narx_poly_nterms(int d, int degree);

/**
 * Expand a linear regressor vector into a polynomial basis.
 *
 * Example for d=2, p=2:
 *   phi = [x1, x2] -> basis = [1, x1, x2, x1*x1, x1*x2, x2*x2]
 *
 * @param phi      Linear regressors [d]
 * @param d        Number of linear regressors
 * @param degree   Polynomial degree (1, 2, or 3)
 * @param basis    Output: expanded basis vector [narx_poly_nterms(d, degree)]
 * @param n_basis  Output: actual number of basis terms
 * @return         0 on success, -1 on error
 */
int sysid_narx_poly_expand(const double *phi, int d, int degree,
                            double *basis, int *n_basis);

/**
 * Estimate polynomial NARX model from I/O data via least squares.
 *
 * Builds the expanded regression matrix and solves:
 *   theta = argmin sum_k (y(k) - sum_j theta_j * basis_j(phi(k)))^2
 *
 * @param u         Input data [N]
 * @param y         Output data [N]
 * @param N         Number of samples
 * @param na        AR order
 * @param nb        Input order
 * @param nk        Input delay
 * @param degree    Polynomial degree (1, 2, or 3)
 * @param theta     Output: polynomial coefficients [*n_terms]
 * @param n_terms   Output: number of terms in expansion
 * @param sigma2    Output: noise variance estimate (NULL to skip)
 * @return          0 on success, -1 on error
 */
int sysid_narx_poly_estimate(const double *u, const double *y, int N,
                              int na, int nb, int nk, int degree,
                              double *theta, int *n_terms, double *sigma2);

/**
 * Predict output using polynomial NARX model.
 *
 * @param u      Input data [N]
 * @param y      Output data [N] (provides past values for prediction)
 * @param N      Number of samples
 * @param na     AR order
 * @param nb     Input order
 * @param nk     Input delay
 * @param degree Polynomial degree
 * @param theta  Polynomial coefficients [narx_poly_nterms(na+nb, degree)]
 * @param ypred  Output: predicted output [N]
 * @return       0 on success, -1 on error
 */
int sysid_narx_poly_predict(const double *u, const double *y, int N,
                             int na, int nb, int nk, int degree,
                             const double *theta, double *ypred);

/*---------------------------------------------------------------------------
 * L8: Sigmoid Network NARX model
 *
 * A 1-hidden-layer feedforward neural network used as the nonlinear
 * function f in the NARX structure:
 *
 *   f(phi) = sum_{j=1}^{n_hidden} W2_j * sigma(W1_j^T * phi + b1_j) + b2
 *
 * where sigma(x) = 1 / (1 + e^{-x}) is the logistic sigmoid.
 *
 * This is a universal approximator (Hornik et al., 1989) and can represent
 * a much wider class of nonlinearities than polynomial models.
 *---------------------------------------------------------------------------*/

/** Sigmoid network model structure */
typedef struct {
    int     na, nb, nk;    /**< Model orders */
    int     n_hidden;       /**< Number of hidden neurons */
    int     n_input;        /**< Number of linear regressors (na + nb) */
    double *W1;             /**< Input-to-hidden weights [n_hidden x n_input] */
    double *b1;             /**< Hidden biases [n_hidden] */
    double *W2;             /**< Hidden-to-output weights [n_hidden] */
    double  b2;             /**< Output bias */
    double  Ts;             /**< Sampling period */
} sysid_sig_net_t;

/**
 * Initialize sigmoid network with random small weights (Xavier init).
 *
 * @param net       Network structure to initialize (must be freed by caller)
 * @param na        AR order
 * @param nb        Input order
 * @param nk        Input delay
 * @param n_hidden  Number of hidden neurons
 * @param Ts        Sampling period
 * @param seed      Random seed for weight initialization
 * @return          0 on success, -1 on error
 */
int sysid_sig_net_alloc(sysid_sig_net_t *net, int na, int nb, int nk,
                         int n_hidden, double Ts, unsigned int seed);

/** Free sigmoid network memory */
void sysid_sig_net_free(sysid_sig_net_t *net);

/**
 * Forward pass: compute f(phi) for a given regressor vector.
 *
 * @param net  Trained network
 * @param phi  Regressor vector [n_input]
 * @return     Network output f(phi)
 */
double sysid_sig_net_forward(const sysid_sig_net_t *net, const double *phi);

/**
 * Train sigmoid network using stochastic gradient descent with backpropagation.
 *
 * Minimizes sum_k (y(k) - f(phi(k)))^2 via online SGD.
 *
 * Each sample update:
 *   1. Forward pass: compute predicted output
 *   2. Backward pass: propagate error through output and hidden layers
 *   3. Update weights via gradient descent: w = w + lr * dL/dw
 *
 * @param net     Network to train (modified in-place)
 * @param u       Input data [N]
 * @param y       Target output data [N]
 * @param N       Number of samples
 * @param epochs  Number of training epochs over the data
 * @param lr      Learning rate (typically 0.001 to 0.1)
 * @return        0 on success, -1 on error
 */
int sysid_sig_net_train(sysid_sig_net_t *net, const double *u, const double *y,
                         int N, int epochs, double lr);

/**
 * Predict using trained sigmoid network NARX model.
 *
 * @param net   Trained network
 * @param u     Input data [N]
 * @param y     Output data [N] (past values for regressor construction)
 * @param N     Number of samples
 * @param ypred Output: predicted output [N]
 * @return      0 on success, -1 on error
 */
int sysid_sig_net_predict(const sysid_sig_net_t *net, const double *u,
                           const double *y, int N, double *ypred);

/*===========================================================================
 * L8: Hammerstein Model
 *
 * Structure: static nonlinearity -> linear dynamics
 *
 *   v(k) = g(u(k); alpha)         [static nonlinearity]
 *   y(k) = [B(q)/F(q)] v(k) + e(k) [linear dynamics]
 *
 * Identification (iterative Narendra-Gallman, 1966):
 *   1. Fix linear part, estimate nonlinearity via LS
 *   2. Fix nonlinearity, estimate linear part via PEM/LS
 *   3. Iterate until convergence
 *
 * The static nonlinearity is parameterized as a polynomial:
 *   g(u) = alpha_1 * u + alpha_2 * u^2 + ... + alpha_p * u^p
 *
 * This block structure is natural for systems where the sensor or actuator
 * has a static nonlinear characteristic (saturation, dead zone, etc.).
 *===========================================================================*/

/**
 * Estimate Hammerstein model parameters.
 *
 * @param u          Input data [N]
 * @param y          Output data [N]
 * @param N          Number of samples
 * @param nb         FIR length for linear dynamics
 * @param nk         Input delay
 * @param poly_order Polynomial order for static nonlinearity
 * @param alpha      Output: polynomial coefficients [poly_order]
 * @param b          Output: FIR coefficients [nb]
 * @param sigma2     Output: noise variance estimate (NULL to skip)
 * @param max_iter   Maximum iterations for Narendra-Gallman algorithm
 * @return           0 on success, -1 on error
 */
int sysid_hammerstein_estimate(const double *u, const double *y, int N,
                                int nb, int nk, int poly_order,
                                double *alpha, double *b,
                                double *sigma2, int max_iter);

/*===========================================================================
 * L8: Wiener Model
 *
 * Structure: linear dynamics -> static nonlinearity
 *
 *   w(k) = [B(q)/F(q)] u(k)        [linear dynamics]
 *   y(k) = h(w(k); beta) + e(k)    [static nonlinearity]
 *
 * This block structure models systems where the output measurement has
 * a static nonlinear characteristic (quantization, sensor saturation).
 *
 * Identification: iterative LS (identity NL init -> FIR -> output NL -> ...)
 *===========================================================================*/

/**
 * Estimate Wiener model parameters.
 *
 * @param u          Input data [N]
 * @param y          Output data [N]
 * @param N          Number of samples
 * @param nb         FIR length for linear dynamics
 * @param nk         Input delay
 * @param poly_order Polynomial order for output static nonlinearity
 * @param b          Output: FIR coefficients [nb]
 * @param beta       Output: polynomial coefficients [poly_order+1]
 *                    (includes constant term beta[0])
 * @param sigma2     Output: noise variance estimate (NULL to skip)
 * @return           0 on success, -1 on error
 */
int sysid_wiener_estimate(const double *u, const double *y, int N,
                           int nb, int nk, int poly_order,
                           double *b, double *beta, double *sigma2);

#endif /* SYSID_NONLINEAR_H */
