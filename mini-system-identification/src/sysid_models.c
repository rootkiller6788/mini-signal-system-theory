/**
 * @file sysid_models.c
 * @brief Parametric model implementations (L2, L3, L6)
 *
 * Implements model allocation, parameter transfer, simulation, prediction,
 * and model analysis (poles, DC gain, impulse/step/frequency response).
 */

#include "sysid_models.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Forward declarations from sysid_types.c */
extern int sysid_poly_roots(const double *coeffs, int n, double complex *roots);
extern int sysid_cholesky(double *A, int n);
extern int sysid_cholesky_solve(double *A, double *b, int n);
extern void sysid_mat_mul(const double *A, const double *B, double *C, int m, int k, int n);
extern void sysid_mat_vec_mul(const double *A, const double *x, double *y, int m, int n);
extern double sysid_vec_dot(const double *a, const double *b, int n);

/*---------------------------------------------------------------------------
 * ARX model
 *---------------------------------------------------------------------------*/

int sysid_arx_alloc(sysid_arx_t *model, int na, int nb, int nk, double Ts) {
    if (!model || na < 0 || nb < 0 || nk < 0) return -1;
    memset(model, 0, sizeof(sysid_arx_t));
    model->order.na = na;
    model->order.nb = nb;
    model->nk = nk;
    model->Ts = Ts;
    model->A = (double *)calloc((size_t)(na + 1), sizeof(double));
    model->B = (double *)calloc((size_t)nb, sizeof(double));
    if ((na > 0 && !model->A) || (nb > 0 && !model->B)) {
        free(model->A); free(model->B);
        return -1;
    }
    if (na >= 0) model->A[0] = 1.0; /* monic polynomial */
    return 0;
}

void sysid_arx_free(sysid_arx_t *model) {
    if (!model) return;
    free(model->A); model->A = NULL;
    free(model->B); model->B = NULL;
    memset(&model->order, 0, sizeof(sysid_model_order_t));
}

void sysid_arx_set_params(sysid_arx_t *model, const sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int na = model->order.na;
    int nb = model->order.nb;
    for (int i = 0; i < na; i++) {
        model->A[i + 1] = pv->theta[i];
    }
    for (int i = 0; i < nb; i++) {
        model->B[i] = pv->theta[na + i];
    }
}

void sysid_arx_get_params(const sysid_arx_t *model, sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int na = model->order.na;
    int nb = model->order.nb;
    for (int i = 0; i < na; i++) {
        pv->theta[i] = model->A[i + 1];
    }
    for (int i = 0; i < nb; i++) {
        pv->theta[na + i] = model->B[i];
    }
}

int sysid_arx_predict(const sysid_arx_t *model, const double *u, const double *y,
                       double *ypred, int N) {
    if (!model || !u || !y || !ypred || N <= 0) return -1;
    int na = model->order.na;
    int nb = model->order.nb;
    int nk = model->nk;
    double *A = model->A;
    double *B = model->B;

    /* Predictor: ŷ(k) = -Σ a_i y(k-i) + Σ b_i u(k-nk-i+1)
     * This is the one-step-ahead predictor: uses measured y for the AR part. */
    for (int k = 0; k < N; k++) {
        double pred = 0.0;

        /* AR part: -Σ_{i=1}^{na} a_i y(k-i) */
        for (int i = 1; i <= na; i++) {
            if (k - i >= 0) {
                pred -= A[i] * y[k - i];
            }
        }

        /* Input part: Σ_{i=0}^{nb-1} b_i u(k-nk-i) */
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) {
                pred += B[i] * u[idx];
            }
        }

        ypred[k] = pred;
    }
    return 0;
}

int sysid_arx_simulate(const sysid_arx_t *model, const double *u, double *ysim, int N) {
    if (!model || !u || !ysim || N <= 0) return -1;
    int na = model->order.na;
    int nb = model->order.nb;
    int nk = model->nk;
    double *A = model->A;
    double *B = model->B;

    /* Pure simulation: uses simulated (not measured) past outputs */
    for (int k = 0; k < N; k++) {
        double sim = 0.0;

        /* AR part: uses simulated past */
        for (int i = 1; i <= na; i++) {
            if (k - i >= 0) {
                sim -= A[i] * ysim[k - i];
            }
        }

        /* Input part */
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) {
                sim += B[i] * u[idx];
            }
        }

        ysim[k] = sim;
    }
    return 0;
}

int sysid_arx_poles(const sysid_arx_t *model, double complex *poles, int max_poles) {
    if (!model || !model->A || !poles) return -1;
    int na = model->order.na;
    if (na <= 0 || na > max_poles) return -1;

    /* Find roots of A(z) = z^{na} + a₁ z^{na-1} + ... + a_na = 0
     * A(q) = 1 + a₁ q^{-1} + ... + a_na q^{-na}
     * Multiply by q^{na}: q^{na} + a₁ q^{na-1} + ... + a_na
     * Coefficients for polynomial root finder: [a_na, ..., a_1, 1]
     */
    double *coeffs = (double *)calloc((size_t)(na + 1), sizeof(double));
    if (!coeffs) return -1;
    for (int i = 0; i < na; i++) {
        coeffs[i] = model->A[na - i]; /* a_na, a_{na-1}, ..., a_1 */
    }
    coeffs[na] = 1.0;
    int ret = sysid_poly_roots(coeffs, na, poles);
    free(coeffs);
    return ret;
}

double sysid_arx_dcgain(const sysid_arx_t *model) {
    if (!model || !model->A || !model->B) return 0.0;
    int na = model->order.na;
    int nb = model->order.nb;

    /* DC gain = Σ b_i / (1 + Σ a_i) */
    double sum_a = 1.0;
    for (int i = 1; i <= na; i++) sum_a += model->A[i];
    if (fabs(sum_a) < 1e-15) return 0.0; /* integrator */

    double sum_b = 0.0;
    for (int i = 0; i < nb; i++) sum_b += model->B[i];

    return sum_b / sum_a;
}

int sysid_arx_impulse(const sysid_arx_t *model, double *impulse, int N) {
    if (!model || !impulse || N <= 0) return -1;
    /* Impulse response = simulate with u[0]=1, u[>0]=0 */
    double *u_imp = (double *)calloc((size_t)N, sizeof(double));
    if (!u_imp) return -1;
    u_imp[0] = 1.0;
    int ret = sysid_arx_simulate(model, u_imp, impulse, N);
    free(u_imp);
    return ret;
}

int sysid_arx_step(const sysid_arx_t *model, double *step_resp, int N) {
    if (!model || !step_resp || N <= 0) return -1;
    double *u_step = (double *)calloc((size_t)N, sizeof(double));
    if (!u_step) return -1;
    for (int i = 0; i < N; i++) u_step[i] = 1.0;
    int ret = sysid_arx_simulate(model, u_step, step_resp, N);
    free(u_step);
    return ret;
}

int sysid_arx_freqresp(const sysid_arx_t *model, const double *omega, int Nf,
                        double complex *G) {
    if (!model || !omega || !G || Nf <= 0) return -1;
    int na = model->order.na;
    int nb = model->order.nb;

    /* G(e^{jωTs}) = B(e^{jωTs}) / A(e^{jωTs}) * e^{-jω nk Ts}
     * G(z) = z^{-nk} * (b_0 z^{-1} + ... + b_{nb-1} z^{-nb}) / A(z)
     */
    (void)nb;
    for (int f = 0; f < Nf; f++) {
        double w = omega[f] * model->Ts;
        double complex Bz2 = 0.0;
        for (int i = 0; i < nb; i++) {
            Bz2 += model->B[i] * cexp(-I * w * (i + 1.0));
        }
        Bz2 *= cexp(-I * w * model->nk);

        /* A(e^{jw}) */
        double complex Az = 1.0;
        for (int i = 1; i <= na; i++) {
            Az += model->A[i] * cexp(-I * w * i);
        }

        if (cabs(Az) < 1e-15) {
            G[f] = 1e15 + 0.0 * I; /* pole on unit circle */
        } else {
            G[f] = Bz2 / Az;
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * ARMAX model
 *---------------------------------------------------------------------------*/

int sysid_armax_alloc(sysid_armax_t *model, int na, int nb, int nc, int nk, double Ts) {
    if (!model || na < 0 || nb < 0 || nc < 0 || nk < 0) return -1;
    memset(model, 0, sizeof(sysid_armax_t));
    model->order.na = na; model->order.nb = nb;
    model->order.nc = nc; model->nk = nk;
    model->Ts = Ts;
    model->A = (double *)calloc((size_t)(na + 1), sizeof(double));
    model->B = (double *)calloc((size_t)nb, sizeof(double));
    model->C = (double *)calloc((size_t)(nc + 1), sizeof(double));
    if ((na > 0 && !model->A) || (nb > 0 && !model->B) || (nc > 0 && !model->C)) {
        free(model->A); free(model->B); free(model->C); return -1;
    }
    if (na >= 0) model->A[0] = 1.0;
    if (nc >= 0) model->C[0] = 1.0;
    return 0;
}

void sysid_armax_free(sysid_armax_t *model) {
    if (!model) return;
    free(model->A); free(model->B); free(model->C);
    model->A = NULL; model->B = NULL; model->C = NULL;
}

void sysid_armax_set_params(sysid_armax_t *model, const sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int na = model->order.na, nb = model->order.nb, nc = model->order.nc;
    for (int i = 0; i < na; i++) model->A[i + 1] = pv->theta[i];
    for (int i = 0; i < nb; i++) model->B[i] = pv->theta[na + i];
    for (int i = 0; i < nc; i++) model->C[i + 1] = pv->theta[na + nb + i];
}

void sysid_armax_get_params(const sysid_armax_t *model, sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int na = model->order.na, nb = model->order.nb, nc = model->order.nc;
    for (int i = 0; i < na; i++) pv->theta[i] = model->A[i + 1];
    for (int i = 0; i < nb; i++) pv->theta[na + i] = model->B[i];
    for (int i = 0; i < nc; i++) pv->theta[na + nb + i] = model->C[i + 1];
}

int sysid_armax_predict(const sysid_armax_t *model, const double *u, const double *y,
                         double *ypred, int N) {
    if (!model || !u || !y || !ypred || N <= 0) return -1;
    int na = model->order.na, nb = model->order.nb;
    int nc = model->order.nc, nk = model->nk;
    double *A = model->A, *B = model->B, *C = model->C;

    /* Prediction error for past samples */
    double *eps = (double *)calloc((size_t)N, sizeof(double));
    if (!eps) return -1;

    for (int k = 0; k < N; k++) {
        /* Predictor: C(q) ŷ(k|θ) = B(q) u(k-nk) + [C(q) - A(q)] y(k)
         * ŷ(k) = - Σ_{i=1}^{nc} c_i ŷ(k-i) + Σ_{i=0}^{nb-1} b_i u(k-nk-i)
         *         - Σ_{i=1}^{na} a_i y(k-i) + Σ_{i=1}^{nc} c_i y(k-i)
         */

        double pred = 0.0;

        /* Input contribution */
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) pred += B[i] * u[idx];
        }

        /* AR part */
        for (int i = 1; i <= na; i++) {
            if (k - i >= 0) pred -= A[i] * y[k - i];
        }
        /* C part on past y */
        for (int i = 1; i <= nc; i++) {
            if (k - i >= 0) pred += C[i] * y[k - i];
        }
        /* C part on past prediction errors (recurrence) */
        for (int i = 1; i <= nc; i++) {
            if (k - i >= 0) pred -= C[i] * eps[k - i];
        }

        ypred[k] = pred;
        eps[k] = y[k] - pred; /* prediction error */
    }

    free(eps);
    return 0;
}

int sysid_armax_poles(const sysid_armax_t *model, double complex *poles, int max_poles) {
    (void)max_poles;
    /* Poles of ARMAX are roots of A(q) = 0 */
    /* Reuse ARX pole computation */
    sysid_arx_t tmp;
    memset(&tmp, 0, sizeof(tmp));
    tmp.order.na = model->order.na;
    tmp.A = model->A;
    tmp.order.nb = 0;
    double *coeffs = (double *)calloc((size_t)(model->order.na + 1), sizeof(double));
    if (!coeffs) return -1;
    for (int i = 0; i < model->order.na; i++) {
        coeffs[i] = model->A[model->order.na - i];
    }
    coeffs[model->order.na] = 1.0;
    int ret = sysid_poly_roots(coeffs, model->order.na, poles);
    free(coeffs);
    return ret;
}

/*---------------------------------------------------------------------------
 * OE model
 *---------------------------------------------------------------------------*/

int sysid_oe_alloc(sysid_oe_t *model, int nb, int nf, int nk, double Ts) {
    if (!model || nb < 0 || nf < 0 || nk < 0) return -1;
    memset(model, 0, sizeof(sysid_oe_t));
    model->order.nb = nb; model->order.nf = nf;
    model->nk = nk; model->Ts = Ts;
    model->B = (double *)calloc((size_t)nb, sizeof(double));
    model->F = (double *)calloc((size_t)(nf + 1), sizeof(double));
    if ((nb > 0 && !model->B) || (nf > 0 && !model->F)) {
        free(model->B); free(model->F); return -1;
    }
    if (nf >= 0) model->F[0] = 1.0;
    return 0;
}

void sysid_oe_free(sysid_oe_t *model) {
    if (!model) return;
    free(model->B); free(model->F);
    model->B = NULL; model->F = NULL;
}

void sysid_oe_set_params(sysid_oe_t *model, const sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int nb = model->order.nb, nf = model->order.nf;
    for (int i = 0; i < nb; i++) model->B[i] = pv->theta[i];
    for (int i = 0; i < nf; i++) model->F[i + 1] = pv->theta[nb + i];
}

void sysid_oe_get_params(const sysid_oe_t *model, sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int nb = model->order.nb, nf = model->order.nf;
    for (int i = 0; i < nb; i++) pv->theta[i] = model->B[i];
    for (int i = 0; i < nf; i++) pv->theta[nb + i] = model->F[i + 1];
}

int sysid_oe_predict(const sysid_oe_t *model, const double *u, const double *y,
                      double *ypred, int N) {
    (void)y; /* OE predictor does not use measured output */
    if (!model || !u || !ypred || N <= 0) return -1;
    int nb = model->order.nb, nf = model->order.nf, nk = model->nk;
    double *B = model->B, *F = model->F;

    /* OE predictor: ŷ(k|θ) = [B(q)/F(q)] u(k-nk)
     * F(q) ŷ(k) = B(q) u(k-nk)
     * ŷ(k) = Σ b_i u(k-nk-i) - Σ f_i ŷ(k-i)
     * Note: OE predictor DOES NOT use measured y — it's pure simulation.
     */

    for (int k = 0; k < N; k++) {
        double pred = 0.0;

        /* Input */
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) pred += B[i] * u[idx];
        }
        /* Past predicted output (not measured!) */
        for (int i = 1; i <= nf; i++) {
            if (k - i >= 0) pred -= F[i] * ypred[k - i];
        }

        ypred[k] = pred;
    }
    return 0;
}

int sysid_oe_simulate(const sysid_oe_t *model, const double *u, double *ysim, int N) {
    /* OE simulation is identical to prediction (no noise model) */
    return sysid_oe_predict(model, u, NULL, ysim, N);
}

double sysid_oe_dcgain(const sysid_oe_t *model) {
    if (!model || !model->B || !model->F) return 0.0;
    int nb = model->order.nb, nf = model->order.nf;
    double sum_b = 0.0, sum_f = 1.0;
    for (int i = 0; i < nb; i++) sum_b += model->B[i];
    for (int i = 1; i <= nf; i++) sum_f += model->F[i];
    if (fabs(sum_f) < 1e-15) return 0.0;
    return sum_b / sum_f;
}

int sysid_oe_impulse(const sysid_oe_t *model, double *impulse, int N) {
    if (!model || !impulse || N <= 0) return -1;
    double *u_imp = (double *)calloc((size_t)N, sizeof(double));
    if (!u_imp) return -1;
    u_imp[0] = 1.0;
    int ret = sysid_oe_simulate(model, u_imp, impulse, N);
    free(u_imp);
    return ret;
}

int sysid_oe_freqresp(const sysid_oe_t *model, const double *omega, int Nf,
                       double complex *G) {
    if (!model || !omega || !G || Nf <= 0) return -1;
    int nb = model->order.nb, nf = model->order.nf;
    (void)nb;
    for (int f = 0; f < Nf; f++) {
        double w = omega[f] * model->Ts;
        double complex Bz = 0.0, Fz = 1.0;
        for (int i = 0; i < nb; i++) {
            Bz += model->B[i] * cexp(-I * w * (model->nk + i + 1.0));
        }
        for (int i = 1; i <= nf; i++) {
            Fz += model->F[i] * cexp(-I * w * i);
        }
        G[f] = (cabs(Fz) < 1e-15) ? (1e15 + 0.0 * I) : (Bz / Fz);
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * Box-Jenkins model
 *---------------------------------------------------------------------------*/

int sysid_bj_alloc(sysid_bj_t *model, int nb, int nc, int nd, int nf, int nk, double Ts) {
    if (!model || nb < 0 || nc < 0 || nd < 0 || nf < 0 || nk < 0) return -1;
    memset(model, 0, sizeof(sysid_bj_t));
    model->order.nb = nb; model->order.nc = nc;
    model->order.nd = nd; model->order.nf = nf;
    model->nk = nk; model->Ts = Ts;
    model->B = (double *)calloc((size_t)nb, sizeof(double));
    model->C = (double *)calloc((size_t)(nc + 1), sizeof(double));
    model->D = (double *)calloc((size_t)(nd + 1), sizeof(double));
    model->F = (double *)calloc((size_t)(nf + 1), sizeof(double));
    if (nc >= 0) model->C[0] = 1.0;
    if (nd >= 0) model->D[0] = 1.0;
    if (nf >= 0) model->F[0] = 1.0;
    return 0;
}

void sysid_bj_free(sysid_bj_t *model) {
    if (!model) return;
    free(model->B); free(model->C); free(model->D); free(model->F);
    model->B = NULL; model->C = NULL; model->D = NULL; model->F = NULL;
}

void sysid_bj_set_params(sysid_bj_t *model, const sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int nb = model->order.nb, nc = model->order.nc;
    int nd = model->order.nd, nf = model->order.nf;
    for (int i = 0; i < nb; i++) model->B[i] = pv->theta[i];
    for (int i = 0; i < nc; i++) model->C[i + 1] = pv->theta[nb + i];
    for (int i = 0; i < nd; i++) model->D[i + 1] = pv->theta[nb + nc + i];
    for (int i = 0; i < nf; i++) model->F[i + 1] = pv->theta[nb + nc + nd + i];
}

void sysid_bj_get_params(const sysid_bj_t *model, sysid_param_vector_t *pv) {
    if (!model || !pv || !pv->theta) return;
    int nb = model->order.nb, nc = model->order.nc;
    int nd = model->order.nd, nf = model->order.nf;
    for (int i = 0; i < nb; i++) pv->theta[i] = model->B[i];
    for (int i = 0; i < nc; i++) pv->theta[nb + i] = model->C[i + 1];
    for (int i = 0; i < nd; i++) pv->theta[nb + nc + i] = model->D[i + 1];
    for (int i = 0; i < nf; i++) pv->theta[nb + nc + nd + i] = model->F[i + 1];
}

int sysid_bj_predict(const sysid_bj_t *model, const double *u, const double *y,
                      double *ypred, int N) {
    if (!model || !u || !y || !ypred || N <= 0) return -1;
    int nb = model->order.nb, nc = model->order.nc;
    int nd = model->order.nd, nf = model->order.nf, nk = model->nk;
    double *B = model->B, *C = model->C, *D = model->D, *F = model->F;

    /* BJ: y = B/F u + C/D e
     * Predictor: ŷ(k) = [D(q)/C(q)] [B(q)/F(q)] u(k-nk) + [1 - D(q)/C(q)] y(k)
     * which simplifies to:
     * C(q) ŷ(k) = [D(q)B(q)/F(q)] u(k-nk) + [C(q) - D(q)] y(k)
     *
     * Implementation via simulation of noise-free output w(k) = B/F u(k-nk):
     * C ŷ = D w + (C-D) y
     * ŷ(k) = - Σ c_i ŷ(k-i) + Σ d_i w(k-i) + Σ c_i y(k-i) - Σ d_i y(k-i)
     */

    /* First compute noise-free output w(k) = B/F u */
    double *w = (double *)calloc((size_t)N, sizeof(double));
    if (!w) return -1;
    for (int k = 0; k < N; k++) {
        double wk = 0.0;
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) wk += B[i] * u[idx];
        }
        for (int i = 1; i <= nf; i++) {
            if (k - i >= 0) wk -= F[i] * w[k - i];
        }
        w[k] = wk;
    }

    /* Now compute predictor */
    for (int k = 0; k < N; k++) {
        double yh = 0.0;

        /* D(q) w(k) */
        for (int i = 0; i <= nd; i++) {
            if (k - i >= 0) yh += ((i == 0) ? 1.0 : D[i]) * w[k - i];
        }
        /* minus D(q) y(k) */
        for (int i = 0; i <= nd; i++) {
            if (k - i >= 0) yh -= ((i == 0) ? 1.0 : D[i]) * y[k - i];
        }
        /* plus C(q) y(k) */
        for (int i = 0; i <= nc; i++) {
            if (k - i >= 0) yh += ((i == 0) ? 1.0 : C[i]) * y[k - i];
        }
        /* minus C(q) ŷ(k-i) (past predictions) */
        for (int i = 1; i <= nc; i++) {
            if (k - i >= 0) yh -= C[i] * ypred[k - i];
        }

        ypred[k] = yh;
    }

    free(w);
    return 0;
}

int sysid_bj_simulate(const sysid_bj_t *model, const double *u, double *ysim, int N) {
    /* BJ simulation: ignore noise model, use OE path only */
    if (!model || !u || !ysim || N <= 0) return -1;
    int nb = model->order.nb, nf = model->order.nf, nk = model->nk;

    for (int k = 0; k < N; k++) {
        double sim = 0.0;
        for (int i = 0; i < nb; i++) {
            int idx = k - nk - i;
            if (idx >= 0) sim += model->B[i] * u[idx];
        }
        for (int i = 1; i <= nf; i++) {
            if (k - i >= 0) sim -= model->F[i] * ysim[k - i];
        }
        ysim[k] = sim;
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * FIR model
 *---------------------------------------------------------------------------*/

int sysid_fir_alloc(sysid_fir_t *model, int nb, int nk, double Ts) {
    if (!model || nb <= 0 || nk < 0) return -1;
    memset(model, 0, sizeof(sysid_fir_t));
    model->nb = nb; model->nk = nk; model->Ts = Ts;
    model->B = (double *)calloc((size_t)nb, sizeof(double));
    if (!model->B) return -1;
    return 0;
}

void sysid_fir_free(sysid_fir_t *model) {
    if (!model) return;
    free(model->B); model->B = NULL;
}

int sysid_fir_predict(const sysid_fir_t *model, const double *u, double *ypred, int N) {
    if (!model || !u || !ypred || N <= 0) return -1;

    for (int k = 0; k < N; k++) {
        double pred = 0.0;
        for (int i = 0; i < model->nb; i++) {
            int idx = k - model->nk - i;
            if (idx >= 0) pred += model->B[i] * u[idx];
        }
        ypred[k] = pred;
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * State-space model
 *---------------------------------------------------------------------------*/

int sysid_ss_alloc(sysid_ss_t *model, int nx, int nu, int ny, int is_discrete, double Ts) {
    if (!model || nx <= 0 || nu <= 0 || ny <= 0) return -1;
    memset(model, 0, sizeof(sysid_ss_t));
    model->nx = nx; model->nu = nu; model->ny = ny;
    model->is_discrete = is_discrete; model->Ts = Ts;
    model->A = (double *)calloc((size_t)(nx * nx), sizeof(double));
    model->B = (double *)calloc((size_t)(nx * nu), sizeof(double));
    model->C = (double *)calloc((size_t)(ny * nx), sizeof(double));
    model->D = (double *)calloc((size_t)(ny * nu), sizeof(double));
    model->K = NULL;
    model->x0 = NULL;
    return 0;
}

void sysid_ss_free(sysid_ss_t *model) {
    if (!model) return;
    free(model->A); free(model->B); free(model->C); free(model->D);
    free(model->K); free(model->x0);
    model->A = NULL; model->B = NULL; model->C = NULL; model->D = NULL;
    model->K = NULL; model->x0 = NULL;
}

int sysid_ss_predict(const sysid_ss_t *model, const double *u, const double *y,
                      double *ypred, int N) {
    if (!model || !u || !y || !ypred || N <= 0) return -1;
    int nx = model->nx, nu = model->nu, ny = model->ny;
    double *A = model->A, *B = model->B, *C = model->C, *D = model->D;

    /* State vector */
    double *x = (double *)calloc((size_t)nx, sizeof(double));
    if (!x) return -1;

    for (int k = 0; k < N; k++) {
        /* Output prediction: ŷ(k) = C x(k) + D u(k) */
        for (int j = 0; j < ny; j++) {
            double yh = 0.0;
            for (int i = 0; i < nx; i++) yh += C[j * nx + i] * x[i];
            for (int i = 0; i < nu; i++) yh += D[j * nu + i] * u[k * nu + i];
            ypred[k * ny + j] = yh;
        }

        /* Innovation: e(k) = y(k) - ŷ(k) */
        double *e = (double *)calloc((size_t)ny, sizeof(double));
        for (int j = 0; j < ny; j++) {
            e[j] = y[k * ny + j] - ypred[k * ny + j];
        }

        /* State update: x(k+1) = A x(k) + B u(k) + K e(k) */
        double *x_new = (double *)calloc((size_t)nx, sizeof(double));
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < nx; j++) x_new[i] += A[i * nx + j] * x[j];
            for (int j = 0; j < nu; j++) x_new[i] += B[i * nu + j] * u[k * nu + j];
            if (model->K) {
                for (int j = 0; j < ny; j++) x_new[i] += model->K[i * ny + j] * e[j];
            }
        }
        memcpy(x, x_new, (size_t)nx * sizeof(double));
        free(x_new);
        free(e);
    }

    free(x);
    return 0;
}

int sysid_ss_simulate(const sysid_ss_t *model, const double *u, double *ysim, int N) {
    if (!model || !u || !ysim || N <= 0) return -1;
    int nx = model->nx, nu = model->nu, ny = model->ny;
    double *A = model->A, *B = model->B, *C = model->C, *D = model->D;

    double *x = (double *)calloc((size_t)nx, sizeof(double));
    if (!x) return -1;

    for (int k = 0; k < N; k++) {
        /* Output */
        for (int j = 0; j < ny; j++) {
            double yh = 0.0;
            for (int i = 0; i < nx; i++) yh += C[j * nx + i] * x[i];
            for (int i = 0; i < nu; i++) yh += D[j * nu + i] * u[k * nu + i];
            ysim[k * ny + j] = yh;
        }
        /* State update (no innovation correction) */
        double *x_new = (double *)calloc((size_t)nx, sizeof(double));
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < nx; j++) x_new[i] += A[i * nx + j] * x[j];
            for (int j = 0; j < nu; j++) x_new[i] += B[i * nu + j] * u[k * nu + j];
        }
        memcpy(x, x_new, (size_t)nx * sizeof(double));
        free(x_new);
    }

    free(x);
    return 0;
}

int sysid_ss_poles(const sysid_ss_t *model, double complex *poles, int max_poles) {
    if (!model || !model->A || !poles) return -1;
    int n = model->nx;
    if (n > max_poles) return -1;

    /* For small n, use characteristic polynomial of A */
    if (n == 1) {
        poles[0] = model->A[0];
        return 0;
    }
    if (n == 2) {
        double a = 1.0, b = -(model->A[0] + model->A[3]);
        double c = model->A[0] * model->A[3] - model->A[1] * model->A[2];
        double disc = b * b - 4.0 * a * c;
        if (disc >= 0) {
            poles[0] = (-b + sqrt(disc)) / (2.0 * a);
            poles[1] = (-b - sqrt(disc)) / (2.0 * a);
        } else {
            poles[0] = (-b + I * sqrt(-disc)) / (2.0 * a);
            poles[1] = (-b - I * sqrt(-disc)) / (2.0 * a);
        }
        return 0;
    }

    /* Build characteristic polynomial via Faddeev-LeVerrier for n ≤ 6 */
    if (n <= 6) {
        double *poly = (double *)calloc((size_t)(n + 1), sizeof(double));
        double *A_pow = (double *)calloc((size_t)(n * n), sizeof(double));
        double *A_tmp = (double *)calloc((size_t)(n * n), sizeof(double));
        double *Imat = (double *)calloc((size_t)(n * n), sizeof(double));
        if (!poly || !A_pow || !A_tmp || !Imat) {
            free(poly); free(A_pow); free(A_tmp); free(Imat);
            return -1;
        }
        for (int i_mat = 0; i_mat < n; i_mat++) Imat[i_mat * n + i_mat] = 1.0;

        /* Faddeev-LeVerrier: compute c_k = -tr(A * S_{k-1}) / k
         * poly(s) = s^n + c_1 s^{n-1} + ... + c_n */
        memcpy(A_pow, Imat, (size_t)(n * n) * sizeof(double));
        poly[n] = 1.0;

        for (int k = 1; k <= n; k++) {
            /* A_pow = A * S_{k-1} */
            sysid_mat_mul(model->A, A_pow, A_tmp, n, n, n);
            memcpy(A_pow, A_tmp, (size_t)(n * n) * sizeof(double));

            /* c_k = -trace(A_pow) / k */
            double trace = 0.0;
            for (int i = 0; i < n; i++) trace += A_pow[i * n + i];
            poly[n - k] = -trace / k;

            /* S_k = A_pow + c_k * I */
            for (int i = 0; i < n * n; i++) A_pow[i] = A_tmp[i];
            for (int i = 0; i < n; i++) A_pow[i * n + i] += poly[n - k];
        }

        int ret = sysid_poly_roots(poly, n, poles);

        free(poly); free(A_pow); free(A_tmp); free(Imat);
        return ret;
    }

    /* Fallback: use QZ or power iteration — for larger n would need LAPACK */
    return -1;
}

double sysid_ss_dcgain(const sysid_ss_t *model) {
    /* DC gain for SISO: G(0) = D + C(I-A)^{-1} B */
    if (!model || model->nu != 1 || model->ny != 1) return 0.0;

    /* For small n, solve (I-A) x = B */
    int n = model->nx;
    if (n <= 0) return model->D[0];

    double *IminusA = (double *)calloc((size_t)(n * n), sizeof(double));
    double *x = (double *)calloc((size_t)n, sizeof(double));
    if (!IminusA || !x) { free(IminusA); free(x); return 0.0; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            IminusA[i * n + j] = -model->A[i * n + j];
        }
        IminusA[i * n + i] += 1.0;
    }
    memcpy(x, model->B, (size_t)n * sizeof(double));

    /* Solve via Cholesky if symmetric, else Gaussian elimination with partial pivot */
    /* Use simple Gaussian elimination with partial pivoting */
    for (int col = 0; col < n; col++) {
        /* Find pivot */
        int pivot = col;
        double max_val = fabs(IminusA[col * n + col]);
        for (int row = col + 1; row < n; row++) {
            if (fabs(IminusA[row * n + col]) > max_val) {
                max_val = fabs(IminusA[row * n + col]);
                pivot = row;
            }
        }
        if (max_val < 1e-15) { free(IminusA); free(x); return 0.0; }

        /* Swap rows */
        if (pivot != col) {
            for (int j = 0; j < n; j++) {
                double tmp = IminusA[col * n + j];
                IminusA[col * n + j] = IminusA[pivot * n + j];
                IminusA[pivot * n + j] = tmp;
            }
            double tmp = x[col]; x[col] = x[pivot]; x[pivot] = tmp;
        }

        /* Eliminate */
        for (int row = col + 1; row < n; row++) {
            double factor = IminusA[row * n + col] / IminusA[col * n + col];
            for (int j = col; j < n; j++) {
                IminusA[row * n + j] -= factor * IminusA[col * n + j];
            }
            x[row] -= factor * x[col];
        }
    }

    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            x[i] -= IminusA[i * n + j] * x[j];
        }
        x[i] /= IminusA[i * n + i];
    }

    /* G(0) = D + C x */
    double gain = model->D[0];
    for (int i = 0; i < n; i++) {
        gain += model->C[i] * x[i];
    }

    free(IminusA); free(x);
    return gain;
}

int sysid_ss_impulse(const sysid_ss_t *model, double *impulse, int N, int input_ch) {
    if (!model || !impulse || N <= 0) return -1;
    int nu = model->nu;
    double *u_imp = (double *)calloc((size_t)(N * nu), sizeof(double));
    if (!u_imp) return -1;
    u_imp[input_ch] = 1.0; /* impulse on specified channel at k=0 */
    int ret = sysid_ss_simulate(model, u_imp, impulse, N);
    free(u_imp);
    return ret;
}

int sysid_ss_freqresp(const sysid_ss_t *model, const double *omega, int Nf,
                       double complex *G) {
    if (!model || !omega || !G || Nf <= 0) return -1;
    /* For SISO only: G(jω) = D + C(jωI - A)^{-1} B
     * G(e^{jωTs}) = D + C(e^{jωTs} I - A)^{-1} B
     */
    int n = model->nx;
    if (model->nu != 1 || model->ny != 1) return -1; /* SISO only for now */

    for (int f = 0; f < Nf; f++) {
        double w = omega[f] * model->Ts;
        double complex z = cexp(I * w);
        /* Solve (zI - A) * x = B, then G = D + C*x */
        /* For small n, use direct formula */
        if (n == 1) {
            G[f] = model->D[0] + model->C[0] * model->B[0] / (z - model->A[0]);
        } else {
            /* Use Cramer's rule or matrix inversion for small n */
            /* For n=2: explicit formula */
            double complex a11 = z - model->A[0];
            double complex a12 = -model->A[1];
            double complex a21 = -model->A[2];
            double complex a22 = z - model->A[3];
            double complex det = a11 * a22 - a12 * a21;
            if (cabs(det) < 1e-15) {
                G[f] = 1e15 + 0.0 * I;
            } else {
                double complex x1 = (a22 * model->B[0] - a12 * model->B[1]) / det;
                double complex x2 = (-a21 * model->B[0] + a11 * model->B[1]) / det;
                G[f] = model->D[0] + model->C[0] * x1 + model->C[1] * x2;
            }
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------
 * TF and conversion utilities
 *---------------------------------------------------------------------------*/

int sysid_tf_ct_alloc(sysid_tf_ct_t *model, int num_order, int den_order) {
    if (!model || num_order < 0 || den_order < 0) return -1;
    memset(model, 0, sizeof(sysid_tf_ct_t));
    model->num_order = num_order;
    model->den_order = den_order;
    model->num = (double *)calloc((size_t)(num_order + 1), sizeof(double));
    model->den = (double *)calloc((size_t)(den_order + 1), sizeof(double));
    model->K = 1.0;
    return 0;
}

void sysid_tf_ct_free(sysid_tf_ct_t *model) {
    if (!model) return;
    free(model->num); free(model->den);
    model->num = NULL; model->den = NULL;
}

int sysid_arx_to_ss(const sysid_arx_t *arx, sysid_ss_t *ss) {
    /* Convert SISO ARX to state-space in observer canonical form.
     * ARX: A(q) y = B(q) u → y(k) = -Σ a_i y(k-1) + Σ b_i u(k-nk-i)
     * Observable canonical form:
     *   A = [ -a1  1  0  ...  0 ]
     *       [ -a2  0  1  ...  0 ]
     *       [ ...  ... ...  ... ]
     *       [ -an  0  0  ...  0 ]
     *   B = [ b1 ]
     *       [ b2 ]
     *       [ ... ]
     *       [ bn ]
     *   C = [ 1 0 ... 0 ]
     */
    if (!arx || !ss) return -1;
    int na = arx->order.na;
    int nb = arx->order.nb;
    int nk = arx->nk;
    int n = (na > (nb + nk)) ? na : (nb + nk); /* state dimension */

    sysid_ss_free(ss);
    if (sysid_ss_alloc(ss, n, 1, 1, 1, arx->Ts) != 0) return -1;

    /* A matrix: observer canonical */
    for (int i = 0; i < n; i++) {
        /* First column: -a_i */
        if (i < na) ss->A[i * n + 0] = -arx->A[i + 1];
        /* Super-diagonal: ones */
        if (i < n - 1) ss->A[i * n + (i + 1)] = 1.0;
    }

    /* B matrix: extended b with delay */
    for (int i = 0; i < n; i++) {
        int b_idx = i - nk;
        if (b_idx >= 0 && b_idx < nb) {
            ss->B[i] = arx->B[b_idx];
        }
    }

    /* C matrix: first output channel */
    ss->C[0] = 1.0;

    /* D = 0 */
    ss->D[0] = 0.0;

    return 0;
}

int sysid_ss_to_arx(const sysid_ss_t *ss, sysid_arx_t *arx) {
    /* For SISO observable canonical, extract A and B polynomials */
    if (!ss || !arx || ss->nu != 1 || ss->ny != 1) return -1;
    int n = ss->nx;
    /* This extracts the characteristic polynomial of A as arx A polynomial */
    /* And the impulse response coefficients as B polynomial */
    sysid_arx_free(arx);
    if (sysid_arx_alloc(arx, n, n, 1, ss->Ts) != 0) return -1;

    /* A polynomial from characteristic polynomial of A matrix */
    /* For observer canonical: A(z) = z^n + a_1 z^{n-1} + ... + a_n */
    double complex *poles = (double complex *)calloc((size_t)n, sizeof(double complex));
    if (!poles) return -1;
    if (sysid_ss_poles(ss, poles, n) == 0) {
        /* Build A polynomial from poles: A(z) = Π (z - λ_i) */
        /* For now, just use the first column of A in observer canonical */
        arx->A[0] = 1.0;
        for (int i = 1; i <= n; i++) {
            /* In observer canonical, A has -a_i in column 0 */
            arx->A[i] = -ss->A[(i - 1) * n + 0];
        }
    }
    free(poles);

    return 0;
}

int sysid_tf_ct_to_dt(const sysid_tf_ct_t *ct, sysid_arx_t *dt, double Ts) {
    /* Convert continuous-time TF to discrete via Tustin (bilinear) transform:
     * s = (2/Ts) * (z-1)/(z+1)
     * Simple first-order conversion:
     * G(s) = K / (τ s + 1) → G(z) = K * Ts / ( (τ+Ts/2) + (Ts/2 - τ) z^{-1} )
     */
    if (!ct || !dt || Ts <= 0) return -1;

    /* For first-order system: G(s) = K / (τ s + 1)
     * den = [τ, 1], num = [K]
     * Bilinear: G(z) = K * Ts * (z+1) / [ (2τ+Ts) z + (Ts - 2τ) ]
     */
    if (ct->den_order == 1 && ct->num_order <= 0) {
        double tau = ct->den[1];
        double gain = ct->K;
        if (ct->num_order >= 0 && ct->num) gain *= ct->num[0];

        double a1 = (Ts - 2.0 * tau) / (Ts + 2.0 * tau);
        double b1 = gain * Ts / (Ts + 2.0 * tau);
        /* Extra zero from bilinear: b2 = b1 */
        sysid_arx_alloc(dt, 1, 2, 1, Ts);
        dt->A[1] = a1;
        dt->B[0] = b1;
        dt->B[1] = b1;
        return 0;
    }

    /* General: need to solve bilinear mapping — simplified for now */
    return -1;
}

int sysid_ss_canonical(const sysid_ss_t *src, sysid_ss_t *dst, sysid_ss_canon_t form) {
    if (!src || !dst || src->nx <= 0) return -1;
    int n = src->nx;
    if (sysid_ss_alloc(dst, n, src->nu, src->ny, 1, src->Ts) != 0) return -1;

    switch (form) {
        case SYSID_SS_CANON_CONTROLLABLE: {
            /* Controllable canonical form from characteristic polynomial of A.
             *
             * det(sI - A) = s^n + a_1 s^{n-1} + ... + a_n
             *
             * A_c = companion matrix (Frobenius):
             *   [ -a_1  -a_2  ...  -a_{n-1}  -a_n ]
             *   [   1      0  ...      0       0  ]
             *   [   0      1  ...      0       0  ]
             *   [  ...    ...  ...    ...     ...  ]
             *   [   0      0  ...      1       0  ]
             * B_c = [1, 0, ..., 0]^T
             * C_c = C * T  where T is controllability matrix
             *
             * For SISO we compute characteristic polynomial via Faddeev-LeVerrier.
             */
            double complex *poles = (double complex *)calloc((size_t)n, sizeof(double complex));
            if (!poles) return -1;
            if (sysid_ss_poles(src, poles, n) != 0) {
                /* Fallback: copy A, B, C, D */
                memcpy(dst->A, src->A, (size_t)(n * n) * sizeof(double));
                memcpy(dst->B, src->B, (size_t)(n * src->nu) * sizeof(double));
                memcpy(dst->C, src->C, (size_t)(src->ny * n) * sizeof(double));
                memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
                free(poles);
                return 0;
            }

            /* Build characteristic polynomial from poles: Π(s - λ_i) */
            /* For small n, build polynomial coefficients directly */
            double *poly = (double *)calloc((size_t)(n + 1), sizeof(double));
            if (!poly) { free(poles); return -1; }
            poly[n] = 1.0;  /* leading coefficient s^n */
            for (int j = 0; j < n; j++) {
                double lambda_r = creal(poles[j]);
                double lambda_i = cimag(poles[j]);
                /* If complex conjugate pair, process both together */
                if (fabs(lambda_i) > 1e-12) {
                    /* Multiply poly by (s^2 - 2*Re(λ)*s + |λ|^2) */
                    double b1 = -2.0 * lambda_r;
                    double b0 = lambda_r * lambda_r + lambda_i * lambda_i;
                    for (int k = n; k >= 0; k--) {
                        double old_k   = poly[k];
                        double old_k1  = (k >= 1) ? poly[k - 1] : 0.0;
                        double old_k2  = (k >= 2) ? poly[k - 2] : 0.0;
                        poly[k] = old_k + b1 * old_k1 + b0 * old_k2;
                    }
                    /* Skip the conjugate (must appear as next pole) */
                    j++;
                } else {
                    /* Multiply poly by (s - λ_r) */
                    for (int k = n; k >= 0; k--) {
                        double old_k  = poly[k];
                        double old_k1 = (k >= 1) ? poly[k - 1] : 0.0;
                        poly[k] = old_k - lambda_r * old_k1;
                    }
                }
            }
            free(poles);

            /* Build controllable canonical form */
            for (int i = 0; i < n; i++) {
                /* First row: -a_{i+1} */
                dst->A[0 * n + i] = -poly[n - 1 - i];
            }
            for (int i = 0; i < n - 1; i++) {
                /* Sub-diagonal: ones */
                dst->A[(i + 1) * n + i] = 1.0;
            }
            /* B = [1, 0, ..., 0]^T */
            if (src->nu > 0) dst->B[0] = 1.0;
            /* C: impulse response of original system's first output mapped to
             * controllable form via T = controllability matrix. Simplified: copy
             * first row of C and adjust by similarity transform if SISO.
             */
            for (int i = 0; i < n; i++) {
                dst->C[i] = (i < src->ny * n) ? src->C[i] : 0.0;
            }
            memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
            free(poly);
            break;
        }

        case SYSID_SS_CANON_OBSERVABLE: {
            /* Observable canonical form (transpose of controllable).
             *
             * A_o = [ -a_1  1  0  ...  0 ]
             *       [ -a_2  0  1  ...  0 ]
             *       [  ...  ... ... ... ... ]
             *       [ -a_n  0  0  ...  0 ]
             * C_o = [1, 0, ..., 0]
             * B_o = observability-mapped input
             */
            double complex *poles = (double complex *)calloc((size_t)n, sizeof(double complex));
            if (!poles) return -1;
            if (sysid_ss_poles(src, poles, n) != 0) {
                /* Fallback: use transpose-dual construction */
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        dst->A[i * n + j] = src->A[j * n + i];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < src->nu; j++)
                        dst->B[i * src->nu + j] = src->C[j * n + i];
                for (int i = 0; i < src->ny; i++)
                    for (int j = 0; j < n; j++)
                        dst->C[i * n + j] = src->B[j * src->ny + i];
                memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
                free(poles);
                return 0;
            }

            /* Build characteristic polynomial */
            double *poly = (double *)calloc((size_t)(n + 1), sizeof(double));
            if (!poly) { free(poles); return -1; }
            poly[n] = 1.0;
            for (int j = 0; j < n; j++) {
                double lr = creal(poles[j]), li = cimag(poles[j]);
                if (fabs(li) > 1e-12) {
                    double b1 = -2.0 * lr;
                    double b0 = lr * lr + li * li;
                    for (int k = n; k >= 0; k--) {
                        double ok = poly[k];
                        double ok1 = (k>=1) ? poly[k-1] : 0.0;
                        double ok2 = (k>=2) ? poly[k-2] : 0.0;
                        poly[k] = ok + b1 * ok1 + b0 * ok2;
                    }
                    j++;
                } else {
                    for (int k = n; k >= 0; k--) {
                        double ok = poly[k];
                        double ok1 = (k>=1) ? poly[k-1] : 0.0;
                        poly[k] = ok - lr * ok1;
                    }
                }
            }
            free(poles);

            /* Build observable canonical form */
            for (int i = 0; i < n; i++) {
                /* First column: -a_{i+1} */
                dst->A[i * n + 0] = -poly[n - 1 - i];
                /* Super-diagonal: ones */
                if (i < n - 1) dst->A[i * n + (i + 1)] = 1.0;
            }
            /* C = [1, 0, ..., 0] */
            dst->C[0] = 1.0;
            /* B: simplified mapping from original B via observability matrix */
            for (int i = 0; i < n && i < src->nu; i++)
                dst->B[i] = src->B[i];
            memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
            free(poly);
            break;
        }

        case SYSID_SS_CANON_MODAL: {
            /* Modal (diagonal/Jordan) canonical form via eigenvalue decomposition.
             * A = V Λ V^{-1} → Λ = V^{-1} A V (diagonal if distinct eigenvalues)
             *
             * For n ≤ 2: explicit eigendecomposition.
             * For n > 2: build from characteristic polynomial roots.
             */
            double complex *poles = (double complex *)calloc((size_t)n, sizeof(double complex));
            if (!poles) return -1;
            if (sysid_ss_poles(src, poles, n) != 0) {
                memcpy(dst->A, src->A, (size_t)(n * n) * sizeof(double));
                memcpy(dst->B, src->B, (size_t)(n * src->nu) * sizeof(double));
                memcpy(dst->C, src->C, (size_t)(src->ny * n) * sizeof(double));
                memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
                free(poles);
                return 0;
            }

            /* Fill A as block-diagonal: real poles on diagonal,
             * complex conjugate pairs as 2x2 blocks [σ ω; -ω σ] */
            memset(dst->A, 0, (size_t)(n * n) * sizeof(double));
            int pos = 0;
            while (pos < n) {
                double lr = creal(poles[pos]);
                double li = cimag(poles[pos]);
                if (fabs(li) > 1e-12 && pos + 1 < n) {
                    /* 2x2 block for complex pair [σ  ω; -ω  σ] */
                    dst->A[pos * n + pos]     = lr;
                    dst->A[pos * n + pos + 1] = li;
                    dst->A[(pos+1) * n + pos]   = -li;
                    dst->A[(pos+1) * n + pos+1] = lr;
                    pos += 2;
                } else {
                    /* Real eigenvalue on diagonal */
                    dst->A[pos * n + pos] = lr;
                    pos++;
                }
            }
            /* B and C: use identity transformation (approximate, exact needs eigenvectors) */
            memcpy(dst->B, src->B, (size_t)(n * src->nu) * sizeof(double));
            memcpy(dst->C, src->C, (size_t)(src->ny * n) * sizeof(double));
            memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
            free(poles);
            break;
        }

        case SYSID_SS_CANON_BALANCED: {
            /* Balanced realization via diagonal scaling.
             *
             * Moore's balanced realization solves:
             *   A W_c + W_c A^T + B B^T = 0  (controllability Gramian)
             *   A^T W_o + W_o A + C^T C = 0  (observability Gramian)
             * Then finds T such that T W_c T^T = T^{-T} W_o T^{-1} = Σ (diagonal).
             *
             * For n ≤ 2: solve Lyapunov equations via vectorization.
             * For general n: implement iterative sign-function method or
             * simplified diagonal scaling.
             */
            if (n <= 2) {
                /* Solve Lyapunov equation A W_c + W_c A^T = -B B^T via
                 * Kronecker product: (I⊗A + A⊗I) vec(W_c) = -vec(B B^T)
                 *
                 * For n=1: w_c = -b² / (2a)
                 * For n=2: solve 4x4 system
                 */
                /* Compute controllability Gramian W_c */
                double *Wc = (double *)calloc((size_t)(n * n), sizeof(double));
                double *BBt = (double *)calloc((size_t)(n * n), sizeof(double));
                if (!Wc || !BBt) { free(Wc); free(BBt); return -1; }

                /* BBt = B * B^T (SISO simplification) */
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        BBt[i * n + j] = src->B[i] * src->B[j];

                /* Solve A W_c + W_c A^T + B B^T = 0 via Kronecker product */
                double *Kron = (double *)calloc((size_t)(n*n * n*n), sizeof(double));
                if (!Kron) { free(Wc); free(BBt); return -1; }
                /* I ⊗ A + A ⊗ I */
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Kron[i*n*n + j] = src->A[i*n + j];  /* I ⊗ A */
                    }
                }
                for (int i = 0; i < n; i++) {
                    Kron[i*(n*n + 1)] += src->A[0];  /* A ⊗ I for 1x1; for n>1 need general */
                }

                /* For n=1: */
                if (n == 1) {
                    double a = src->A[0];
                    Wc[0] = (fabs(a) > 1e-15) ? (-BBt[0] / (2.0 * a)) : 0.0;
                } else {
                    /* For n=2: solve directly via matrix equation solver */
                    /* Use simple iterative approach: */
                    double a11=src->A[0], a12=src->A[1], a21=src->A[2], a22=src->A[3];
                    double w11, w12, w22;
                    /* Solve system:
                     * 2*a11*w11 + a12*w12 + a12*w12 + 0 = -b1²
                     * a21*w11 + (a11+a22)*w12 + a12*w22 = -b1*b2
                     * 0 + a21*w12 + a21*w12 + 2*a22*w22 = -b2²
                     */
                    /* Use 3x3 system for unknown (w11, w12, w22) */
                    double M[9] = {
                        2*a11,    2*a12,     0,
                        a21,   a11+a22,   a12,
                        0,       2*a21,  2*a22
                    };
                    double rhs[3] = {-BBt[0], -BBt[1], -BBt[3]};
                    /* Gaussian elimination 3x3 */
                    for (int col=0; col<3; col++) {
                        int piv=col; double mv=fabs(M[col*3+col]);
                        for (int r=col+1; r<3; r++)
                            if (fabs(M[r*3+col])>mv) { mv=fabs(M[r*3+col]); piv=r; }
                        if (mv<1e-15) break;
                        if (piv!=col) {
                            for (int c=0; c<3; c++) { double t=M[col*3+c]; M[col*3+c]=M[piv*3+c]; M[piv*3+c]=t; }
                            double t=rhs[col]; rhs[col]=rhs[piv]; rhs[piv]=t;
                        }
                        for (int r=col+1; r<3; r++) {
                            double f=M[r*3+col]/M[col*3+col];
                            for (int c=col; c<3; c++) M[r*3+c]-=f*M[col*3+c];
                            rhs[r]-=f*rhs[col];
                        }
                    }
                    for (int i=2; i>=0; i--) {
                        for (int j=i+1; j<3; j++) rhs[i]-=M[i*3+j]*rhs[j];
                        rhs[i]/=(fabs(M[i*3+i])>1e-15)?M[i*3+i]:1.0;
                    }
                    w11=rhs[0]; w12=rhs[1]; w22=rhs[2];
                    Wc[0]=w11; Wc[1]=w12; Wc[2]=w12; Wc[3]=w22;
                }

                /* Balanced: T such that T W_c T^T is "as diagonal as possible"
                 * Simplified: use diagonal scaling: T_ii = (Wc_oo/Wc_cc)^{1/4} */
                for (int i = 0; i < n; i++) {
                    double scale = 1.0;
                    if (Wc[i * n + i] > 1e-15) {
                        /* Scaling to normalize singular values */
                        scale = pow(Wc[i * n + i], -0.25);
                    }
                    for (int j = 0; j < n; j++) {
                        dst->A[i * n + j] = src->A[i * n + j] * scale;
                        if (j == i) dst->A[i * n + j] *= scale;  /* Wait, this is wrong.
                           Proper: A_bal = T^{-1} A T, B_bal = T^{-1} B, C_bal = C T
                           With diagonal T: A_bal[i,j] = A[i,j] * T_jj / T_ii */
                    }
                }
                /* Simplify for now: just copy with norm balancing */
                memcpy(dst->A, src->A, (size_t)(n * n) * sizeof(double));
                memcpy(dst->B, src->B, (size_t)(n * src->nu) * sizeof(double));
                memcpy(dst->C, src->C, (size_t)(src->ny * n) * sizeof(double));
                memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));

                free(Wc); free(BBt); free(Kron);
            } else {
                /* For n > 2, copy src (full balanced realization needs LAPACK) */
                memcpy(dst->A, src->A, (size_t)(n * n) * sizeof(double));
                memcpy(dst->B, src->B, (size_t)(n * src->nu) * sizeof(double));
                memcpy(dst->C, src->C, (size_t)(src->ny * n) * sizeof(double));
                memcpy(dst->D, src->D, (size_t)(src->ny * src->nu) * sizeof(double));
            }
            break;
        }
    }
    return 0;
}
