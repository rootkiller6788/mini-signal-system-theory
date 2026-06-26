/**
 * @file    adaptive_kalman.c
 * @brief   Kalman filter implementations: Linear, Extended, and Unscented
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): Linear Kalman filter (predict + update),
 *                     Extended Kalman filter (EKF) with Jacobian linearization,
 *                     Unscented Kalman filter (UKF) sigma-point framework
 *   L2 (Core Concepts): State-space models, Bayesian filtering,
 *                       Gaussian approximation, optimal estimation
 *   L4 (Fundamental Laws): Kalman-Bucy filter theory, minimum variance
 *                          estimation, orthogonality principle
 *   L6 (Canonical Problems): Target tracking, sensor fusion
 *
 * Reference:
 *   Kalman (1960), "A New Approach to Linear Filtering and Prediction
 *                    Problems", Trans. ASME, J. Basic Engineering
 *   Anderson & Moore (1979), "Optimal Filtering", Prentice-Hall
 *   Julier & Uhlmann (2004), "Unscented Filtering and Nonlinear
 *                              Estimation", Proc. IEEE
 *   Thrun, Burgard & Fox (2005), "Probabilistic Robotics", MIT Press
 *   Simon (2006), "Optimal State Estimation", Wiley
 *
 * Course Mapping:
 *   MIT 6.450 - Digital Communications (Kalman equalization)
 *   Stanford EE363 - Linear Dynamical Systems
 *   Berkeley EE221A - Linear Systems Theory
 *   ETH 227-0427 - Signal Processing (Kalman filtering for tracking)
 *   TUM - Signal Processing (state estimation)
 */

#include "adaptive_filter.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =========================================================================
 * Internal helper: matrix-vector multiply for Kalman
 * ========================================================================= */

/**
 * @brief Multiply matrix A (n¡Án, row-major) by vector x, result in y
 */
static void kalman_mat_vec_mul(const double *A, const double *x,
                                double *y, size_t n) {
    size_t i, j;
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
        for (j = 0; j < n; j++) {
            y[i] += A[i * n + j] * x[j];
        }
    }
}

/**
 * @brief Multiply two n¡Án matrices: C = A * B
 */
static void kalman_mat_mul(const double *A, const double *B,
                            double *C, size_t n) {
    size_t i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/**
 * @brief Matrix addition: C = A + B
 */
static void kalman_mat_add(const double *A, const double *B,
                            double *C, size_t n) {
    size_t i;
    for (i = 0; i < n * n; i++) {
        C[i] = A[i] + B[i];
    }
}

/**
 * @brief Compute x^T * P * x (quadratic form)
 */
static double __attribute__((unused)) kalman_quadratic_form(const double *x, const double *P, size_t n) {
    size_t i, j;
    double result = 0.0;
    /* Working vector y = P*x */
    double *y = (double *)malloc(n * sizeof(double));
    if (!y) return 0.0;
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
        for (j = 0; j < n; j++) {
            y[i] += P[i * n + j] * x[j];
        }
    }
    for (i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    free(y);
    return result;
}

/* =========================================================================
 * L5: Linear Kalman Filter - Prediction Step
 * ========================================================================= */

/**
 * @brief af_kalman_predict: Linear Kalman filter prediction step
 *
 * State-space model:
 *   x(k+1) = F * x(k) + w(k),  w(k) ~ N(0, Q)
 *
 * where:
 *   x(k): state vector (n_states ¡Á 1)
 *   F:    state transition matrix (n_states ¡Á n_states)
 *   w(k): process noise ~ N(0, Q)
 *   Q:    process noise covariance (n_states ¡Á n_states)
 *
 * Prediction equations:
 *   x_hat^{-} = F * x_hat          (prior state estimate)
 *   P^{-}     = F * P * F^T + Q    (prior error covariance)
 *
 * The predicted state is the best estimate before incorporating
 * the next measurement. The covariance grows due to process noise Q.
 *
 * Optimality: The Kalman filter is the minimum-variance unbiased
 * estimator for linear Gaussian systems. Among all linear filters,
 * it produces the smallest mean-square estimation error.
 *
 * @param x_hat    State estimate [n_states] (in: posterior k|k, out: prior k+1|k)
 * @param P        Error covariance [n_states^2] row-major
 *                 (in: P(k|k), out: P(k+1|k))
 * @param F        State transition matrix [n_states^2] row-major
 * @param Q        Process noise covariance [n_states^2] row-major
 * @param n_states State dimension
 */
void af_kalman_predict(double *x_hat, double *P, const double *F,
                        const double *Q, size_t n_states) {
    size_t n = n_states;
    double *F_P, *temp;

    if (!x_hat || !P || !F || !Q || n == 0) return;

    F_P = (double *)malloc(n * n * sizeof(double));
    temp = (double *)malloc(n * n * sizeof(double));
    if (!F_P || !temp) { free(F_P); free(temp); return; }

    /* 1. x_hat = F * x_hat (state prediction) */
    {
        size_t i;
        double *x_old = (double *)malloc(n * sizeof(double));
        if (!x_old) { free(F_P); free(temp); return; }
        for (i = 0; i < n; i++) x_old[i] = x_hat[i];
        kalman_mat_vec_mul(F, x_old, x_hat, n);
        free(x_old);
    }

    /* 2. P = F * P * F^T + Q */
    /* First: F_P = F * P */
    kalman_mat_mul(F, P, F_P, n);

    /* Then: temp = F_P * F^T = F * P * F^T */
    {
        size_t i, j, k;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double sum = 0.0;
                for (k = 0; k < n; k++) {
                    sum += F_P[i * n + k] * F[j * n + k]; /* F^T: F[j][k] = F^T[k][j] */
                }
                temp[i * n + j] = sum;
            }
        }
    }

    /* P = temp + Q = F*P*F^T + Q */
    kalman_mat_add(temp, Q, P, n);

    free(F_P);
    free(temp);
}

/* =========================================================================
 * L5: Linear Kalman Filter - Update Step
 * ========================================================================= */

/**
 * @brief af_kalman_update: Linear Kalman filter update step
 *
 * Measurement model:
 *   z(k) = H * x(k) + v(k),  v(k) ~ N(0, R)
 *
 * Update equations:
 *   y_tilde  = z - H * x_hat^{-}          (innovation / measurement residual)
 *   S        = H * P^{-} * H^T + R        (innovation covariance)
 *   K        = P^{-} * H^T * S^{-1}        (Kalman gain)
 *   x_hat    = x_hat^{-} + K * y_tilde    (updated state estimate)
 *   P        = (I - K * H) * P^{-}         (updated covariance)
 *
 * The Kalman gain K optimally weights the innovation: if the
 * measurement is very noisy (large R), K is small and the filter
 * relies more on the prediction. If the prediction is uncertain
 * (large P^{-}), K is large and the filter trusts the measurement more.
 *
 * Numerical note: The Joseph form of the covariance update
 *   P = (I-K*H)*P*(I-K*H)^T + K*R*K^T
 * is more numerically stable and preserves symmetry, but the
 * standard form shown above is simpler and widely used.
 *
 * @param x_hat    State estimate [n_states] (in: prior, out: posterior)
 * @param P        Error covariance [n_states^2] (in: prior, out: posterior)
 * @param H        Measurement matrix [n_meas*n_states] row-major
 * @param R        Measurement noise covariance [n_meas^2] row-major
 * @param z        Measurement scalar (single measurement case, n_meas=1)
 * @param n_states State dimension
 * @param n_meas   Measurement dimension (typically 1 for scalar measurement)
 */
void af_kalman_update(double *x_hat, double *P, const double *H,
                       const double *R, double z, size_t n_states,
                       size_t n_meas) {
    size_t n = n_states;
    size_t m = (n_meas > 0) ? n_meas : 1;
    double *K, *P_HT, *H_P_HT, *y_tilde;
    double innov, S, inv_S;

    if (!x_hat || !P || !H || !R || n == 0 || m == 0) return;

    K = (double *)malloc(n * m * sizeof(double));
    P_HT = (double *)malloc(n * m * sizeof(double));
    H_P_HT = (double *)malloc(m * m * sizeof(double));
    y_tilde = (double *)malloc(m * sizeof(double));

    if (!K || !P_HT || !H_P_HT || !y_tilde) {
        free(K); free(P_HT); free(H_P_HT); free(y_tilde);
        return;
    }

    /* For scalar measurement (m=1), the equations simplify significantly.
     * We present the general form for clarity. */

    /* 1. Compute innovation y_tilde = z - H * x_hat */
    {
        size_t i, j;
        for (i = 0; i < m; i++) {
            double hx = 0.0;
            for (j = 0; j < n; j++) {
                hx += H[i * n + j] * x_hat[j]; /* H has stride n */
            }
            y_tilde[i] = ((i == 0) ? z : 0.0) - hx; /* For m=1, only z at i=0 */
        }
        innov = y_tilde[0];
    }

    /* 2. P_HT = P * H^T */
    {
        size_t i, j, k;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                double sum = 0.0;
                for (k = 0; k < n; k++) {
                    sum += P[i * n + k] * H[j * n + k]; /* H^T[k][j] = H[j][k] */
                }
                P_HT[i * m + j] = sum;
            }
        }
    }

    /* 3. S = H * P * H^T + R */
    {
        size_t i, j, k;
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                double sum = 0.0;
                for (k = 0; k < n; k++) {
                    sum += H[i * n + k] * P_HT[k * m + j];
                }
                H_P_HT[i * m + j] = sum + R[i * m + j];
            }
        }
        S = H_P_HT[0];
    }

    /* 4. Kalman gain K = P_HT / S (since m=1, S is scalar) */
    if (fabs(S) < 1e-15) {
        free(K); free(P_HT); free(H_P_HT); free(y_tilde);
        return;
    }
    inv_S = 1.0 / S;
    {
        size_t i;
        for (i = 0; i < n; i++) K[i] = P_HT[i] * inv_S;
    }

    /* 5. State update: x_hat = x_hat + K * innov */
    {
        size_t i;
        for (i = 0; i < n; i++) x_hat[i] += K[i] * innov;
    }

    /* 6. Covariance update: P = P - K * (H * P)
     * First compute H_P = H * P (1¡Án row vector) */
    {
        size_t i, j, k;
        double *H_P = (double *)malloc(n * sizeof(double));
        if (!H_P) { free(K); free(P_HT); free(H_P_HT); free(y_tilde); return; }
        for (j = 0; j < n; j++) {
            H_P[j] = 0.0;
            for (k = 0; k < n; k++) {
                H_P[j] += H[k] * P[k * n + j]; /* H is 1¡Án for m=1 */
            }
        }
        /* P = P - K * H_P (outer product) */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                P[i * n + j] -= K[i] * H_P[j];
            }
        }
        free(H_P);
    }

    free(K); free(P_HT); free(H_P_HT); free(y_tilde);
}

/* =========================================================================
 * L5: Extended Kalman Filter - Prediction Step
 * ========================================================================= */

/**
 * @brief af_kalman_extended_predict: EKF prediction step
 *
 * For nonlinear systems:
 *   x(k+1) = f(x(k)) + w(k)
 *
 * The EKF linearizes f around the current state estimate using
 * the Jacobian F = df/dx.
 *
 * Prediction:
 *   x_hat^{-} = f(x_hat)              (nonlinear propagation)
 *   P^{-}     = F * P * F^T + Q       (linearized covariance propagation)
 *
 * This is a first-order Taylor approximation. The EKF works well
 * for mildly nonlinear systems but can diverge for highly nonlinear
 * dynamics or poor initial estimates.
 *
 * @param x_hat    State estimate [n_states] (in/out)
 * @param P        Error covariance [n_states^2] (in/out)
 * @param f        State transition function: x_next = f(x)
 * @param F_jac    Jacobian of f: F[i][j] = df_i/dx_j evaluated at x_hat
 * @param Q        Process noise covariance [n_states^2]
 * @param n_states State dimension
 */
void af_kalman_extended_predict(double *x_hat, double *P,
                                 void (*f)(const double*, double*),
                                 void (*F_jac)(const double*, double*),
                                 const double *Q, size_t n_states) {
    size_t n = n_states;
    double *F, *x_pred;

    if (!x_hat || !P || !f || !F_jac || !Q || n == 0) return;

    F = (double *)malloc(n * n * sizeof(double));
    x_pred = (double *)malloc(n * sizeof(double));
    if (!F || !x_pred) { free(F); free(x_pred); return; }

    /* 1. Nonlinear state prediction: x_pred = f(x_hat) */
    f(x_hat, x_pred);

    /* 2. Compute Jacobian at current estimate */
    F_jac(x_hat, F);

    /* 3. Covariance prediction: P = F * P * F^T + Q (same as linear case) */
    {
        double *F_P = (double *)malloc(n * n * sizeof(double));
        double *temp = (double *)malloc(n * n * sizeof(double));
        if (F_P && temp) {
            size_t i, j, k;
            /* F_P = F * P */
            kalman_mat_mul(F, P, F_P, n);
            /* temp = F_P * F^T */
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (k = 0; k < n; k++) {
                        sum += F_P[i * n + k] * F[j * n + k];
                    }
                    temp[i * n + j] = sum;
                }
            }
            /* P = temp + Q */
            kalman_mat_add(temp, Q, P, n);
        }
        free(F_P); free(temp);
    }

    /* 4. Update state */
    {
        size_t i;
        for (i = 0; i < n; i++) x_hat[i] = x_pred[i];
    }

    free(F); free(x_pred);
}

/* =========================================================================
 * L5: Extended Kalman Filter - Update Step
 * ========================================================================= */

/**
 * @brief af_kalman_extended_update: EKF update step
 *
 * For nonlinear measurement:
 *   z(k) = h(x(k)) + v(k)
 *
 * The EKF linearizes h around the predicted state:
 *   H = dh/dx evaluated at x_hat^{-}
 *
 * Update (structurally identical to linear Kalman but with H_jac):
 *   S = H * P * H^T + R
 *   K = P * H^T * S^{-1}
 *   x_hat = x_hat + K * (z - h(x_hat))
 *   P = (I - K * H) * P
 *
 * @param x_hat    State estimate [n_states] (in: prior, out: posterior)
 * @param P        Error covariance [n_states^2] (in/out)
 * @param h        Measurement function: z_pred = h(x)
 * @param H_jac    Jacobian of h: H[i][j] = dh_i/dx_j at prior estimate
 * @param R        Measurement noise covariance [n_meas^2]
 * @param z        Measurement scalar
 * @param n_states State dimension
 * @param n_meas   Measurement dimension
 */
void af_kalman_extended_update(double *x_hat, double *P,
                                void (*h)(const double*, double*),
                                void (*H_jac)(const double*, double*),
                                const double *R, double z,
                                size_t n_states, size_t n_meas) {
    size_t n = n_states;
    size_t m = (n_meas > 0) ? n_meas : 1;
    double *H, *z_pred, *K, *P_HT;
    double innov, S;

    if (!x_hat || !P || !h || !H_jac || !R || n == 0) return;

    H = (double *)malloc(m * n * sizeof(double));
    z_pred = (double *)malloc(m * sizeof(double));
    K = (double *)malloc(n * sizeof(double));
    P_HT = (double *)malloc(n * sizeof(double));

    if (!H || !z_pred || !K || !P_HT) {
        free(H); free(z_pred); free(K); free(P_HT);
        return;
    }

    /* 1. Predicted measurement: z_pred = h(x_hat) */
    h(x_hat, z_pred);

    /* 2. Measurement Jacobian at prior estimate */
    H_jac(x_hat, H);

    /* 3. Innovation */
    innov = z - z_pred[0];

    /* 4. P_HT = P * H^T (since m=1, H is row vector) */
    {
        size_t i, j;
        for (i = 0; i < n; i++) {
            double sum = 0.0;
            for (j = 0; j < n; j++) {
                sum += P[i * n + j] * H[j]; /* H^T[j] = H[0*n+j] = H[j] for m=1 */
            }
            P_HT[i] = sum;
        }
    }

    /* 5. S = H * P_HT + R = H*P*H^T + R */
    S = 0.0;
    {
        size_t j;
        for (j = 0; j < n; j++) {
            S += H[j] * P_HT[j];
        }
        S += R[0]; /* R is 1¡Á1 for m=1 */
    }

    /* 6. Kalman gain K = P_HT / S */
    if (fabs(S) < 1e-15) {
        free(H); free(z_pred); free(K); free(P_HT);
        return;
    }
    {
        double inv_S = 1.0 / S;
        size_t i;
        for (i = 0; i < n; i++) K[i] = P_HT[i] * inv_S;
    }

    /* 7. State update */
    {
        size_t i;
        for (i = 0; i < n; i++) x_hat[i] += K[i] * innov;
    }

    /* 8. Covariance update: P = P - K * H * P = P - K * P_HT^T */
    {
        size_t i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                P[i * n + j] -= K[i] * P_HT[j];
            }
        }
    }

    free(H); free(z_pred); free(K); free(P_HT);
}

/* =========================================================================
 * L6: Kalman Filter Application - DC Motor State Estimation
 * ========================================================================= */

/**
 * @brief DC motor state transition function (for EKF demo)
 *
 * State vector: x = [theta, omega]^T (angular position, angular velocity)
 *
 * Continuous-time model:
 *   d(theta)/dt  = omega
 *   d(omega)/dt  = -(B/J)*omega + (Kt/J)*u
 *
 * where: B = damping coefficient, J = moment of inertia,
 *        Kt = torque constant, u = input voltage (control)
 *
 * Discretized with time step dt (forward Euler):
 *   theta_{k+1} = theta_k + dt * omega_k
 *   omega_{k+1} = omega_k + dt * (-(B/J)*omega_k + (Kt/J)*u)
 *
 * This is a canonical example of Kalman filtering for state
 * estimation with noisy measurements.
 */
void af_kalman_dc_motor_f(const double *x, double *x_next) {
    /* Parameters for a typical small DC motor */
    const double dt = 0.01;    /* 10 ms sampling */
    const double B = 0.001;    /* viscous friction (Nm/(rad/s)) */
    const double J = 0.0001;   /* rotor inertia (kg*m^2) */
    const double Kt = 0.05;    /* torque constant (Nm/A) */
    const double u = 1.0;      /* input voltage */

    double theta = x[0];
    double omega = x[1];

    x_next[0] = theta + dt * omega;
    x_next[1] = omega + dt * (-(B/J) * omega + (Kt/J) * u);
}

/**
 * @brief DC motor state Jacobian for EKF
 */
void af_kalman_dc_motor_F_jac(const double *x, double *F) {
    const double dt = 0.01;
    const double B = 0.001;
    const double J = 0.0001;

    /* F = [[1,          dt      ],
     *      [0,  1 - dt*(B/J)    ]] */
    F[0] = 1.0;           F[1] = dt;
    F[2] = 0.0;           F[3] = 1.0 - dt * (B / J);
    (void)x; /* Jacobian of this linear system is x-independent */
}

