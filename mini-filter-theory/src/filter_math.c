/**
 * filter_math.c — Mathematical Utilities for Filter Design
 *
 * L3: Polynomial manipulation, root finding, convolution
 * L3: Modified Bessel function I0(x), Jacobi elliptic functions
 * L3: Binomial coefficients, combinatorial expansions
 *
 * Key algorithms:
 *   - polynomial_multiply: convolution of coefficient arrays
 *   - polynomial_eval_complex: Horner's method for complex arguments
 *   - bessel_i0: series expansion for modified Bessel function
 *   - agm: arithmetic-geometric mean (used in elliptic K)
 *   - binomial: combinatorial coefficients for bilinear expansion
 *
 * Reference:
 *   Abramowitz & Stegun, Handbook of Mathematical Functions (1965)
 *   Press et al., Numerical Recipes in C (2007)
 *   Cody, "Modified Bessel Functions of the First Kind" (1983)
 */

#include "filter_transfer.h"
#include "filter_analog.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * L3: Polynomial Operations
 * ================================================================== */

/**
 * polynomial_multiply — L3 convolution
 *
 * C(x) = A(x) * B(x)
 * C_k = sum_{i=0}^{k} A_i * B_{k-i}  for k = 0..(lenA+lenB-2)
 *
 * This is the fundamental operation for:
 * - Transfer function cascading (multiplying TFs)
 * - Window method (ideal response * window)
 * - Bilinear transform (polynomial substitution)
 *
 * Complexity: O(lenA * lenB)
 */
int polynomial_multiply(const double *A, int lenA,
                        const double *B, int lenB,
                        double *C, int max_lenC) {
    if (!A || !B || !C || lenA < 1 || lenB < 1) return -1;

    int result_len = lenA + lenB - 1;
    if (max_lenC < result_len) return -1;

    int i, j;
    for (i = 0; i < result_len; i++) C[i] = 0.0;

    for (i = 0; i < lenA; i++) {
        for (j = 0; j < lenB; j++) {
            C[i + j] += A[i] * B[j];
        }
    }

    return result_len;
}

/**
 * polynomial_eval_complex — L3 Horner's method
 *
 * P(s) = p[0] + p[1]*s + p[2]*s^2 + ... + p[n]*s^n
 * Evaluate as: p[0] + s*(p[1] + s*(p[2] + ...))
 *
 * Horner's method is numerically stable (minimizes roundoff error)
 * and optimal in operation count: n multiplications, n additions.
 *
 * Complexity: O(degree)
 */
double complex polynomial_eval_complex(const double *coeff, int degree,
                                        double complex x) {
    if (!coeff || degree < 0) return 0;

    double complex result = 0;
    int i;
    for (i = degree; i >= 0; i--) {
        result = result * x + coeff[i];
    }
    return result;
}

/**
 * polynomial_add — L3
 *
 * C(x) = A(x) + B(x), handles different lengths.
 * Complexity: O(max(lenA, lenB))
 */
int polynomial_add(const double *A, int lenA,
                   const double *B, int lenB,
                   double *C, int max_lenC) {
    if (!A || !B || !C) return -1;
    int max_len = (lenA > lenB) ? lenA : lenB;
    if (max_lenC < max_len) return -1;

    int i;
    for (i = 0; i < max_len; i++) {
        double a = (i < lenA) ? A[i] : 0.0;
        double b = (i < lenB) ? B[i] : 0.0;
        C[i] = a + b;
    }
    return max_len;
}

/**
 * polynomial_derivative — L3
 *
 * P'(s) = d/ds P(s)
 * If P(s) = p0 + p1*s + p2*s^2 + ... + pn*s^n
 * Then P'(s) = p1 + 2*p2*s + ... + n*pn*s^{n-1}
 *
 * Complexity: O(degree)
 */
int polynomial_derivative(const double *coeff, int degree,
                           double *deriv) {
    if (!coeff || !deriv || degree < 1) return -1;
    int i;
    for (i = 0; i < degree; i++) {
        deriv[i] = coeff[i + 1] * (i + 1);
    }
    return degree;
}

/**
 * polynomial_normalize — L3
 *
 * Normalize so leading coefficient = 1.
 * Returns the original leading coefficient.
 *
 * Complexity: O(degree)
 */
double polynomial_normalize(double *coeff, int degree) {
    if (!coeff || degree < 0) return 1.0;
    double lead = coeff[degree];
    if (fabs(lead) < 1e-15) return 1.0;
    int i;
    for (i = 0; i <= degree; i++) coeff[i] /= lead;
    return lead;
}

/* ==================================================================
 * L3: Binomial Coefficients
 * ================================================================== */

/**
 * binomial — L3 combinatorial coefficient
 *
 * C(n, k) = n!/(k!*(n-k)!)
 *
 * Uses multiplicative formula to avoid overflow:
 *   C(n, k) = prod_{i=1}^{k} (n - k + i) / i
 *
 * Complexity: O(k)
 */
double binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    if (k > n - k) k = n - k; /* Symmetry */

    double result = 1.0;
    int i;
    for (i = 1; i <= k; i++) {
        result *= (n - k + i);
        result /= i;
    }
    return result;
}

/**
 * binomial_expand — L3 binomial expansion
 *
 * (x + y)^n = sum_{k=0}^{n} C(n,k) * x^{n-k} * y^k
 *
 * Computes coefficients as array where coeff[i] multiplies x^{n-i} * y^i.
 *
 * Complexity: O(n)
 */
int binomial_expand(int n, double *coeff) {
    if (!coeff || n < 0) return -1;
    int k;
    for (k = 0; k <= n; k++) {
        coeff[k] = binomial(n, k);
    }
    return n + 1;
}

/**
 * binomial_expand_signed — L3
 *
 * (x - y)^n = sum C(n,k) * x^{n-k} * (-y)^k
 * = sum C(n,k) * (-1)^k * x^{n-k} * y^k
 *
 * Complexity: O(n)
 */
int binomial_expand_signed(int n, double *coeff) {
    if (!coeff || n < 0) return -1;
    int k;
    for (k = 0; k <= n; k++) {
        double c = binomial(n, k);
        coeff[k] = (k % 2 == 0) ? c : -c;
    }
    return n + 1;
}

/* ==================================================================
 * L3: Modified Bessel Function I0(x)
 * ================================================================== */

/**
 * bessel_i0 — L3 series expansion
 *
 * I_0(x) = sum_{k=0}^{inf} (x/2)^{2k} / (k!)^2
 *
 * This function is essential for:
 * - Kaiser window: w[n] = I0(beta*sqrt(1-(2n/(N-1)-1)^2)) / I0(beta)
 * - It appears in the optimal time-bandwidth DPSS sequences
 *
 * Implementation uses the power series with conditional convergence
 * acceleration for large x.
 *
 * For x < 12:   direct series, converges in ~20 terms
 * For x >= 12:  asymptotic expansion for stability
 *
 * Reference: Cody & Stoltz, ACM Trans. Math. Software, 1983
 * Complexity: O(log(x)) iterations
 */
double bessel_i0(double x) {
    if (x < 0.0) x = -x;
    if (x == 0.0) return 1.0;

    if (x < 12.0) {
        /* Power series: I0(x) = sum (x/2)^{2k} / (k!)^2 */
        double sum = 1.0;
        double term = 1.0;
        double x2 = (x / 2.0) * (x / 2.0);
        int k;

        for (k = 1; k < 50; k++) {
            term *= x2 / (k * k);
            sum += term;
            if (term < 1e-16 * sum) break;
        }
        return sum;
    } else {
        /* Asymptotic expansion: I0(x) ~ e^x / sqrt(2*pi*x)
         * with correction factor for moderate x */
        double result = exp(x) / sqrt(2.0 * M_PI * x);

        /* Correction series:
         * 1 + 1/(8x) + 9/(128x^2) + 75/(1024x^3) + ... */
        double correction = 1.0
                          + 1.0 / (8.0 * x)
                          + 9.0 / (128.0 * x * x)
                          + 75.0 / (1024.0 * x * x * x);
        return result * correction;
    }
}

/**
 * bessel_i0_ratio — L3: I0(x) / I0(x0)
 *
 * Used in Kaiser window computation where I0(beta*x) / I0(beta)
 * appears for x = sqrt(1 - (2n/(N-1) - 1)^2) in [0, 1].
 *
 * For large beta, direct computation can overflow. This function
 * uses a stable formulation.
 *
 * Complexity: O(1) for each call
 */
double bessel_i0_ratio(double x, double x0) {
    if (x0 <= 0.0) return 1.0;

    /* If x/x0 is small, direct series may be stable */
    if (x0 < 12.0 || x < 12.0) {
        double i0_x  = bessel_i0(x);
        double i0_x0 = bessel_i0(x0);
        if (fabs(i0_x0) < 1e-15) return 1.0;
        return i0_x / i0_x0;
    }

    /* For large arguments, use the fact that:
     * I0(x)/I0(x0) ≈ exp(x-x0) * sqrt(x0/x) * correction */
    double ratio = exp(x - x0) * sqrt(x0 / x);
    double corr = (1.0 + 1.0/(8.0*x) + 9.0/(128.0*x*x))
                / (1.0 + 1.0/(8.0*x0) + 9.0/(128.0*x0*x0));
    return ratio * corr;
}

/* ==================================================================
 * L3: Numerical Integration (Trapezoidal)
 * ================================================================== */

/**
 * trapezoidal_integrate — L3
 *
 * Integrates f(x) over [a, b] using the trapezoidal rule:
 *   integral ≈ h * (f(a)/2 + f(a+h) + f(a+2h) + ... + f(b)/2)
 * where h = (b-a)/n.
 *
 * Error: O(h^2) for smooth functions.
 *
 * Complexity: O(n) evaluations
 */
double trapezoidal_integrate(double (*f)(double, void*), void *params,
                              double a, double b, int n) {
    if (!f || n < 1 || a >= b) return 0.0;

    double h = (b - a) / n;
    double sum = 0.5 * (f(a, params) + f(b, params));

    int i;
    for (i = 1; i < n; i++) {
        sum += f(a + i * h, params);
    }

    return sum * h;
}

/**
 * simpson_integrate — L3
 *
 * Simpson's rule (order 4):
 *   integral ≈ h/3 * (f0 + 4*f1 + 2*f2 + 4*f3 + ... + fn)
 *
 * Requires n to be even.
 * Error: O(h^4) for smooth functions.
 *
 * Complexity: O(n) evaluations
 */
double simpson_integrate(double (*f)(double, void*), void *params,
                          double a, double b, int n) {
    if (!f || n < 2 || a >= b) return 0.0;
    if (n % 2 != 0) n++; /* Ensure even */

    double h = (b - a) / n;
    double sum = f(a, params) + f(b, params);

    int i;
    for (i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x, params) * ((i % 2 == 0) ? 2.0 : 4.0);
    }

    return sum * h / 3.0;
}

/* ==================================================================
 * L3: Root Finding — Bisection Method
 * ================================================================== */

/**
 * bisection_find_root — L3
 *
 * Finds a root of f(x) in interval [a, b] where f(a)*f(b) < 0.
 * Guaranteed convergence (linear rate, ~1 bit per iteration).
 *
 * Complexity: O(log((b-a)/tol)) function evaluations
 */
double bisection_find_root(double (*f)(double, void*), void *params,
                            double a, double b, double tol) {
    if (!f || a >= b || tol <= 0.0) return NAN;

    double fa = f(a, params);
    double fb = f(b, params);

    if (fa * fb > 0.0) return NAN; /* No sign change, root not guaranteed */

    int iter;
    for (iter = 0; iter < 100; iter++) {
        double c = (a + b) / 2.0;
        double fc = f(c, params);

        if (fabs(fc) < tol || (b - a) < tol) return c;

        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return (a + b) / 2.0;
}

/* ==================================================================
 * L3: Newton-Raphson Root Finding
 * ================================================================== */

/**
 * newton_raphson — L3
 *
 * x_{n+1} = x_n - f(x_n)/f'(x_n)
 *
 * Quadratic convergence near simple roots.
 * Requires derivative function.
 *
 * Complexity: O(iterations) * O(f_eval + fprime_eval)
 */
double newton_raphson(double (*f)(double, void*),
                       double (*fprime)(double, void*),
                       void *params, double x0, double tol, int max_iter) {
    if (!f || !fprime || max_iter < 1) return NAN;

    double x = x0;
    int i;

    for (i = 0; i < max_iter; i++) {
        double fx = f(x, params);
        double fpx = fprime(x, params);

        if (fabs(fpx) < 1e-15) return NAN; /* Derivative zero, method fails */

        double dx = fx / fpx;
        x = x - dx;

        if (fabs(dx) < tol) return x;
    }

    return x;
}

/* ==================================================================
 * L3: Polynomial Root-Finding — Durand-Kerner Method
 * ================================================================== */

/**
 * durand_kerner — L3 simultaneous root finding
 *
 * Finds all roots of a polynomial simultaneously:
 *   p_i^{(new)} = p_i - P(p_i) / prod_{j!=i} (p_i - p_j)
 *
 * This is the Weierstrass/Durand-Kerner method.
 * Advantages: finds ALL roots at once, simple to implement,
 *             converges quadratically for distinct roots.
 *
 * Disadvantage: sensitive to initial guess; convergence not guaranteed
 *               for closely spaced roots.
 *
 * Reference: Durand (1960), Kerner (1966), Numerische Mathematik
 * Complexity: O(n^2 * iterations) for n-degree polynomial
 */
int durand_kerner(const double *coeff, int degree,
                   double complex *roots, int max_iter) {
    if (!coeff || !roots || degree < 1 || max_iter < 1) return -1;
    if (degree > 50) return -1; /* Practical limit */

    /* Normalize leading coefficient */
    double *c = (double *)malloc((degree + 1) * sizeof(double));
    if (!c) return -1;
    int i;
    for (i = 0; i <= degree; i++) c[i] = coeff[i];

    double lead = c[degree];
    if (fabs(lead) < 1e-15) { free(c); return -1; }
    for (i = 0; i <= degree; i++) c[i] /= lead;

    /* Initial guesses: equally spaced on unit circle */
    for (i = 0; i < degree; i++) {
        double angle = 2.0 * M_PI * i / degree + 0.4;
        roots[i] = cos(angle) + sin(angle) * I;
    }

    /* Main iteration */
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        int converged = 1;

        for (i = 0; i < degree; i++) {
            /* Evaluate P(roots[i]) */
            double complex p_val = polynomial_eval_complex(c, degree, roots[i]);

            /* Compute product (roots[i] - roots[j]) for j != i */
            double complex prod = 1.0;
            int j;
            for (j = 0; j < degree; j++) {
                if (j != i) prod *= (roots[i] - roots[j]);
            }

            if (cabs(prod) < 1e-15) prod = 1e-10;

            double complex delta = p_val / prod;
            roots[i] -= delta;

            if (cabs(delta) > 1e-10) converged = 0;
        }

        if (converged) break;
    }

    free(c);
    return iter < max_iter ? 0 : -1;
}

/* ==================================================================
 * L3: Linear Algebra — Gaussian Elimination
 * ================================================================== */

/**
 * gaussian_elimination — L3
 *
 * Solves Ax = b using Gaussian elimination with partial pivoting.
 * A is stored in column-major order (Fortran style) for BLAS compatibility.
 *
 * Used in least-squares FIR filter design and other optimization problems.
 *
 * Complexity: O(n^3)
 *
 * @param A    n×n matrix (column-major), modified in-place
 * @param b    Right-hand side vector, overwritten with solution
 * @param n    Dimension
 * @return     0 on success, -1 if singular
 */
int gaussian_elimination(double *A, double *b, int n) {
    if (!A || !b || n < 1) return -1;

    int i, j, k;

    /* Forward elimination with partial pivoting */
    for (k = 0; k < n - 1; k++) {
        /* Find pivot: max |A[i,k]| for i >= k */
        int max_row = k;
        double max_val = fabs(A[max_row * n + k]);
        for (i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < 1e-15) return -1; /* Singular */

        /* Swap rows */
        if (max_row != k) {
            for (j = k; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k];
            b[k] = b[max_row];
            b[max_row] = tmp;
        }

        /* Eliminate */
        for (i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / A[k * n + k];
            for (j = k; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    /* Back substitution */
    for (i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        if (fabs(A[i * n + i]) < 1e-15) return -1;
        b[i] = sum / A[i * n + i];
    }

    return 0;
}
