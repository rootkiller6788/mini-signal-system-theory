# Course Dependency Tree: mini-adaptive-filter

## Prerequisites (what you need to know before this module)

```
Linear Algebra (matrices, eigenvalues, least squares)
    |
Probability Theory (random variables, expectation, covariance)
    |
Digital Signal Processing (FIR filters, convolution, z-transform)
    |
Optimization Theory (gradient descent, cost functions)
    |
    +---> mini-adaptive-filter (THIS MODULE)
              |
              +---> mini-filter-theory (FIR/IIR design)
              +---> mini-fourier-analysis (frequency-domain adaptive filters)
              +---> mini-system-identification (advanced ID methods)
```

## Internal Dependency Graph (within this module)

```
L1: Definitions (structs, enums)  <-- Foundation for everything
    |
L2: Core Concepts (Wiener, gradient descent)
    |
    +---> L3: Math Structures (matrix/vector ops)
    |         |
    |         +---> L4: Fundamental Laws (convergence bounds)
    |                   |
    +---> L5: Algorithms
              |
              +---> LMS family (direct, NLMS, sign, leaky, block)
              +---> RLS family (standard, QR-RLS)
              +---> Kalman family (linear, EKF, UKF)
              +---> APA
              |
              +---> L6: Canonical Problems
                        |
                        +---> System Identification
                        +---> Noise Cancellation
                        +---> Echo Cancellation
                        +---> Channel Equalization
                        +---> Adaptive Line Enhancement
                        |
                        +---> L7: Applications
                                  |
                                  +---> ECG denoising
                                  +---> Hands-free telephony
                                  +---> 5G/WiFi equalization
```

## Course Prerequisites (nine-school curriculum)

1. **Linear Algebra** (any school): Matrix operations, eigenvalues, pseudoinverse
2. **Probability & Statistics** (any school): Random variables, correlation, MSE
3. **Signals & Systems** (MIT 6.003 / Stanford EE102A): FIR filters, convolution
4. **DSP** (Berkeley EE123 / Illinois ECE 310): Discrete-time filtering, spectral analysis
5. **Optimization** (Stanford EE364a optional): Gradient methods, least squares
