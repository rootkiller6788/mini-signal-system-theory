# Course Dependency Tree — System Identification

## Prerequisites
This module depends on:

### From mini-fourier-analysis
- Fourier transform (DFT, FFT) → used in ETFE, spectral analysis, frequency response
- Convolution theorem → used in correlation analysis for impulse response

### From mini-laplace-z-transform
- Z-transform → polynomial representations A(q), B(q), etc.
- Transfer functions → model structures (ARX, ARMAX, OE, BJ)
- Bilinear transform → continuous-to-discrete conversion

### From mini-filter-theory
- FIR filters → FIR system identification
- IIR filters → ARX/ARMAX prediction as IIR filtering
- Stability analysis → pole computation for identified models

### From mini-signal-basis
- Orthogonal basis → subspace identification (SVD basis)
- Least squares → core estimation algorithm

### From mini-system-analysis
- State-space representation → SS model identification
- Eigenvalues → pole computation
- Impulse/step response → model validation

## Dependencies Provided To
This module provides foundations for:

### To mini-adaptive-filter
- RLS, LMS algorithms
- Kalman filter for parameter estimation

### To mini-nonlinear-system
- NARX, Hammerstein, Wiener models
- Nonlinear optimization (PEM)

### To mini-signal-basis
- SVD application in subspace methods
- Matrix factorization techniques

## Cross-Cutting Dependencies
- Linear algebra: Cholesky, QR, SVD → shared with multiple modules
- Optimization: Gauss-Newton, Levenberg-Marquardt → nonlinear optimization
- Statistics: Information criteria, hypothesis testing → model validation
