# Knowledge Graph: mini-adaptive-filter

## L1: Definitions (Complete)
- Adaptive filter (FIR/IIR/Lattice/Subband/Freq-domain/Volterra structures)
- Tap-weight vector (w ¡Ê R^M)
- Filter order M (number of adaptive coefficients)
- Step-size parameter mu (LMS convergence control)
- Forgetting factor lambda (RLS memory, 0 < lambda <= 1)
- Regularization parameter delta (RLS initialization: P(0) = delta^{-1} * I)
- Leakage factor gamma (leaky LMS stabilization, 0 < gamma <= 1)
- Instantaneous error e(n) = d(n) - w^T(n)*x(n)
- Mean-square error MSE(n) = E[|e(n)|^2]
- Excess MSE J_excess = MSE(inf) - J_min
- Misadjustment M = J_excess / J_min
- Wiener solution w_opt = R^{-1} * p
- Minimum MSE J_min = sigma_d^2 - p^T * R^{-1} * p
- Autocorrelation matrix R = E[x(n)*x^T(n)] (Toeplitz)
- Cross-correlation vector p = E[d(n)*x(n)]
- Eigenvalue spread chi(R) = lambda_max / lambda_min
- Convergence time tau_mse (LMS: ¡Ö 1/(2*mu*sigma_x^2))
- Signal statistics (af_statistics_t)
- Performance metrics (af_performance_t)
- Complex adaptive filter (I/Q representation)

## L2: Core Concepts (Complete)
- Wiener filter theory (optimal linear filter in MSE sense)
- Gradient descent optimization
- Stochastic gradient approximation (LMS)
- Mean-square error surface (quadratic in w)
- Convergence analysis (mean and mean-square)
- Tracking capability vs. steady-state error trade-off
- Exponentially weighted least squares (RLS cost function)
- Matrix inversion lemma (Woodbury identity for RLS derivation)
- Kalman filtering as recursive Bayesian estimation
- Orthogonality principle (E[e_opt(n)*x(n-k)] = 0)

## L3: Mathematical Structures (Complete)
- Toeplitz matrix representation (compact O(N) storage)
- Matrix operations: multiplication, transpose, inverse (Gauss-Jordan)
- Cholesky decomposition (L*L^T for symmetric positive definite)
- Vector operations: dot product, norm, scale, add, sub
- Power iteration for eigenvalue estimation
- Gershgorin circle theorem (eigenvalue bounds)
- Rank-1 matrix update (RLS Riccati equation)
- Givens rotation (QR-RLS numerical stability)

## L4: Fundamental Laws (Complete)
- Wiener-Hopf equation: R*w_opt = p
- LMS convergence in mean: 0 < mu < 2/lambda_max
- LMS convergence in mean-square: 0 < mu < 2/(3*tr(R))
- NLMS convergence bound: 0 < mu < 2 (independent of input power)
- Misadjustment formula: M ¡Ö (mu/2)*tr(R) (LMS, small mu)
- RLS convergence: ~2M iterations (independent of eigenvalue spread)
- Orthogonality principle at Wiener optimum
- Kalman-Bucy filter: minimum variance unbiased estimation
- Woodbury matrix identity: (A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1}+VA^{-1}U)^{-1}VA^{-1}
- Riccati equation for P matrix update in RLS

## L5: Algorithms/Methods (Complete)
- LMS (Widrow-Hoff, 1960): w(n+1) = w(n) + mu*e(n)*x(n)
- NLMS (Nagumo & Noda, 1967): normalized step-size
- Sign-Error LMS: sgn(e(n)) * x(n)
- Sign-Data LMS: e(n) * sgn(x(n))
- Sign-Sign LMS: sgn(e(n)) * sgn(x(n))
- Leaky LMS: gamma*w(n) + mu*e(n)*x(n)
- Block LMS (frequency-domain concepts)
- RLS: recursive least squares with forgetting
- QR-RLS: Givens rotation-based (numerically stable)
- Affine Projection Algorithm (APA)
- Linear Kalman filter (predict + update)
- Extended Kalman filter (Jacobian linearization)
- Unscented Kalman filter (sigma-point framework, structure defined)
- Levinson-Durbin recursion (Toeplitz system solver, O(M^2))
- Power iteration for eigenvalue estimation

## L6: Canonical Problems (Complete)
- System identification (unknown plant estimation via LMS/RLS)
- Adaptive noise cancellation (Widrow architecture, 1975)
- Acoustic echo cancellation (NLMS, hands-free telephony)
- Channel equalization (adaptive decision-feedback for digital comms)
- Adaptive line enhancement (ALE: enhance sinusoid in noise)
- DC motor state estimation (Kalman filter application)

## L7: Applications (Partial+)
- ECG 60 Hz power-line noise cancellation (biomedical signal processing)
- Acoustic echo cancellation for hands-free telephony
- Channel equalization for 5G/WiFi digital communications
- GPS receiver multipath mitigation (channel equalization application)

## L8: Advanced Topics (Partial)
- QR-RLS for numerically robust adaptive filtering
- Extended Kalman filter for nonlinear system identification
- Unscented Kalman filter framework (sigma-point selection)
- Adaptive beamforming (foundations via complex LMS)

## L9: Research Frontiers (Partial - documented)
- AI/Deep learning based adaptive filtering
- Adaptive filtering for 6G RIS (Reconfigurable Intelligent Surfaces)
- Quantum adaptive filtering
- Distributed adaptive filtering for IoT sensor networks
