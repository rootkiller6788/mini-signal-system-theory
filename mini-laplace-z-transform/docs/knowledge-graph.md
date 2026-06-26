# Knowledge Graph — Laplace & Z Transform Theory

## L1: Definitions (Complete)
- Laplace transform: F(s) = ∫₀^∞ f(t)e^{-st} dt
- Unilateral vs bilateral Laplace transform
- Z-transform: X(z) = Σ_{n=-∞}^{∞} x[n]z^{-n}
- Region of Convergence (ROC)
- Pole: root of denominator polynomial D(s)=0
- Zero: root of numerator polynomial N(s)=0
- Transfer function H(s) = Y(s)/X(s), H(z) = Y(z)/X(z)
- Impulse response h(t), h[n]
- Step response y_step(t), y_step[n]
- Frequency response H(jω), H(e^{jω})
- Partial fraction expansion terms r/(s-p)^k

## L2: Core Concepts (Complete)
- Convergence of Laplace integral
- ROC properties for causal/anti-causal/two-sided signals
- Relationship Laplace ↔ Fourier (s=jω)
- Relationship Z-transform ↔ DTFT (z=e^{jω})
- Causality condition and ROC
- BIBO stability: all poles in LHP (CT) or inside unit circle (DT)
- Frequency response from pole-zero geometry
- DC gain: H(0) or H(1)
- -3dB cutoff frequency
- Transient response metrics: settling time, overshoot, time constant
- Proper vs strictly proper rational functions

## L3: Mathematical Structures (Complete)
- Complex number algebra (+, -, *, /, mag, phase, conj, exp, sqrt, pow)
- Polynomial representation and algebra (add, sub, mul, div, derivative, integral)
- Rational function representation and algebra
- Polynomial GCD (Euclidean algorithm)
- Pole-zero geometric representation
- Companion matrix eigenvalue method for root-finding
- Newton's method for polynomial roots
- Partial fraction data structure
- Binomial coefficient computation
- Synthetic division / polynomial deflation

## L4: Fundamental Laws (Complete)
- Linearity of Laplace/Z transforms
- Time shifting (Laplace): L{f(t-t₀)u(t-t₀)} = e^{-t₀s}F(s)
- Frequency shifting: L{e^{-at}f(t)} = F(s+a)
- Time scaling: L{f(at)} = (1/a)F(s/a)
- Convolution theorem: L{f*g} = F(s)G(s)
- Time differentiation: L{f'} = sF(s) - f(0⁻)
- Time integration: L{∫₀^t f(τ)dτ} = F(s)/s
- Initial Value Theorem: f(0⁺) = lim_{s→∞} sF(s)
- Final Value Theorem: f(∞) = lim_{s→0} sF(s)
- Parseval's relation (Laplace)
- Z-transform time shift (right/left)
- Z-transform initial/final value theorems
- Z-transform Parseval
- Routh-Hurwitz stability criterion
- Jury stability test
- Series/parallel/feedback interconnection laws

## L5: Algorithms/Methods (Complete)
- Partial fraction decomposition (Laplace): Heaviside cover-up + derivative method
- Partial fraction decomposition (Z-domain)
- Residue computation for simple poles: r = N(p)/D'(p)
- Residue computation for repeated poles: derivative formula
- Polynomial root-finding via Newton's method with deflation
- Routh array construction
- Jury table construction
- Bilinear transform (Tustin): s = (2/T)(1-z⁻¹)/(1+z⁻¹)
- Bilinear transform for 1st and 2nd-order systems
- Inverse Laplace via partial fractions + table lookup
- Inverse Z via partial fractions
- Inverse Z via power series (long division)
- Inverse Z via residue method
- Numerical Laplace inversion (Gaver-Stehfest)
- Chirp Z-Transform (CZT)
- Difference equation solution (Direct Form I)

## L6: Canonical Problems (Complete)
- First-order system: step/impulse/frequency response, settling time, cutoff
- Second-order system: damping ratio, natural frequency, overshoot, settling time
- Digital Butterworth lowpass filter design via bilinear transform
- Pre-warping for bilinear transform
- Stability analysis of feedback systems
- Characteristic equation from closed-loop transfer function

## L7: Applications (Partial+)
- DC motor speed model (1st-order, Detroit manufacturing)
- Quadrotor altitude dynamics (2nd-order, mass-spring-damper)
- Automotive suspension quarter-car model (Toyota class vehicles)
- GPS L1 C/A code anti-aliasing filter design (GPS receiver application)

## L8: Advanced Topics (Partial)
- Numerical inverse Laplace transform (Gaver-Stehfest method, Stehfest 1970)
- Chirp Z-Transform for zoomed frequency analysis (Rabiner et al. 1969)
- Lyapunov stability analysis for linear systems (Faddeev-Leverrier algorithm)

## L9: Research Frontiers (Partial)
- Fractional-order Laplace transforms
- Multi-dimensional Z-transforms
- Wavelet-domain transforms as generalizations
