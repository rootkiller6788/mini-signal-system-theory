# mini-laplace-z-transform — Laplace & Z Transform Theory

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Partial+ (DC motor, quadrotor dynamics, automotive suspension, GPS anti-aliasing filter)
- **L8**: Partial (Gaver-Stehfest numerical inversion, Chirp Z-Transform)
- **L9**: Partial (documented)
- **Code**: 3450 lines (include/ + src/)

---

## Nine-Level Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | **Complete** | Laplace transform, Z-transform, ROC, poles/zeros, transfer function, impulse/step response, partial fractions |
| **L2** | Core Concepts | **Complete** | Convergence, ROC properties, causality-stability relationship, frequency response, BIBO stability, DC gain, cutoff frequency |
| **L3** | Math Structures | **Complete** | Complex numbers (7 ops), polynomial algebra (11 ops), rational functions (12 ops), pole-zero geometry, partial fraction expansion |
| **L4** | Fundamental Laws | **Complete** | Linearity, time/frequency shifting, convolution theorem, differentiation/integration, Initial Value Theorem, Final Value Theorem, Parseval, Routh-Hurwitz criterion, Jury stability test |
| **L5** | Algorithms | **Complete** | PFE (Heaviside cover-up, derivative method for repeated poles), power series inverse Z, Routh-Hurwitz, Jury test, bilinear transform, numerical Laplace inversion (Gaver-Stehfest) |
| **L6** | Canonical Problems | **Complete** | First-order system analysis, second-order system analysis (damping, overshoot, settling time), digital filter design via bilinear transform, difference equation solving |
| **L7** | Applications | **Partial+** | DC motor speed model, quadrotor altitude dynamics, automotive suspension (quarter-car), GPS anti-aliasing filter design |
| **L8** | Advanced Topics | **Partial** | Numerical inverse Laplace (Gaver-Stehfest), Chirp Z-Transform (CZT), Lyapunov stability for 2x2 systems |
| **L9** | Research Frontiers | **Partial** | Documented: fractional-order transforms, multi-dimensional transforms |

---

## Core Definitions (L1)

| Definition | Symbol | Formula |
|-----------|--------|---------|
| Laplace Transform | F(s) | ∫₀^∞ f(t)e⁻ˢᵗ dt |
| Z-Transform | X(z) | Σ x[n]z⁻ⁿ |
| Region of Convergence | ROC | {s : |F(s)| < ∞} or {z : |X(z)| < ∞} |
| Pole | pₖ | Root of denominator D(s)=0 |
| Zero | zₖ | Root of numerator N(s)=0 |
| Transfer Function | H(s) | Y(s)/X(s) |
| Impulse Response | h(t) | L⁻¹{H(s)} |
| Step Response | y_step(t) | L⁻¹{H(s)/s} |
| Partial Fraction | r/(s-p)ᵏ | Residue/pole decomposition |

## Core Theorems (L4)

| Theorem | Formula | Reference |
|---------|---------|-----------|
| Convolution | L{f∗g} = F(s)G(s) | O&W §9.5.7 |
| Time Differentiation | L{f'} = sF(s) − f(0⁻) | O&W §9.5.5 |
| Time Integration | L{∫f} = F(s)/s | O&W §9.5.6 |
| Frequency Shift | L{e⁻ᵃᵗf} = F(s+a) | O&W §9.5.3 |
| Time Shift (Laplace) | L{f(t−t₀)u(t−t₀)} = e⁻ᵗ⁰ˢF(s) | O&W §9.5.2 |
| Initial Value | f(0⁺) = lim_{s→∞} sF(s) | O&W §9.5.9 |
| Final Value | f(∞) = lim_{s→0} sF(s) | O&W §9.5.9 |
| Parseval (Laplace) | ∫₀^∞|f|²dt = (1/2π)∫|F(jω)|²dω | O&W §9.5.11 |
| Time Shift (Z) | Z{x[n−k]} = z⁻ᵏX(z) | O&S §3.4.2 |
| Parseval (Z) | Σ|x[n]|² = (1/2π)∫|X(eʲω)|²dω | O&S §2.6.5 |
| Routh-Hurwitz | All LHP roots ↔ 1st column sign const | Ogata §5-6 |
| Jury Test | All |z|<1 ↔ Jury conditions | Jury (1964) |

## Core Algorithms (L5)

| Algorithm | Complexity | Key Feature |
|-----------|-----------|-------------|
| PFE Simple Poles | O(n³) | Heaviside cover-up + residue formula |
| PFE Repeated Poles | O(n³) | Derivative method |
| Routh-Hurwitz | O(n²) | First-column sign changes |
| Jury Stability Test | O(n²) | Table construction + constraints |
| Power Series Inverse Z | O(N·deg) | Synthetic division / deconvolution |
| Gaver-Stehfest Inversion | O(M) | Numerical Laplace inversion |
| Bilinear Transform (1st/2nd order) | O(1) | Tustin's method s = (2/T)(1−z⁻¹)/(1+z⁻¹) |

## Canonical Problems (L6)

1. **First-order system analysis** — Settling time, cutoff frequency, step/impulse response
2. **Second-order system analysis** — Damping ratio, natural frequency, overshoot, settling time
3. **Digital filter design via bilinear transform** — Anti-aliasing filter, Butterworth lowpass

## Nine-School Course Mapping

| School | Course | Topic |
|--------|--------|-------|
| **MIT** | 6.003 Signal Processing | Laplace transform, poles/zeros, Bode plots |
| **Stanford** | EE102A Signal Processing | Z-transform, DTFT relationship, filter design |
| **Berkeley** | EE123 Digital Signal Processing | Bilinear transform, IIR filter structures |
| **Illinois** | ECE 310 DSP | Z-transform, difference equations |
| **Michigan** | EECS 351 DSP | Frequency response, stability |
| **Georgia Tech** | ECE 4270 DSP | Digital filter design, structures |
| **TU Munich** | Signal Processing | Laplace transform, control applications |
| **ETH** | 227-0427 Signal Processing | Advanced transform theory |
| **清华** | 信号与系统 | 拉普拉斯变换，Z变换 |

## Key References

- Oppenheim & Willsky, *Signals and Systems* (1997)
- Oppenheim & Schafer, *Discrete-Time Signal Processing* (2010)
- Proakis & Manolakis, *Digital Signal Processing* (2007)
- Ogata, *Modern Control Engineering* (2010)
- Jury, *Theory and Application of the Z-Transform Method* (1964)

---

## File Structure

```
mini-laplace-z-transform/
├── Makefile              # make test compiles and runs all tests
├── README.md             # This file
├── include/              # 6 header files
│   ├── laplace_z_transform_core.h   # Core types (complex, polynomial, rational, PZ, PF)
│   ├── laplace_transform.h          # Laplace transform pairs and properties
│   ├── z_transform.h                # Z-transform pairs and properties
│   ├── partial_fraction.h           # Partial fraction decomposition
│   ├── transfer_function.h          # Transfer function analysis
│   └── stability_analysis.h         # Stability criteria (Routh-Hurwitz, Jury)
├── src/                  # 7 C implementation files (3450 lines total)
│   ├── laplace_z_transform_core.c   # Complex ops, polynomials, rational funcs, PZ, PF
│   ├── laplace_transform.c          # Laplace pairs, properties, theorems
│   ├── z_transform.c                # Z-transform pairs, properties, ROC
│   ├── partial_fraction.c           # PFE: Heaviside, derivative method, pole finding
│   ├── transfer_function.c          # Time/freq response, interconnections, metrics
│   ├── stability_analysis.c         # Routh-Hurwitz, Jury, gain/phase margin, Lyapunov
│   └── inverse_transform.c          # Analytic + numerical inverse transforms, CZT
├── tests/                # 4 test files (all passing)
│   ├── test_core.c                   # Complex, polynomial, rational, PZ tests
│   ├── test_laplace.c                # Laplace transform tests
│   ├── test_z.c                      # Z-transform tests
│   └── test_stability.c             # Stability analysis tests
├── examples/             # 3 end-to-end examples
│   ├── ex1_first_order.c             # First-order system (DC motor)
│   ├── ex2_second_order.c            # Second-order system (quadrotor, suspension)
│   └── ex3_bilinear_filter.c         # Digital filter design (GPS anti-aliasing)
└── docs/                 # Knowledge documentation
```

## Build and Run

```bash
make          # Compile and run all tests + examples
make test     # Run all tests
make examples # Run all examples
make clean    # Clean build artifacts
```
