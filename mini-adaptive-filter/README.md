# mini-adaptive-filter

Adaptive Filter Theory — Implementation of LMS, NLMS, RLS, QR-RLS, Kalman, and Wiener filters in C with Lean 4 formalization.

## Module Status: COMPLETE ?

- **L1-L6**: Complete
- **L7**: Partial (4 applications: ECG noise cancel, echo cancel, channel equalization, GPS multipath)
- **L8**: Partial (QR-RLS, EKF, UKF, beamforming foundations)
- **L9**: Partial (documented in knowledge-graph.md)

| Level | Status | Score |
|-------|--------|-------|
| L1 Definitions | Complete | 2 |
| L2 Core Concepts | Complete | 2 |
| L3 Math Structures | Complete | 2 |
| L4 Fundamental Laws | Complete | 2 |
| L5 Algorithms/Methods | Complete | 2 |
| L6 Canonical Problems | Complete | 2 |
| L7 Applications | Partial | 1 |
| L8 Advanced Topics | Partial | 1 |
| L9 Research Frontiers | Partial | 1 |
| **Total** | | **15/18** |

**Line Count**: include/ + src/ = 3874 lines (threshold: 3000) ?

---

## Core Definitions

### Adaptive Filter Types
| Type | Symbol | Definition |
|------|--------|-----------|
| FIR Transversal | w ∈ R^M | Tapped-delay-line with M tap weights |
| Step-size | μ | Controls LMS convergence speed vs. steady-state error |
| Forgetting factor | λ | RLS memory: 0 < λ ≤ 1 |
| Instantaneous error | e(n) | d(n) - w^T(n)·x(n) |
| Mean-square error | MSE(n) | E[‖e(n)‖2] |
| Excess MSE | J_excess | MSE(∞) - J_min |
| Misadjustment | M | J_excess / J_min |

## Core Theorems

### Wiener-Hopf Equation
```
R · w_opt = p
J_min = σ_d2 - p^T · R?1 · p
```
where R = E[x(n)·x^T(n)], p = E[d(n)·x(n)].

### LMS Convergence (Widrow-Hoff, 1960)
```
Convergence in mean:     0 < μ < 2/λ_max(R)
Convergence in mean-sq:  0 < μ < 2/(3·tr(R))
Misadjustment:           M ≈ (μ/2)·tr(R)  (small μ)
```

### NLMS Convergence
```
Convergence in mean-sq:  0 < μ < 2  (independent of input power)
```

### RLS Convergence
```
RLS converges in ~2M iterations regardless of eigenvalue spread.
Misadjustment: M_RLS ≈ (1-λ)·M/2  (for λ close to 1)
```

### Kalman-Bucy Filter
```
x?? = F·x?                    (state prediction)
P? = F·P·F^T + Q            (covariance prediction)
K  = P?·H^T·(H·P?·H^T+R)?1  (Kalman gain)
x?  = x?? + K·(z - H·x??)     (state update)
P  = (I - K·H)·P?           (covariance update)
```

## Core Algorithms

| Algorithm | Update Rule | Complexity | Reference |
|-----------|------------|------------|-----------|
| LMS | w(n+1)=w(n)+μ·e(n)·x(n) | O(M) | Widrow & Hoff (1960) |
| NLMS | w(n+1)=w(n)+(μ/(ε+‖x‖2))·e(n)·x(n) | O(M) | Nagumo & Noda (1967) |
| Sign-Error LMS | w(n+1)=w(n)+μ·sgn(e(n))·x(n) | O(M) | Widrow & Stearns (1985) |
| Sign-Data LMS | w(n+1)=w(n)+μ·e(n)·sgn(x(n)) | O(M) | Widrow & Stearns (1985) |
| Sign-Sign LMS | w(n+1)=w(n)+μ·sgn(e(n))·sgn(x(n)) | O(M) | Widrow & Stearns (1985) |
| Leaky LMS | w(n+1)=γ·w(n)+μ·e(n)·x(n) | O(M) | Widrow & Stearns (1985) |
| RLS | w(n)=w(n-1)+K(n)·e(n) | O(M2) | Haykin §9.2 |
| QR-RLS (Givens) | Givens rotations on Cholesky factor | O(M2) | Haykin §10.4 |
| APA | w=w+μ·X·(εI+X^TX)?1·e | O(P3+MP2) | Ozeki & Umeda (1984) |
| Kalman Filter | Predict+Update (see theorem above) | O(N3) | Kalman (1960) |
| EKF | Nonlinear f/h + Jacobian linearization | O(N3) | Anderson & Moore (1979) |
| Wiener (Levinson) | Levinson-Durbin recursion | O(M2) | Haykin §2.6 |

## Classic Problems

1. **System Identification**: Estimate unknown system h from (x, d) pairs
2. **Noise Cancellation**: Extract clean signal s[n] from s[n]+n0[n] using reference n1[n]
3. **Echo Cancellation**: Cancel acoustic echo in hands-free telephony
4. **Channel Equalization**: Compensate for multipath distortion in digital comms
5. **Adaptive Line Enhancement**: Detect sinusoid in broadband noise
6. **DC Motor State Estimation**: Kalman filtering for position/velocity estimation

## Course Mapping

| School | Courses | Topics |
|--------|---------|--------|
| MIT | 6.003, 6.450 | Gradient descent, adaptive equalization |
| Stanford | EE264, EE359, EE363 | Full LMS/RLS/Kalman, wireless equalization |
| Berkeley | EE225D, EE123, EE221A | Adaptive signal processing, optimal estimation |
| Illinois | ECE 310, ECE 459 | DSP adaptive filtering, comm equalizer |
| Michigan | EECS 351, EECS 455 | LMS algorithm, channel estimation |
| Georgia Tech | ECE 4270, ECE 6601 | LMS/RLS comparison, adaptive receivers |
| TU Munich | Signal Processing, Communications | Kalman, adaptive algorithms |
| ETH Zurich | 227-0427, 227-0436 | LMS/RLS theory, wireless receivers |
| Tsinghua | 信号与系统, 通信原理 | Wiener filtering, adaptive equalization |

## Build & Test

```bash
make          # Build library and test binary
make test     # Build and run all tests
make examples # Build all 3 examples
make clean    # Remove artifacts
```

## File Structure

```
mini-adaptive-filter/
├── Makefile
├── README.md                          ← This file (COMPLETE ?)
├── include/
│   └── adaptive_filter.h              (440 lines) - Core types, enums, API
├── src/
│   ├── adaptive_filter_core.c         (1467 lines) - Lifecycle, matrix, Wiener, Toeplitz
│   ├── adaptive_lms.c                 (465 lines) - LMS, NLMS, sign variants, leaky
│   ├── adaptive_rls.c                 (478 lines) - RLS, QR-RLS, Kalman gain, Riccati
│   ├── adaptive_kalman.c              (570 lines) - Linear KF, EKF, UKF, DC motor demo
│   ├── adaptive_applications.c        (454 lines) - System ID, noise cancel, echo cancel
│   └── adaptive_formal.lean           (190 lines) - Lean 4 formalization
├── tests/
│   └── test_adaptive_filter.c         (743 lines) - 13 test groups
├── examples/
│   ├── example_system_id.c            - System identification demo
│   ├── example_noise_cancel.c         - Adaptive noise cancellation demo
│   └── example_echo_cancel.c          - Acoustic echo cancellation demo
├── docs/
│   ├── knowledge-graph.md             - L1-L9 knowledge coverage
│   ├── coverage-report.md             - Coverage assessment
│   ├── gap-report.md                  - Missing items & priorities
│   ├── course-alignment.md            - Nine-school course mapping
│   └── course-tree.md                 - Prerequisite dependency tree
├── demos/
├── benches/
└── SKILL.md (symlink to parent)
```

## References

- Haykin, S. (2014). *Adaptive Filter Theory*, 5th ed. Pearson.
- Sayed, A.H. (2008). *Adaptive Filters*. Wiley-IEEE Press.
- Widrow, B. & Stearns, S.D. (1985). *Adaptive Signal Processing*. Prentice-Hall.
- Widrow, B. & Hoff, M.E. (1960). "Adaptive Switching Circuits". IRE WESCON.
- Kalman, R.E. (1960). "A New Approach to Linear Filtering and Prediction Problems". Trans. ASME.
- Widrow, B. et al. (1975). "Adaptive Noise Cancelling: Principles and Applications". Proc. IEEE.
- Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins.
