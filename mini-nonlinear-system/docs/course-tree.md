# Course Tree — mini-nonlinear-system

## Prerequisites (Dependencies)
```
Linear Algebra → Eigenvalues/Eigenvectors → Lyapunov Indirect Method
    ↓
Differential Equations → Phase Plane → Bifurcation Theory
    ↓
Linear Systems Theory → Transfer Functions → Circle/Popov Criteria
    ↓
Fourier Analysis → Describing Functions → Harmonic Balance
    ↓
Probability Theory → Stochastic Processes → Chaos Analysis
```

## Knowledge Dependency Graph
```
L1: Definitions
 ├─ Nonlinearity types, state-space, Lyapunov function
 │
 ├─► L2: Core Concepts
 │    ├─ Phase plane, limit cycles, bifurcations, chaos
 │    │
 │    ├─► L3: Math Structures
 │    │    ├─ Jacobian, eigenvalues, GFRF
 │    │    │
 │    │    ├─► L4: Fundamental Laws
 │    │    │    ├─ Lyapunov theorems
 │    │    │    ├─ Hartman-Grobman
 │    │    │    ├─ Poincaré-Bendixson
 │    │    │    ├─ Circle/Popov criteria
 │    │    │    └─ Bendixson/Dulac criteria
 │    │    │
 │    │    └─► L5: Algorithms
 │    │         ├─ RK4 / Dormand-Prince
 │    │         ├─ QR eigensolver
 │    │         ├─ DF analytical / numerical
 │    │         ├─ Newton-Raphson / continuation
 │    │         ├─ Lyapunov exponents (Benettin)
 │    │         ├─ Correlation dimension
 │    │         └─ Volterra evaluation / identification
 │    │
 │    └─► L6: Canonical Problems
 │         ├─ Van der Pol (limit cycle)
 │         ├─ Lorenz (chaos)
 │         ├─ PLL (lock range)
 │         ├─ Duffing (nonlinear resonance)
 │         └─ Chua's circuit (chaos)
 │
 └─► L7: Applications
      ├─ GPS carrier tracking
      ├─ 5G NR PHY
      ├─ Amplifier saturation
      └─ Sigma-delta ADC
      
L8: Advanced Topics
 ├─ Kaplan-Yorke dimension
 ├─ Surrogate data testing
 ├─ Delay embedding
 └─ 0-1 chaos test

L9: Research Frontiers
 ├─ 6G RIS nonlinear optimization
 ├─ Quantum nonlinear systems
 └─ Semantic communication
