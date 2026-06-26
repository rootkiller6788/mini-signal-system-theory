# Course Dependency Tree - mini-system-analysis

## Prerequisites

```
Calculus (MIT 18.01/02)
  |
  +-- Differential Equations (MIT 18.03)
  |     |
  |     +-- Laplace Transform
  |     +-- Z-Transform
  |     +-- Matrix Exponential
  |
  +-- Linear Algebra (MIT 18.06)
        |
        +-- Matrix operations (multiply, invert, transpose)
        +-- Eigenvalues and eigenvectors
        +-- Rank, determinant, trace
        +-- Positive definiteness
```

## Module Dependency Tree

```
mini-system-analysis
  |
  L1: Definitions
     |-- Signal types (CT/DT)
     |-- System classification (LTI/LTV/NL)
     |-- Transfer function representation
     |-- Pole-zero representation
     |-- Frequency response data structures
  |
  L2: Core Concepts
     |-- Convolution (CT/DT/circular/Fast/block)
     |   depends_on: L1 (signal types)
     |
     |-- Correlation (cross/auto/normalized)
     |   depends_on: L1, L2-convolution
     |
     |-- Frequency Response (magnitude/phase/bandwidth/Q)
     |   depends_on: L1, L3 (TF evaluation)
     |
     |-- Step/Impulse Relationship
     |   depends_on: L1, L2-convolution
  |
  L3: Math Structures
     |-- Transfer Function Evaluation
     |   depends_on: L1 (TF struct)
     |
     |-- Pole-Zero Analysis
     |   depends_on: L1, L3-TF, L5 (Newton)
     |
     |-- CT-to-DT Transforms
     |   depends_on: L1, L3-TF, L4 (stability)
     |
     |-- State-Space Representation
     |   depends_on: L1, Linear Algebra prereq
     |
     |-- Controllability/Observability
     |   depends_on: L3-state_space
     |
     |-- Matrix Exponential
     |   depends_on: L3-state_space, L5 (Pade)
     |
     |-- Lyapunov Equation
     |   depends_on: L3-state_space, L5 (Gauss-Jordan)
  |
  L4: Fundamental Laws
     |-- BIBO Stability
     |   depends_on: L1, L3 (pole-zero)
     |
     |-- Routh-Hurwitz Criterion
     |   depends_on: L1 (TF coefficient access)
     |
     |-- Jury Stability Test
     |   depends_on: L1 (DT TF)
     |
     |-- Nyquist Criterion
     |   depends_on: L3 (frequency response)
     |
     |-- Gain/Phase Margins
     |   depends_on: L3 (frequency response)
     |
     |-- Root Locus
     |   depends_on: L3 (pole finding), L5 (Newton)
     |
     |-- Lyapunov Direct Method
     |   depends_on: L3 (state-space), L5 (Gauss-Jordan)
     |
     |-- Feedback Stability
     |   depends_on: L4 (Routh-Hurwitz)
  |
  L5: Algorithms
     |-- Horner Method (O(n))
     |-- Newton-Raphson (O(n*iter))
     |-- Synthetic Deflation (O(n))
     |-- Faddeev-LeVerrier (O(n^4))
     |-- Gauss-Jordan Elimination (O(n^3))
     |-- Pade Approximation (O(n^3))
     |-- Runge-Kutta 4 (O(n^2)/step)
     |-- DFT/IDFT (O(N^2))
     |-- Matrix Rank (Gaussian elim, O(n^3))
  |
  L6: Canonical Problems
     |-- ex1: Step Response (1st/2nd order)
     |   depends_on: L1, L2, L3
     |
     |-- ex2: Stability Analysis
     |   depends_on: L1, L3, L4
     |
     |-- ex3: State-Space Simulation
     |   depends_on: L1, L3, L4, L5
     |
     |-- demo: Bode/Filter Analysis
     |   depends_on: L1, L3, L2 (freq. resp.)
  |
  L7: Applications (Partial+)
     |-- Control system stability margins
     |-- Mass-spring-damper vibration
     |-- 60 Hz notch filter
  |
  L8: Advanced Topics (Partial+)
     |-- Lyapunov stability (P solution)
     |-- Matrix exponential (scaling+squaring)
  |
  L9: Research Frontiers (Partial)
     |-- Compressed sensing
     |-- 6G RIS
     |-- Quantum control
```

## Build Order (Topological Sort)

1. system_defs.c (L1 — no dependencies)
2. convolution.c (L2 — depends on L1)
3. transfer_function.c (L3 — depends on L1)
4. frequency_response.c (L2/L5 — depends on L1, L3)
5. state_space.c (L3 — depends on L1)
6. stability.c (L4 — depends on L1, L3, L4-state_space)
7. tests/ (depends on all src/)
8. examples/ (depends on all src/)
9. demos/ (depends on all src/)
10. benches/ (depends on all src/)
