# Course Dependency Tree — mini-fourier-analysis

## Prerequisites

This submodule assumes knowledge from:
- **Calculus**: Integration, differentiation, complex numbers, infinite series
- **Linear Algebra**: Vector spaces, inner products, matrix multiplication, eigenvalues
- **Probability**: Random variables, expectation, variance (for spectral estimation)
- **Signals & Systems basics**: LTI systems, impulse response, frequency response

## Internal Dependency Graph

```
┌─────────────────────────────────────────────────────┐
│                  fourier_core.h/.c                   │
│   L1: Fourier series, DFT, IDFT, amplitude/phase    │
│   L3: Complex basis, DFT matrix                      │
│   L4: Parseval, Riemann-Lebesgue (verification)     │
└──────────┬──────────────────────┬───────────────────┘
           │                      │
           ▼                      ▼
┌──────────────────────┐  ┌──────────────────────┐
│  fourier_fft.h/.c    │  │ fourier_window.h/.c  │
│  L5: DIT/DIF/split-  │  │ L5: 12 window types │
│      radix/Goertzel  │  │ L2: Window metrics   │
│      Bluestein/Zoom  │  │ L4: Dirichlet kernel │
│  L6: Fast conv       │  └──────────┬───────────┘
└──────────┬───────────┘             │
           │                         │
           ▼                         ▼
┌─────────────────────────────────────────────────────┐
│              fourier_spectrum.h/.c                   │
│   L6: Welch, Bartlett, Blackman-Tukey, AR, MUSIC    │
│   L5: Levinson-Durbin, Burg MEM                     │
│   L6: Coherence, cepstrum, bandwidth                │
└──────────┬──────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────┐
│           fourier_convolution.h/.c                   │
│   L3: Linear/circular convolution, correlation      │
│   L4: Convolution theorem, Wiener-Khinchin           │
│   L6: Wiener deconv, Hilbert transform, analytic    │
│   L8: Fractional Fourier Transform (FrFT)           │
└──────────┬──────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────┐
│          fourier_applications.h/.c                   │
│   L7: Audio, radar, vibration, GPS, OFDM, power     │
│   L8: Sparse Fourier Transform                      │
└─────────────────────────────────────────────────────┘

                    ┌──────────────────┐
                    │ fourier_lean.lean│
                    │ L1-L6 formalized │
                    │ Pure Lean 4 core │
                    └──────────────────┘
```

## Module Dependencies Within Parent

```
0. mini-signal-system-theory/
├── mini-fourier-analysis/     ← THIS MODULE
│   └── (foundation: Fourier methods)
├── mini-laplace-transform/     ← depends on this (complex analysis)
├── mini-z-transform/            ← depends on this (discrete-time)
└── mini-system-response/       ← depends on this (convolution, freq response)
```

## Knowledge Flow

1. **L1 Definitions** → **L2 Concepts** → **L3 Math Structures**
   - Understanding definitions leads to grasping concepts
   - Concepts require mathematical formalization

2. **L4 Theorems** ← **L3 Math Structures**
   - Theorems are expressed in terms of structures
   - Verification bridges math ↔ code

3. **L5 Algorithms** ← **L4 Theorems** + **L3 Structures**
   - Algorithms exploit theorem properties (e.g., convolution theorem → FFT convolution)
   - Algorithms manipulate math structures efficiently

4. **L6 Canonical Problems** ← **L5 Algorithms** + **L2 Concepts**
   - Problems solved by combining algorithms
   - Domain knowledge (concepts) guides problem formulation

5. **L7 Applications** ← **L6 Problems** + domain context
   - Applications are instantiations of canonical problems
   - Real-world constraints add complexity

6. **L8 Advanced Topics** ← independent innovations
   - Build on foundations but introduce new principles
   - May require additional math (e.g., compressed sensing needs convex optimization)

7. **L9 Research Frontiers** ← speculative extensions
   - Future directions, not yet standardized
   - Documented for awareness
