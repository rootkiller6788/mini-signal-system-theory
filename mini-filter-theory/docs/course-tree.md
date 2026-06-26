# Course Tree — mini-filter-theory

## Prerequisite Dependency Graph

```
mini-signal-basis (complex numbers, sinusoids)
    │
    ├── mini-fourier-analysis (frequency domain, spectra)
    │       │
    │       └── mini-filter-theory ◄── YOU ARE HERE
    │               │
    │               ├── mini-laplace-z-transform (s-plane, z-plane)
    │               ├── mini-system-analysis (stability, response)
    │               └── mini-adaptive-filter (LMS, RLS, Kalman)
    │
    └── mini-signal-system-theory (system classification)
```

## Internal Dependency Tree (this module)

```
filter_defs.h (L1: core types)
    │
    ├── filter_transfer.h (L3: polynomial TF, L4: stability)
    │       │
    │       └── filter_math.c (L3: Bessel I0, elliptic K, root finding)
    │
    ├── filter_analog.h (L5: analog prototypes)
    │       │
    │       └── filter_analog.c (Butterworth, Chebyshev I/II, Elliptic, Bessel)
    │
    ├── filter_digital.h (L2: FIR/IIR structures, L5: CIC, half-band)
    │       │
    │       └── filter_digital.c (streaming processing, window FIR design)
    │
    ├── filter_design.h (L5: complete design methods)
    │       │
    │       └── filter_design.c (windows, freq sampling, IIR design, quantization)
    │
    └── filter_analysis.h (L2: response, L3: group delay, L6: metrics)
            │
            └── filter_analysis.c (freq resp, stability, step/impulse)

filter.lean (L4: formal theorems — Routh-Hurwitz, Jury, bilinear, Parseval)
```

## Knowledge Progression

1. **L1 Definitions** → Understand filter types, data structures
2. **L2 Concepts** → Grasp direct forms, cascade/parallel, linear phase
3. **L3 Math** → Master polynomials, Bessel functions, elliptic integrals
4. **L4 Laws** → Apply Routh-Hurwitz, Jury, bilinear stability
5. **L5 Algorithms** → Design Butterworth, Chebyshev, elliptic, window FIR
6. **L6 Problems** → Design complete filters, analyze performance
7. **L7 Applications** → Audio crossover, ECG, communications
8. **L8 Advanced** → Coefficient sensitivity, CIC compensation
9. **L9 Frontiers** → AI-designed filters, quantum, fractional-order
