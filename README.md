# Mini Signal & System Theory

A collection of **from-scratch, zero-dependency C implementations** of university-level signal processing and system theory. Each module maps to MIT, Stanford, and Berkeley courses, bridging theory and practice by translating textbook equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-adaptive-filter](mini-adaptive-filter/) | LMS, RLS, Kalman filter, adaptive FIR/IIR/Lattice/Subband structures | MIT 6.341, Stanford EE264 |
| [mini-filter-theory](mini-filter-theory/) | Butterworth, Chebyshev I/II, Elliptic, Bessel analog/digital filter design | MIT 6.003, Stanford EE102A, Berkeley EE105 |
| [mini-fourier-analysis](mini-fourier-analysis/) | Cooley-Tukey FFT, DFT, convolution, correlation, spectral estimation, window functions | MIT 6.003, Stanford EE261 |
| [mini-laplace-z-transform](mini-laplace-z-transform/) | Laplace transform pairs, Z-transform, rational functions, s-domain to z-domain mapping | MIT 6.003, Stanford EE102A |
| [mini-nonlinear-system](mini-nonlinear-system/) | Bifurcation analysis, chaos, Lyapunov stability, describing functions, MEMS, GPS PLL | MIT 6.450, Stanford EE359, Berkeley EE222 |
| [mini-signal-basis](mini-signal-basis/) | CT/DT signal types, elementary signals, decomposition, DTMF, ECG QRS detection | Stanford EE359, Michigan EECS 455 |
| [mini-system-analysis](mini-system-analysis/) | CT/DT convolution, frequency response, stability analysis, impulse/step response | MIT 6.003, Stanford EE102A |
| [mini-system-identification](mini-system-identification/) | LS/WLS/PEM/ML/IV estimation, RLS, LMS, NLMS, Kalman filter, EKF | MIT 6.435, Stanford EE263 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes
- **Practical demos** — adaptive noise cancellation, filter designer, spectrum analyzer, ECG detector, system identifier, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-adaptive-filter
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-signal-system-theory/
├── mini-adaptive-filter/       # Adaptive filtering: LMS, RLS, Kalman, lattice structures
├── mini-filter-theory/         # Analog/digital filter design: Butterworth, Chebyshev, Elliptic
├── mini-fourier-analysis/      # Fourier analysis: FFT, DFT, spectral estimation, windows
├── mini-laplace-z-transform/   # Laplace & Z transforms: transform pairs, s-to-z mapping
├── mini-nonlinear-system/      # Nonlinear systems: bifurcation, chaos, Lyapunov stability
├── mini-signal-basis/          # Signal fundamentals: types, decompositions, applications
├── mini-system-analysis/       # LTI system analysis: convolution, frequency response, stability
└── mini-system-identification/ # System ID: LS, PEM, ML, RLS, Kalman for parameter estimation
```

## License

MIT
