# Course Alignment - mini-system-analysis

| School | Course | Chapters | Key Topics Covered |
|--------|--------|----------|-------------------|
| **MIT** | 6.003 Signal Processing | Ch.1-3, Ch.6, Ch.9-12 | LTI systems, convolution, Fourier/Laplace, stability, state-space |
| **Stanford** | EE102A Signal Processing | Ch.1-5, Ch.7-10 | CT/DT signals, frequency response, Z-transform, filter design |
| **Berkeley** | EE123 Digital Signal Processing | Ch.1-8 | DT convolution, DFT/FFT, filter structures, sampling |
| **Illinois** | ECE 310 DSP | Ch.1-6 | Discrete-time systems, z-transform, FIR/IIR design |
| **Michigan** | EECS 351 DSP | Ch.1-5 | Signal representation, convolution, sampling theory |
| **Georgia Tech** | ECE 4270 DSP | Ch.1-7 | Time/frequency analysis, DFT, filter banks |
| **TU Munich** | Signal Processing | Ch.1-6 | CT/DT systems, Laplace/Z transform, stability |
| **ETH** | 227-0427 Signal Processing | Ch.1-8 | Signal spaces, transforms, sampling, filter design |
| **清华** | 信号与系统 (Signals & Systems) | Ch.1-8 | Linear systems, convolution, Fourier/Laplace/Z, state-space |

## Chapter-by-Chapter Mapping

| Chapter | Topic | Implementation |
|---------|-------|---------------|
| Signals & Classification | CT/DT signals, energy/power | system_defs.c (L1) |
| LTI Systems | Impulse response, convolution | convolution.c (L2) |
| Fourier Series/Transform | Frequency domain representation | transfer_function.c (L3) |
| Laplace Transform | Transfer functions, poles/zeros, stability | transfer_function.c, stability.c (L3-L4) |
| Z-Transform | DT transfer functions, DT stability | transfer_function.c, stability.c (L3-L4) |
| Frequency Response | Bode plots, filter types, bandwidth | frequency_response.c (L2-L5) |
| Stability Analysis | Routh-Hurwitz, Jury, Nyquist, Lyapunov | stability.c (L4) |
| State-Space Methods | Controllability, observability, simulation | state_space.c (L3) |
| Sampling | A/D conversion, aliasing, reconstruction | system_defs.c (L1) |
| DT Filter Design | Bilinear, impulse invariance, ZOH | transfer_function.c (L3) |
| Feedback Systems | Closed-loop stability, margins | stability.c (L4) |
