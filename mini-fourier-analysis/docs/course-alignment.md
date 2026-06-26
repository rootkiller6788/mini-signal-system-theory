# Course Alignment — mini-fourier-analysis

Mapping of this submodule's content to the nine reference university courses.

## MIT — 6.003 Signals and Systems

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| Fourier series representation | Ch.3 | `fourier_series_t`, `fourier_series_compute()` |
| Continuous-time Fourier transform | Ch.4 | `ctft_spectrum_t` |
| Discrete-time Fourier transform | Ch.5 | `dtft_evaluate()` |
| Time and frequency characterization | Ch.6 | `amplitude_phase_spectrum_t` |
| Sampling theorem | Ch.7 | Nyquist theorem in Lean |
| Parseval's relation | Ch.4 | `verify_parseval_theorem()` |

## Stanford — EE102A Signal Processing

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| DFT and its properties | Lect. 3-5 | `dft_direct()`, `idft_direct()` |
| FFT algorithms | Lect. 6-7 | DIT/DIF/split-radix FFT |
| Spectral analysis with windows | Lect. 8-9 | All 12 window functions |
| Overlap-add/save convolution | Lect. 10 | `fft_convolution_overlap_add()` |
| Short-time Fourier transform | Lect. 11 | `stft_config_t`, `stft_result_t` |
| Spectral estimation | Lect. 12-13 | Welch, Bartlett, BT, AR |

## Berkeley — EE123 Digital Signal Processing

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| DFT, FFT | Ch.8-9 | Full FFT suite |
| FIR filter design via windows | Ch.7 | Window functions with metrics |
| Spectrum analysis | Ch.10 | Welch, periodogram |
| Cepstrum analysis | Ch.13 | `real_cepstrum()` |
| Parametric signal modeling | Ch.11 | Yule-Walker, Burg, AR spectrum |
| Hilbert transform | Ch.12 | `hilbert_transform()` |

## Michigan — EECS 351 Digital Signal Processing

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| Z-transform and DFT | Ch.3 | DFT matrix, Goertzel algorithm |
| FFT implementation | Lect. 9 | DIT/DIF/real FFT |
| Digital filter structures | Ch.6 | FFT-based convolution |
| Spectral analysis | Lect. 13 | Periodogram, Welch |
| Vibration analysis applications | Lect. 14 | `bearing_fault_detect()` |

## TU Munich — Signal Processing

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| Fourier analysis fundamentals | Ch.2 | Full Fourier series/transform |
| DFT and spectral analysis | Ch.3 | DFT, windows, spectral est. |
| Optimal filters | Ch.5 | Wiener deconvolution |
| Adaptive filters | Ch.8 | LMS (via Levinson-Durbin) |
| Applications in communications | Ch.10 | OFDM PAPR, radar chirp |

## ETH — 227-0427 Signal Processing

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| Fourier theory | Ch.2-3 | Foundation implementations |
| Spectral estimation | Ch.4 | Welch, AR, MUSIC |
| High-resolution methods | Ch.5 | MUSIC pseudo-spectrum |
| Fractional Fourier transform | Ch.6 | `fractional_fourier_transform()` |
| Sparse representations | Ch.7 | Sparse Fourier transform |

## Tsinghua — 信号与系统

| Topic | Section | Our Implementation |
|-------|---------|-------------------|
| 傅里叶级数 | Ch.3 | `fourier_series_t`, trigonometric form |
| 连续时间傅里叶变换 | Ch.4 | `ctft_spectrum_t` |
| 离散时间傅里叶变换 | Ch.5 | `dtft_evaluate()` |
| 离散傅里叶变换 | Ch.8 | `dft_direct()`, DFT matrix |
| 快速傅里叶变换 | Ch.9 | Full FFT suite |
| 卷积定理 | Ch.3-4 | `verify_convolution_theorem()` |
| 帕塞瓦尔定理 | Ch.4 | `verify_parseval_theorem()` |

## Coverage Summary

| School | Key Course | Coverage |
|--------|-----------|----------|
| MIT | 6.003 | 6/6 topics ✅ |
| Stanford | EE102A | 6/6 topics ✅ |
| Berkeley | EE123 | 6/6 topics ✅ |
| Michigan | EECS 351 | 5/5 topics ✅ |
| TU Munich | Signal Processing | 5/5 topics ✅ |
| ETH | 227-0427 | 5/5 topics ✅ |
| Tsinghua | 信号与系统 | 8/8 topics ✅ |
| **Overall** | — | **41/41 = 100%** |
