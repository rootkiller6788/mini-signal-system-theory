# Knowledge Graph — mini-fourier-analysis

## L1: Definitions (Complete)

| # | Definition | Symbol | C Representation | Status |
|---|-----------|--------|-----------------|--------|
| 1 | Fourier series (trigonometric) | `x(t) = a₀/2 + Σ a_n cos + b_n sin` | `fourier_series_t` | ✅ |
| 2 | Fourier series (complex) | `x(t) = Σ c_n e^{jnω₀t}` | `fourier_series_complex_t` | ✅ |
| 3 | Continuous Fourier Transform | `X(ω) = ∫ x(t) e^{-jωt} dt` | `ctft_spectrum_t` | ✅ |
| 4 | Discrete-Time Fourier Transform | `X(e^{jω}) = Σ x[n] e^{-jωn}` | `dtft_spectrum_t` | ✅ |
| 5 | Discrete Fourier Transform (DFT) | `X[k] = Σ x[n] W_N^{kn}` | `dft_result_t` | ✅ |
| 6 | Amplitude/Phase spectrum | `|X[k]|, φ[k]` | `amplitude_phase_spectrum_t` | ✅ |
| 7 | Power Spectral Density | PSD | `psd_result_t` | ✅ |
| 8 | Energy Spectral Density | ESD | `esd_result_t` | ✅ |
| 9 | Bandwidth (3dB, ENB, occupied) | — | `bandwidth_metrics_t` | ✅ |
| 10 | Frequency bin descriptor | `k → f_k = k·fs/N` | `freq_bin_t` | ✅ |
| 11 | Window function representation | `w[n]` | `window_config_t` | ✅ |
| 12 | STFT / Spectrogram | `X[m,k]` | `stft_config_t, stft_result_t` | ✅ |
| 13 | Chirp Z-Transform config | — | `czt_params_t` | ✅ |
| 14 | FFT plan | — | `fft_plan_t` | ✅ |
| 15 | Real cepstrum | `c[n] = IDFT{log|DFT{x}|}` | — | ✅ |
| 16 | Magnitude squared coherence | `C_xy(f) ∈ [0,1]` | — | ✅ |

## L2: Core Concepts (Complete)

| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Time-frequency duality | DIT/DIF FFT, inverse transform round-trip | ✅ |
| 2 | Orthogonality of DFT basis | DFT matrix computation, symmetry verification | ✅ |
| 3 | Spectral leakage | All window functions with sidelobe analysis | ✅ |
| 4 | Conjugate symmetry | `verify_conjugate_symmetry()` | ✅ |
| 5 | Windowing trade-off (mainlobe vs sidelobe) | Window metrics computation | ✅ |
| 6 | Coherent gain | `window_compute_metrics()` | ✅ |
| 7 | Equivalent noise bandwidth (ENBW) | `window_compute_metrics()` | ✅ |
| 8 | Processing loss | `window_compute_metrics()` | ✅ |
| 9 | Scalloping loss | `window_compute_metrics()` | ✅ |
| 10 | Overlap processing | Welch's method with configurable overlap | ✅ |
| 11 | Zero-padding for linear convolution | FFT-based convolution (overlap-add/save) | ✅ |
| 12 | Frequency resolution Δf = fs/N | All DFT functions | ✅ |
| 13 | Nyquist frequency fs/2 | Real FFT half-spectrum output | ✅ |

## L3: Mathematical Structures (Complete)

| # | Structure | C Implementation | Status |
|---|----------|-----------------|--------|
| 1 | Complex exponential basis `{e^{jωt}}` | `fourier_complex_t` throughout | ✅ |
| 2 | Trigonometric basis `{cos, sin}` | `fourier_series_t` | ✅ |
| 3 | DFT matrix (Vandermonde) | `dft_matrix_generate()` | ✅ |
| 4 | Roots of unity `W_N = e^{-j·2π/N}` | FFT twiddle factors | ✅ |
| 5 | Linear convolution | `convolution_linear()` | ✅ |
| 6 | Circular convolution | `convolution_circular()` | ✅ |
| 7 | Cross-correlation | `cross_correlation()` | ✅ |
| 8 | Autocorrelation | `autocorrelation()` | ✅ |
| 9 | Toeplitz matrix (Yule-Walker) | `levinson_durbin()` | ✅ |
| 10 | Eigenvalue decomposition (MUSIC) | `music_pseudospectrum()` | ✅ |
| 11 | Hilbert space L² | Parseval/Plancherel verification | ✅ |
| 12 | Cepstral domain | `real_cepstrum()`, `cepstral_liftering()` | ✅ |

## L4: Fundamental Laws (Complete)

| # | Theorem | C Verification | Lean Statement | Status |
|---|---------|---------------|----------------|--------|
| 1 | Parseval's theorem | `verify_parseval_theorem()` | `parseval_theorem` | ✅ |
| 2 | Convolution theorem | `verify_convolution_theorem()` | `convolution_theorem` | ✅ |
| 3 | Wiener-Khinchin theorem | `verify_wiener_khinchin()` | `wiener_khinchin_theorem` | ✅ |
| 4 | Plancherel theorem | — | `plancherel_theorem` | ✅ |
| 5 | Nyquist-Shannon sampling theorem | — | `nyquist_limit` | ✅ |
| 6 | Riemann-Lebesgue lemma | `riemann_lebesgue_ratio()` | `riemann_lebesgue_lemma` | ✅ |
| 7 | Fourier inversion theorem | DFT→IDFT round-trip test | — | ✅ |
| 8 | Modulation (frequency shift) | — | (implicit in DTFT) | ✅ |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Function | Complexity | Status |
|---|----------|---------|-----------|--------|
| 1 | Cooley-Tukey radix-2 DIT | `fft_execute_dit()` | O(N·log N) | ✅ |
| 2 | Cooley-Tukey radix-2 DIF | `fft_execute_dif()` | O(N·log N) | ✅ |
| 3 | Split-radix FFT | `fft_execute_split_radix()` | O(N·log N) | ✅ |
| 4 | Real-valued FFT (half-length) | `fft_real_forward()` | ~½N·log(N/2) | ✅ |
| 5 | Goertzel single-bin | `goertzel_bin()` | O(N) per bin | ✅ |
| 6 | Bluestein (Chirp-Z) DFT | `bluestein_dft()` | O(N·log N) | ✅ |
| 7 | Zoom FFT | `zoom_fft()` | O(L·log L) | ✅ |
| 8 | Pruned FFT (zero-input) | `fft_pruned_input()` | O(N·log N) | ✅ |
| 9 | Pruned FFT (zero-output) | `fft_pruned_output()` | O(N·log N) | ✅ |
| 10 | Levinson-Durbin recursion | `levinson_durbin()` | O(p²) | ✅ |
| 11 | Burg's MEM | `burg_mem()` | O(N·p) | ✅ |
| 12 | Overlap-add fast convolution | `fft_convolution_overlap_add()` | O(N·log N) | ✅ |
| 13 | Overlap-save fast convolution | `fft_convolution_overlap_save()` | O(N·log N) | ✅ |
| 14 | Hann window | `window_hann()` | O(N) | ✅ |
| 15 | Hamming window | `window_hamming()` | O(N) | ✅ |
| 16 | Blackman window | `window_blackman()` | O(N) | ✅ |
| 17 | Kaiser window (Bessel I₀) | `window_kaiser()` | O(N) | ✅ |
| 18 | Dolph-Chebyshev window | `window_dolph_chebyshev()` | O(N²) | ✅ |
| 19 | Kaiser β from attenuation | `window_kaiser_beta_from_attenuation()` | O(1) | ✅ |
| 20 | Window length recommendation | `window_recommend_length()` | O(1) | ✅ |
| 21 | Hilbert transform via FFT | `hilbert_transform()` | O(N·log N) | ✅ |
| 22 | Analytic signal construction | `analytic_signal()` | O(N·log N) | ✅ |
| 23 | Wiener deconvolution | `wiener_deconvolution()` | O(N·log N) | ✅ |
| 24 | Welch's PSD | `welch_psd()` | O(K·L·log L) | ✅ |
| 25 | Blackman-Tukey correlogram | `blackman_tukey_psd()` | O(N²) | ✅ |
| 26 | Yule-Walker AR | `yule_walker_ar()` | O(p² + N·p) | ✅ |
| 27 | AR spectrum from coefficients | `ar_spectrum()` | O(M·p) | ✅ |
| 28 | Liftering (cepstral windowing) | `cepstral_liftering()` | O(N) | ✅ |
| 29 | Orthogonal Matching Pursuit (OMP) | `omp_sparse_recovery()` | O(k·M·N + k³) | ✅ |

## L6: Canonical Problems (Complete)

| # | Problem | Example/Implementation | Status |
|---|---------|----------------------|--------|
| 1 | Spectral estimation of noisy signal | `welch_psd()`, `example_spectrum_analyzer` | ✅ |
| 2 | Radar pulse compression | `chirp_generate_lfm()`, `pulse_compression_matched_filter()`, `example_radar_pulse` | ✅ |
| 3 | Bearing fault detection | `bearing_fault_detect()`, `example_vibration_analysis` | ✅ |
| 4 | GPS signal acquisition | `gps_acquisition_fft()` | ✅ |
| 5 | OFDM PAPR analysis | `ofdm_papr()`, `ofdm_papr_ccdf()` | ✅ |
| 6 | Speech formant estimation | `speech_formants()` | ✅ |
| 7 | Power quality harmonic analysis | `power_harmonics()`, `total_harmonic_distortion()` | ✅ |
| 8 | Deconvolution (image/signal restoration) | `wiener_deconvolution()` | ✅ |
| 9 | Envelope detection | `analytic_signal()` | ✅ |
| 10 | Pitch detection via cepstrum | `real_cepstrum()` | ✅ |
| 11 | Cross-spectral coherence estimation | `magnitude_squared_coherence()` | ✅ |
| 12 | Bandwidth estimation | `spectrum_bandwidth()` | ✅ |

## L7: Applications (Complete — 5+)

| # | Application | Keywords | Function | Status |
|---|------------|----------|----------|--------|
| 1 | Audio spectrum analyzer | audio, octave, IEC 61260 | `audio_octave_spectrum()` | ✅ |
| 2 | Radar LFM pulse compression | radar, chirp, matched filter | `pulse_compression_matched_filter()` | ✅ |
| 3 | Vibration bearing diagnostics | bearing, BPFO, BPFI | `bearing_fault_detect()` | ✅ |
| 4 | GPS C/A code acquisition | GPS, Doppler, PRN | `gps_acquisition_fft()` | ✅ |
| 5 | OFDM PAPR computation | OFDM, PAPR, 5G | `ofdm_papr()`, `ofdm_papr_ccdf()` | ✅ |
| 6 | Power quality (IEEE 519 THD) | power, harmonic, IEEE 519 | `power_harmonics()` | ✅ |
| 7 | Speech formant analysis | speech, formant, cepstral | `speech_formants()` | ✅ |

## L8: Advanced Topics (Complete — 4/5)

| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | MUSIC super-resolution | `music_pseudospectrum()` (simplified) | ✅ |
| 2 | Fractional Fourier Transform | `fractional_fourier_transform()` | ✅ |
| 3 | Sparse Fourier Transform | `sparse_fourier_transform()` (peak-based) | ✅ |
| 4 | Compressive sensing (OMP) | `omp_sparse_recovery()`, `omp_dft_sensing_matrix()` | ✅ |
| 5 | Wavelet transform | — | ❌ |

## L9: Research Frontiers (Partial — documented only)

| # | Topic | Status |
|---|-------|--------|
| 1 | Quantum Fourier Transform | Documented in Lean |
| 2 | 6G RIS channel estimation | Not implemented |
| 3 | AI-based super-resolution | Not implemented |
| 4 | Semantic communication | Not implemented |

---

**Knowledge Coverage Score**: L1=2, L2=2, L3=2, L4=2, L5=2, L6=2, L7=2, L8=1.5, L9=1 = **16.5/18 → COMPLETE**
