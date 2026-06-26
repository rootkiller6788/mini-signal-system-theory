/**
 * @file    fourier_fft.h
 * @brief   Fast Fourier Transform Algorithms — L5 Algorithms + L4 Laws
 *
 * @details Provides industrial-grade FFT implementations: Cooley-Tukey
 *          radix-2 decimation-in-time (DIT), radix-2 decimation-in-frequency
 *          (DIF), split-radix, real-valued optimized FFT (RFFT), and the
 *          inverse transforms.  Also includes the Goertzel algorithm for
 *          single-bin evaluation and the Bluestein (chirp Z) algorithm for
 *          arbitrary-length DFTs.
 *
 * Knowledge Mapping:
 *   L5 - Algorithms:
 *     - Cooley-Tukey radix-2 DIT FFT (1965)
 *     - Cooley-Tukey radix-2 DIF FFT
 *     - Split-radix FFT (Duhamel & Hollmann 1984, min flops for radix-2)
 *     - Real-valued FFT (packing two N-point real FFTs into one N-point complex FFT)
 *     - Goertzel algorithm (single-bin DFT, O(N) per bin)
 *     - Bluestein / Chirp-Z transform (arbitrary-length DFT via convolution)
 *     - Pruned FFT (zero-input / zero-output pruning)
 *   L4 - Theorems:
 *     - DFT ↔ circular convolution duality
 *     - Convolution theorem: FFT-based fast convolution
 *     - Circular shift property of DFT
 *   L6 - Problems:
 *     - Fast convolution / fast correlation via FFT
 *     - Spectral leakage analysis and mitigation
 *
 * Reference: Cooley & Tukey, "An algorithm for the machine calculation
 *            of complex Fourier series", Math. Comp. 19 (1965), pp.297-301.
 *            Oppenheim & Schafer, Ch.9.
 */

#ifndef FOURIER_FFT_H
#define FOURIER_FFT_H

#include "fourier_core.h"

/* ─── L5: FFT Plan Lifecycle ────────────────────────────────────────────── */

/**
 * @brief  Create an FFT plan for Cooley-Tukey radix-2 DIT.
 *
 *         Precomputes twiddle factors W_N^k = exp(-j·2π·k/N) and
 *         bit-reversal permutation table.  The plan can be reused
 *         across multiple transforms of the same size.
 *
 * @param  N           Transform length (must be power of 2, N ≥ 2)
 * @param  is_inverse  0 for forward (exp(-j)), 1 for inverse (exp(+j))
 * @return             Initialized FFT plan; NULL if N is not power-of-2
 *
 * Complexity: O(N) precomputation time, O(N) memory for twiddle + bitrev
 */
fft_plan_t*  fft_plan_create(size_t N, int is_inverse);

/**
 * @brief  Destroy an FFT plan and free all associated memory.
 */
void  fft_plan_destroy(fft_plan_t *plan);

/* ─── L5: Core FFT Transforms ───────────────────────────────────────────── */

/**
 * @brief  Execute Cooley-Tukey radix-2 DIT FFT in-place.
 *
 *         Algorithm:
 *           1. Bit-reversal permutation of input array
 *           2. log₂(N) stages, each stage combines pairs via butterfly:
 *              X'[k]    = X[k] + W_N^r · X[k + step]
 *              X'[k+step] = X[k] - W_N^r · X[k + step]
 *
 *         This is the foundational O(N·log N) DFT algorithm.
 *
 * @param  plan  FFT plan with precomputed twiddle + bitrev tables
 * @param  data  Complex input/output array (length plan->N), transformed in-place
 *
 * Complexity: O(N·log N) time, O(1) additional space
 * Theorem:     This computes exactly the DFT (to within floating-point precision)
 */
void  fft_execute_dit(fft_plan_t *plan, fourier_complex_t *data);

/**
 * @brief  Execute Cooley-Tukey radix-2 DIF FFT in-place.
 *
 *         Decimation-in-frequency applies the butterfly before the
 *         recursive subdivision.  Natural order input, bit-reversed output.
 *
 * @param  plan  FFT plan
 * @param  data  Complex input/output array (length plan->N)
 */
void  fft_execute_dif(fft_plan_t *plan, fourier_complex_t *data);

/**
 * @brief  Execute split-radix FFT (in-place).
 *
 *         The split-radix algorithm decomposes the DFT into one radix-2
 *         and two radix-4 sub-transforms, achieving the lowest known
 *         arithmetic complexity for radix-2 lengths:
 *           - Real multiplications: (N/3)·log₂N (approximately)
 *           - Real additions: N·log₂N
 *
 * @param  plan  FFT plan (must have been initialized with compat. twiddles)
 * @param  data  Complex array of length N (power of 2)
 */
void  fft_execute_split_radix(fft_plan_t *plan, fourier_complex_t *data);

/**
 * @brief  Compute forward FFT of a real-valued signal, returning only
 *         the non-redundant half (0 ≤ k ≤ N/2).
 *
 *         Uses the packing trick: two N-point real sequences x1, x2
 *         are combined into one N-point complex sequence, then one
 *         complex FFT yields both spectra via symmetry unfolding.
 *
 *         For a single real sequence, we split even/odd samples into
 *         two half-length real sequences, reducing cost by ~2×.
 *
 * @param  x         Real-valued input (length N, N must be even)
 * @param  N         Number of samples
 * @param  fs        Sampling rate (Hz)
 * @param  spectrum  Output: complex spectrum X[0..N/2] (length N/2+1)
 * @return           dft_result_t with real-only half-spectrum
 *
 * Complexity: ~N/2 · log₂(N/2) complex ops (≈ ½ of full complex FFT)
 */
dft_result_t  fft_real_forward(const double *x, size_t N, double fs);

/**
 * @brief  Compute inverse FFT to reconstruct real signal from half-spectrum.
 *
 *         Inverse of fft_real_forward.  Takes the non-redundant half
 *         (0..N/2), reconstructs the full conjugate-symmetric spectrum,
 *         then performs inverse FFT to obtain real-valued time domain.
 *
 * @param  X_half  Half-spectrum X[0..N/2] (length N/2+1)
 * @param  N       Full transform length (must be even)
 * @param  x_out   Output real signal (length N, caller-allocated)
 */
void  fft_real_inverse(const fourier_complex_t *X_half, size_t N, double *x_out);

/* ─── L5: Goertzel Algorithm ────────────────────────────────────────────── */

/**
 * @brief  Goertzel algorithm: compute a single DFT bin X[k] in O(N) time.
 *
 *         Derived from the filter-bank interpretation of the DFT:
 *         each bin k is the output of a resonator with pole at z = e^{j·2π·k/N}
 *         evaluated at n = N.
 *
 *         Goertzel recurrence (second-order, real arithmetic for real input):
 *           s[n] = x[n] + 2·cos(2π·k/N)·s[n-1] - s[n-2]
 *           X[k] = s[N] - e^{-j·2π·k/N}·s[N-1]
 *
 *         Advantage over full FFT: O(N) per bin vs O(N·log N) for all bins.
 *         Useful for DTMF tone detection, where only 8 bins are needed.
 *
 * @param  x   Input signal (length N)
 * @param  N   Signal length
 * @param  k   DFT bin index to evaluate (0 ≤ k < N)
 * @return     Complex DFT coefficient X[k]
 *
 * Complexity: O(N) time, O(1) space per bin
 */
fourier_complex_t  goertzel_bin(const double *x, size_t N, size_t k);

/**
 * @brief  Goertzel algorithm for multiple bins.
 *
 * @param  x      Input signal
 * @param  N      Signal length
 * @param  bins   Array of bin indices to evaluate
 * @param  nbins  Number of bins
 * @param  X_out  Output complex coefficients (length nbins, caller-allocated)
 */
void  goertzel_multi_bin(const double *x, size_t N,
                         const size_t *bins, size_t nbins,
                         fourier_complex_t *X_out);

/* ─── L6: FFT-Based Fast Convolution ────────────────────────────────────── */

/**
 * @brief  Fast convolution via FFT (overlap-add method).
 *
 *         Computes y = x * h where x is a long signal and h is a
 *         filter kernel.  Uses overlap-add with FFT block size L.
 *
 *         Theorem (Convolution):  DFT{x * h} = DFT{x} · DFT{h}
 *
 *         Algorithm:
 *           1. Segment x into blocks of length L
 *           2. Zero-pad each block and h to length N = L + M - 1 (power of 2)
 *           3. Compute X = FFT{x_block}, H = FFT{h}
 *           4. Y = X · H  (element-wise multiply)
 *           5. y_block = IFFT{Y}
 *           6. Overlap-add consecutive blocks
 *
 * @param  x         Long input signal
 * @param  x_len     Length of x
 * @param  h         Filter kernel (FIR)
 * @param  h_len     Length of h
 * @param  block_len Processing block length
 * @param  y_out     Output signal (length x_len + h_len - 1, caller-allocated)
 *
 * Complexity: O((x_len + h_len)·log(block_len)) via FFT
 */
void  fft_convolution_overlap_add(const double *x, size_t x_len,
                                   const double *h, size_t h_len,
                                   size_t block_len, double *y_out);

/**
 * @brief  Fast convolution via FFT (overlap-save method).
 *
 *         Alternative to overlap-add, discarding contaminated samples
 *         from circular convolution wraparound.
 *
 * @param  x         Long input signal
 * @param  x_len     Length of x
 * @param  h         Filter kernel
 * @param  h_len     Length of h
 * @param  block_len Processing block length (FFT size)
 * @param  y_out     Output signal (length x_len + h_len - 1)
 */
void  fft_convolution_overlap_save(const double *x, size_t x_len,
                                    const double *h, size_t h_len,
                                    size_t block_len, double *y_out);

/* ─── L5: Bluestein / Chirp-Z Transform ─────────────────────────────────── */

/**
 * @brief  Bluestein's algorithm: compute DFT of arbitrary length N via
 *         convolution (not restricted to power-of-2).
 *
 *         Reformulation:  X[k] = W_N^{-k²/2} · Σ_{n=0}^{N-1} (x[n]·W_N^{-n²/2}) · W_N^{(k-n)²/2}
 *
 *         The sum is a convolution, computable via FFT in O(N·log N).
 *         This is the foundation of the Chirp-Z Transform.
 *
 * @param  x   Input signal (length N, any N ≥ 1)
 * @param  N   Arbitrary length
 * @param  fs  Sampling rate
 * @return     DFT result
 */
dft_result_t  bluestein_dft(const double *x, size_t N, double fs);

/* ─── L6: Zoom FFT ──────────────────────────────────────────────────────── */

/**
 * @brief  Zoom FFT: compute high-resolution DFT over a narrow frequency band.
 *
 *         Uses the chirp Z-transform with:
 *           A = e^{j·2π·f_start/fs}  (starting frequency)
 *           W = e^{-j·2π·(f_end - f_start)/(M·fs)}  (fine step)
 *
 *         Enables millihertz resolution in a narrow band without computing
 *         a full N-point FFT at the same resolution (which would require
 *         N = fs/Δf points, potentially millions).
 *
 * @param  x          Input signal
 * @param  N          Signal length
 * @param  fs         Sampling rate (Hz)
 * @param  f_start    Start frequency (Hz)
 * @param  f_end      End frequency (Hz)
 * @param  M          Number of output bins
 * @param  X_zoom     Output: zoom spectrum (length M, caller-allocated)
 * @param  freq_axis  Output: frequency grid (length M, caller-allocated)
 */
void  zoom_fft(const double *x, size_t N, double fs,
               double f_start, double f_end, size_t M,
               fourier_complex_t *X_zoom, double *freq_axis);

/* ─── L5: Pruned FFT ────────────────────────────────────────────────────── */

/**
 * @brief  Zero-input pruned FFT: compute DFT when only the first L samples
 *         are non-zero (L ≪ N).  Prunes unnecessary butterfly operations.
 *
 * @param  x          Input signal (first L samples non-zero, rest zero)
 * @param  L          Number of non-zero input samples
 * @param  N          Total FFT length (power of 2)
 * @param  X_out      Output spectrum (length N, caller-allocated)
 */
void  fft_pruned_input(const double *x, size_t L, size_t N,
                       fourier_complex_t *X_out);

/**
 * @brief  Zero-output pruned FFT: compute only the first M DFT bins
 *         (M ≪ N).  Prunes butterflies that don't contribute to output.
 *
 * @param  x          Input signal (length N)
 * @param  N          FFT length (power of 2)
 * @param  M          Number of output bins to compute
 * @param  X_out      Output: first M DFT bins (length M)
 */
void  fft_pruned_output(const double *x, size_t N, size_t M,
                        fourier_complex_t *X_out);

#endif /* FOURIER_FFT_H */
