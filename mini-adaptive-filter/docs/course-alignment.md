# Course Alignment: mini-adaptive-filter

## MIT
- **6.003 Signals and Systems**: Gradient descent, convergence analysis
- **6.450 Digital Communications**: Adaptive equalization, LMS/RLS for channel estimation
- **6.630 Electromagnetic Waves**: Adaptive beamforming foundations

## Stanford
- **EE264 Adaptive Filtering**: Complete LMS/RLS/Kalman coverage (full course mapping)
- **EE359 Wireless Communications**: Adaptive equalization for multipath channels
- **EE363 Linear Dynamical Systems**: Kalman filter theory

## Berkeley
- **EE225D Adaptive Signal Processing**: LMS/RLS algorithms, convergence analysis
- **EE123 Digital Signal Processing**: Adaptive filter structure and applications
- **EE221A Linear Systems Theory**: Optimal estimation, Kalman filtering

## Illinois (UIUC)
- **ECE 310 Digital Signal Processing**: Adaptive filtering basics (LMS)
- **ECE 459 Communications**: Adaptive equalization for digital comms

## Michigan
- **EECS 351 Digital Signal Processing**: LMS algorithm
- **EECS 455 Communications**: Adaptive channel estimation
- **EECS 411 Microwave Circuits**: Adaptive impedance matching

## Georgia Tech
- **ECE 4270 Digital Signal Processing**: LMS/RLS algorithm comparison
- **ECE 6601 Communications**: Adaptive receiver structures

## TU Munich
- **Signal Processing**: Kalman filter, adaptive algorithms
- **Communications**: Adaptive equalization, OFDM channel estimation

## ETH Zurich
- **227-0427 Signal Processing**: LMS/RLS theory and applications
- **227-0436 Communications**: Adaptive filtering for wireless receivers

## Tsinghua (清华)
- **信号与系统**: Gradient descent, Wiener filtering
- **通信原理**: Adaptive equalization
- **数字信号处理**: LMS/RLS algorithm implementation

## Chapter Mapping (Haykin, "Adaptive Filter Theory" 5th ed.)
| Chapter | Topic | Implementation |
|---------|-------|---------------|
| Ch 2 | Wiener Filters | af_wiener_solve(), Levinson-Durbin |
| Ch 4 | Method of Steepest Descent | Gradient computation in LMS |
| Ch 5 | LMS Algorithm | af_lms_direct_update() |
| Ch 6 | Normalized LMS | af_lms_normalized_update() |
| Ch 7 | Transform-Domain LMS | Block LMS gradient |
| Ch 9 | RLS Algorithm | af_rls_standard_update() |
| Ch 10 | QR-RLS | af_rls_qr_givens_update() |
| Ch 11 | Robust Adaptive Filters | Sign LMS variants |
| Ch 14 | Kalman Filters | af_kalman_predict/update |
| Ch 16 | Applications | Adaptive noise/echo cancellation |
