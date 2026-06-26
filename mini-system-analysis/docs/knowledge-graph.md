# Knowledge Graph - mini-system-analysis

## L1: Definitions (Complete - 14 items)
13 core struct/enum types defined in system_defs.h:

| # | Definition | Type | File |
|---|-----------|------|------|
| D01 | ct_signal_t | Continuous-time signal with sampling info | system_defs.h |
| D02 | dt_signal_t | Discrete-time signal | system_defs.h |
| D03 | dt_complex_signal_t | Complex-valued DT signal | system_defs.h |
| D04 | system_type_t | LTI/LTV/NL_TI/NL_TV/HYBRID/STOCHASTIC enum | system_defs.h |
| D05 | system_property_t | 14-bit property mask (causal, stable, FIR, etc.) | system_defs.h |
| D06 | ct_impulse_response_t | CT h(t) with energy, peak, settling time | system_defs.h |
| D07 | dt_impulse_response_t | DT h[n] with energy, abs_sum | system_defs.h |
| D08 | ct_step_response_t | CT s(t) with rise time, overshoot, delay | system_defs.h |
| D09 | ct_transfer_function_t | H(s) = N(s)/D(s) rational function | system_defs.h |
| D10 | dt_transfer_function_t | H(z) = N(z)/D(z) rational function | system_defs.h |
| D11 | pole_zero_t | Single pole/zero with damping, wn | system_defs.h |
| D12 | pole_zero_collection_t | Complete PZ map with stability flag | system_defs.h |
| D13 | frequency_response_t | Bode data (magnitude, phase, group delay) | system_defs.h |
| D14 | system_order_info_t | Order, type number, relative degree | system_defs.h |

## L2: Core Concepts (Complete - 16 items)

| # | Concept | Implementation | File |
|---|---------|---------------|------|
| C01 | CT convolution integral (trapezoidal) | ct_convolution_direct | convolution.c |
| C02 | DT linear convolution O(N^2) | dt_convolution_direct | convolution.c |
| C03 | DT circular convolution | dt_circular_convolution_direct | convolution.c |
| C04 | FFT-based convolution (conv. theorem) | dt_convolution_fft | convolution.c |
| C05 | Overlap-add block convolution | dt_overlap_add | convolution.c |
| C06 | Overlap-save block convolution | dt_overlap_save | convolution.c |
| C07 | Frequency-domain deconvolution | dt_deconvolution_freq | convolution.c |
| C08 | Cross-correlation R_xy | dt_cross_correlation | convolution.c |
| C09 | Auto-correlation R_xx | dt_auto_correlation | convolution.c |
| C10 | Normalized cross-correlation | dt_cross_correlation_normalized | convolution.c |
| C11 | Step/impulse relationship | dt_step_from_impulse / dt_impulse_from_step | convolution.c |
| C12 | Convolution commutativity check | dt_convolution_commutativity_check | convolution.c |
| C13 | Magnitude response |H(jw)| compute_magnitude_response | frequency_response.c |
| C14 | Phase response arg(H(jw)) | compute_phase_response | frequency_response.c |
| C15 | -3 dB bandwidth detection | find_3db_bandwidth / find_3db_cutoff_frequency | frequency_response.c |
| C16 | Filter type classification | is_lowpass/highpass/bandpass/bandstop | frequency_response.c |

## L3: Math Structures (Complete - 27 items)

| # | Structure | Implementation | File |
|---|-----------|---------------|------|
| M01 | CT transfer function evaluation H(s) | ct_tf_evaluate (Horner) | transfer_function.c |
| M02 | DT transfer function evaluation H(z) | dt_tf_evaluate (Horner) | transfer_function.c |
| M03 | Vectorized TF evaluation | ct_tf_evaluate_vector / dt_tf_evaluate_vector | transfer_function.c |
| M04 | CT frequency response H(jw) sweep | ct_frequency_response | transfer_function.c |
| M05 | DT frequency response H(e^{jw}) sweep | dt_frequency_response | transfer_function.c |
| M06 | Pole finding (Newton+deflation) | ct_tf_find_poles | transfer_function.c |
| M07 | Zero finding | ct_tf_find_zeros | transfer_function.c |
| M08 | Full pole-zero analysis | ct_tf_pole_zero_analysis | transfer_function.c |
| M09 | Partial fraction expansion (residues) | ct_tf_partial_fraction | transfer_function.c |
| M10 | Bilinear (Tustin) transform | ct_to_dt_bilinear | transfer_function.c |
| M11 | Impulse invariance discretization | ct_to_dt_impulse_invariance | transfer_function.c |
| M12 | Zero-order hold (ZOH) discretization | ct_to_dt_zoh | transfer_function.c |
| M13 | Matched Z-transform | ct_to_dt_matched_z | transfer_function.c |
| M14 | System order analysis | ct_tf_order_analysis | transfer_function.c |
| M15 | ZPK gain extraction | ct_tf_zpk_gain | transfer_function.c |
| M16 | Bode plot generation | ct_tf_bode_plot | transfer_function.c |
| M17 | Nyquist plot generation | ct_tf_nyquist_plot | transfer_function.c |
| M18 | Group delay computation | ct_tf_group_delay | transfer_function.c |
| M19 | CT state-space struct (A,B,C,D) | ct_ss_alloc | state_space.c |
| M20 | Controllability matrix (Kalman) | ct_ss_controllability_matrix | state_space.c |
| M21 | Observability matrix | ct_ss_observability_matrix | state_space.c |
| M22 | Matrix exponential e^{At} (Pade) | ct_ss_matrix_exponential | state_space.c |
| M23 | CT simulation (Euler/RK4) | ct_ss_simulate_euler / ct_ss_simulate_rk4 | state_space.c |
| M24 | Similarity transform | ct_ss_similarity_transform | state_space.c |
| M25 | Canonical forms (CCF/OCF) | ct_ss_to_controllable/observable_canonical | state_space.c |
| M26 | Lyapunov equation A'P+PA+Q=0 | ct_ss_lyapunov | state_space.c |
| M27 | TF from state-space | ct_ss_to_transfer_function | state_space.c |

## L4: Fundamental Laws (Complete - 14 items)

| # | Theorem/Law | C Implementation | Lean Formalization |
|---|------------|-----------------|-------------------|
| T01 | BIBO stability (CT: LHP, DT: unit circle) | ct_is_bibo_stable, dt_is_bibo_stable | ctBIBOStable, dtBIBOStable |
| T02 | Routh-Hurwitz criterion | routh_hurwitz | countSignChanges |
| T03 | Jury stability test | jury_stability_test | juryF conditions |
| T04 | Nyquist criterion (N=Z-P) | nyquist_stability | -- (encirclement counting) |
| T05 | Cayley-Hamilton (A^2 - tr(A)A + det(A)I = 0) | -- | cayley_hamilton_2x2 |
| T06 | Lyapunov direct method (A'P+PA=-Q) | ct_lyapunov_stability | isPositiveDefinite2x2 |
| T07 | Settling time t_s = 4/(zeta*wn) | ct_settling_time_estimate | settlingTime2pct |
| T08 | Peak overshoot M_p = exp(-pi*zeta/sqrt(1-zeta^2)) | ct_overshoot_estimate | peakOvershootPct |
| T09 | Convolution theorem | dt_convolution_fft | convSum |
| T10 | Kalman controllability | ct_ss_is_controllable | isControllable2x2 |
| T11 | Kalman observability | ct_ss_is_observable | -- |
| T12 | Gain/Phase margin criterion | compute_stability_margins | -- |
| T13 | Feedback interconnection stability | feedback_stability_check | -- |
| T14 | Root locus (Evans) | root_locus | -- |

## L5: Algorithms (Complete - 9 items)

| # | Algorithm | Complexity | File |
|---|-----------|-----------|------|
| A01 | Horner polynomial evaluation | O(n) | transfer_function.c |
| A02 | Newton-Raphson root finding | O(n*iter) per root | transfer_function.c |
| A03 | Synthetic polynomial deflation | O(n) | transfer_function.c |
| A04 | Faddeev-LeVerrier (char. poly) | O(n^4) | state_space.c |
| A05 | Gauss-Jordan elimination (pivoting) | O(n^3) | state_space.c |
| A06 | Pade approximation (matrix exp) | O(n^3) | state_space.c |
| A07 | Runge-Kutta 4 (RK4) | O(n^2) per step | state_space.c |
| A08 | DFT/IDFT direct O(N^2) | O(N^2) | convolution.c |
| A09 | Matrix rank via Gaussian elim. | O(n^3) | state_space.c |

## L6: Canonical Problems (Complete - 4 items)

| # | Problem | Example/Demo |
|---|---------|-------------|
| P01 | Step response of 1st/2nd-order systems | ex1_step_response.c |
| P02 | Stability analysis (Routh-Hurwitz, margins) | ex2_stability.c |
| P03 | State-space simulation (mass-spring-damper) | ex3_state_space.c |
| P04 | Bode plot analysis and filter characterization | demo_bode_analysis.c |

## L7: Applications (Complete - 5 items)

| # | Application | Implementation | Keyword |
|---|------------|---------------|---------|
| A01 | GPS position-velocity tracking (8-state ECEF) | kalman_filter.c (kf_configure_gps_tracker + kf_gps_update) | GPS |
| A02 | Sensor fusion (information filter) | kalman_filter.c (kf_to_information_form) | sensor fusion |
| A03 | Control system stability margin analysis | ex2_stability.c (GM/PM) | --
| A04 | Mass-spring-damper vibration analysis | ex3_state_space.c (RK4) | mass-spring-damper |
| A05 | 60 Hz notch filter (power line filtering) | demo_bode_analysis.c | 60 Hz |

## L8: Advanced Topics (Complete - 5 items)

| # | Topic | Implementation |
|---|-------|---------------|
| A01 | Stochastic state estimation (Kalman filter) | kalman_filter.c (predict/update/NIS/log-likelihood) |
| A02 | Extended Kalman filter (EKF) for nonlinear systems | kalman_filter.c (ekf_predict/ekf_update) |
| A03 | Rauch-Tung-Striebel (RTS) fixed-interval smoother | kalman_filter.c (kf_rts_smooth) |
| A04 | Information filter (canonical Y=P⁻¹ form) | kalman_filter.c (kf_to_information_form) |
| A05 | Lyapunov stability direct method with P solution | stability.c + state_space.c + system_analysis.lean |

## L9: Research Frontiers (Partial - documented)

| # | Topic | Status |
|---|-------|--------|
| R01 | Compressed sensing for system identification | Documented in gap-report.md |
| R02 | 6G RIS (Reconfigurable Intelligent Surfaces) | Documented in gap-report.md |
| R03 | Quantum control and measurement | Documented in gap-report.md |
