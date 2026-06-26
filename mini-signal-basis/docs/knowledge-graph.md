# Knowledge Graph - mini-signal-basis
## L1: Definitions (Complete)
signal_value_t, signal_time_t, signal_freq_t, signal_phase_t, signal_index_t - type aliases
signal_domain_t - CT vs DT
signal_periodicity_t - periodic / aperiodic
signal_causality_t - causal / anticausal / noncausal
signal_energy_class_t - energy / power / neither
signal_randomness_t - deterministic / stochastic
signal_symmetry_t - even / odd / none
signal_ct_t - continuous-time signal struct (samples, t_start, t_end, dt, num_samples, domain)
signal_dt_t - discrete-time signal struct (samples, n_start, n_end, length, domain)
Elementary signal param structs: sinusoid, exponential, cisoid, rect, tri, gaussian, sinc
Signal properties: energy, power, RMS, peak, mean, variance, SNR

## L2: Core Concepts (Complete)
Memory management: alloc(ct/dt), free(ct/dt), copy(ct/dt)
10 elementary signal generators (step, delta, ramp, sinusoid, cisoid, exponential, rect, tri, gaussian, sinc)
Time-domain operations: shift, reverse, scale, circular shift
Amplitude operations: gain, offset
Signal arithmetic: add, subtract, multiply, divide
Even/odd decomposition
Inner product, norm, distance
Zero-crossing rate, envelope detection
Numerical derivative and integral, interpolation, decimation
Random signal generation: PRNG, uniform, gaussian, AWGN, colored noise
