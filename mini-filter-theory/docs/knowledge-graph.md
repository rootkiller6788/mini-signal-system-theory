# Knowledge Graph — mini-filter-theory

## L1: Definitions (Complete)
- filter_type_t: 10 canonical filter types (LP, HP, BP, BS, AP, Notch, Peaking, LowShelf, HighShelf, Custom)
- approx_type_t: 8 approximation families (Butterworth, Chebyshev I/II, Elliptic, Bessel, Gaussian, Legendre, Synchronous)
- filter_topology_t: 8 analog circuit topologies (Sallen-Key, MFB, State-Variable, Biquad, GIC, FDNR, Gyrator, Ladder)
- digital_structure_t: 8 digital structures (DF-I, DF-II, Transposed I/II, Cascade, Parallel, Lattice, Wave Ladder)
- fir_linear_phase_type_t: 4 linear phase types
- filter_spec_t: Complete design specification
- biquad_section_t: Digital SOS
- analog_biquad_t: Analog SOS
- pz_pair_t / pz_map_t: Pole-zero representation
- fir_filter_t / iir_filter_t: Complete filter descriptors
- filter_state_t: Streaming state
- freq_resp_point_t / freq_resp_t: Frequency response
- window_type_t: 10 window functions
- filter_error_t: 12 error codes

## L2: Core Concepts (Complete)
- Frequency-selective filtering (LP/HP/BP/BS)
- Direct Form I vs II (canonical vs non-canonical)
- Transposed forms (zero-input limit cycle behavior)
- Cascade vs Parallel SOS (sensitivity vs noise)
- Lattice structures (orthogonal, adaptive)
- Linear phase symmetry conditions (4 types)
- Frequency transformations (LP→LP, LP→HP, LP→BP, LP→BS)
- BIBO stability concept
- Impulse response vs frequency response duality
- Group delay and phase delay
- Passband ripple, stopband attenuation, transition band
- Filter order and roll-off rate (20N dB/decade for Butterworth)

## L3: Mathematical Structures (Complete)
- Polynomial transfer functions (continuous and discrete)
- Horner's method for complex evaluation
- Polynomial convolution (cascade multiplication)
- Partial fraction expansion (PFE)
- Bessel polynomials (recurrence and explicit formula)
- Binomial coefficients and combinatorial expansion
- Complete elliptic integral K(k) via AGM
- Modified Bessel function I0(x) (series + asymptotic)
- Durand-Kerner simultaneous root finding
- Gaussian elimination with partial pivoting
- Numerical integration (trapezoidal, Simpson)
- Root finding (bisection, Newton-Raphson)

## L4: Fundamental Laws (Complete)
- Routh-Hurwitz stability criterion (continuous-time)
- Jury stability test (discrete-time)
- Bilinear transform stability preservation
- Parseval's theorem (energy conservation)
- Linear phase condition (symmetry ↔ constant group delay)
- Biquad stability condition (|a2|<1, |a1|<1+a2)
- Kaiser window optimality (DPSS, maximum energy concentration)

## L5: Algorithms/Methods (Complete)
- Butterworth pole computation (trigonometric)
- Chebyshev I pole computation (hyperbolic functions)
- Chebyshev II pole/zero computation (reciprocal transformation)
- Elliptic K(k) via AGM method
- Bessel pole computation (Durand-Kerner on Bessel polynomials)
- Window function generation (10 types)
- FIR window design method
- FIR frequency sampling design
- IIR bilinear transform design
- Bilinear pre-warping
- Kaiser beta optimization formula
- FIR order estimation (Kaiser/Bellanger)
- IIR order estimation (Butterworth, Chebyshev, Elliptic, Bessel)
- Coefficient quantization simulation

## L6: Canonical Problems (Complete)
- Design a Butterworth lowpass filter from specification
- Design an FIR filter by windowing with performance comparison
- Audio crossover design (Linkwitz-Riley)
- Bessel pulse preservation comparison
- FIR differentiator design
- FIR Hilbert transformer design
- Half-band filter design
- CIC filter design and compensation
- Frequency response analysis (magnitude, phase, group delay)
- Step response analysis (overshoot, settling time, rise time)
- SNR improvement through filtering
- THD analysis after filtering

## L7: Applications (Partial+)
- Audio crossover networks (JAES, Linkwitz 1976)
- ECG/medical waveform preservation (Bessel filters)
- Oscilloscope front-end anti-aliasing
- DAC reconstruction filtering (sinc compensation)
- Communication channel filtering

## L8: Advanced Topics (Partial)
- Coefficient sensitivity analysis (L2 scaling)
- CIC compensation filter design
- Least-squares FIR design
- Spectral factorization
- Digital filter quantization effects

## L9: Research Frontiers (Partial)
- AI-designed filters (documented)
- Quantum filtering (documented)
- Fractional-order filters (documented)
