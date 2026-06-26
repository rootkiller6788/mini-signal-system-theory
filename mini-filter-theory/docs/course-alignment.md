# Course Alignment — mini-filter-theory

## Nine-School Curriculum Mapping

| School | Course | Filter Theory Coverage | Alignment |
|--------|--------|----------------------|-----------|
| **MIT** | 6.003 Signal Processing | FIR/IIR, Butterworth, Chebyshev, bilinear | Ch. 6-7, Oppenheim |
| **Stanford** | EE102A Signal Processing | Filter specs, frequency response, window design | Lecture 14-18 |
| **Berkeley** | EE123 Digital Signal Processing | Digital filter structures, quantization, CIC | Ch. 4-6, Lyons |
| **Illinois** | ECE 310 DSP | FIR design, IIR design, stability | Ch. 5-7 |
| **Michigan** | EECS 351 DSP | Analog prototypes, bilinear, frequency transforms | Lab 4-5 |
| **Georgia Tech** | ECE 4270 DSP | Multi-rate, CIC, half-band filters | Ch. 10 |
| **TU Munich** | Signal Processing | German engineering: analog filter synthesis | Zverev 1967 |
| **ETH** | 227-0427 Signal Processing | Swiss precision: elliptic, Bessel, Cauer | Van Valkenburg 1982 |
| **Tsinghua** | Digital Signal Processing | Chinese EE: comprehensive filter design | Proakis Ch. 8 |

## Textbook Mapping

| Textbook | Chapters Mapped | Key Concepts |
|----------|----------------|--------------|
| Oppenheim & Willsky (1997) | Ch. 6, 9 | CT frequency response, Laplace, filter design |
| Oppenheim & Schafer (2010) | Ch. 5-7 | DT structures, FIR/IIR design, bilinear |
| Proakis & Manolakis (2007) | Ch. 8, 10 | FIR window design, IIR analog-based, optimal |
| Lyons (2010) | Ch. 5-7 | Practical FIR/IIR, CIC, half-band |
| Van Valkenburg (1982) | Ch. 5-12 | Butterworth, Chebyshev, elliptic, Bessel |
| Zverev (1967) | Ch. 3-8 | LC ladder prototypes, frequency transforms |
| Parks & Burrus (1987) | Ch. 3-5 | Parks-McClellan, LS design, windows |

## Implementation Coverage by Course Topic

| Topic | C Implementation | Lean Formalization | Test Coverage |
|-------|-----------------|-------------------|---------------|
| Butterworth | butterworth_poles(), _prototype() | Pole LHP theorem | ✓ |
| Chebyshev I | chebyshev1_poles(), _prototype() | - | ✓ |
| Chebyshev II | chebyshev2_pz(), _prototype() | - | ✓ |
| Elliptic | elliptic_k(), _estimate_order() | - | ✓ |
| Bessel | bessel_poles(), _prototype() | - | ✓ |
| Bilinear transform | bilinear_transform() | Stability preservation | ✓ |
| Routh-Hurwitz | tf_continuous_is_stable() | Sign change theorem | ✓ |
| Jury test | tf_discrete_is_stable() | Biquad stability condition | ✓ |
| Window functions | window_generate() (10 types) | - | ✓ |
| FIR window design | fir_design_by_window() | - | ✓ |
| IIR bilinear design | iir_design_bilinear() | - | ✓ |
| CIC filter | cic_magnitude(), _compensation() | - | ✓ |
| Half-band | fir_design_halfband() | - | ✓ |
| Frequency response | fir/iir_freq_response() | Parseval theorem | ✓ |
| Group delay | group_delay_at_fir/iir() | Linear phase condition | ✓ |
| Step response | fir_step_response() | - | ✓ |
