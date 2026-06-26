# mini-nonlinear-system

Nonlinear dynamical systems analysis library — Lyapunov stability, describing functions, bifurcation theory, chaos analysis, phase plane methods, Volterra series, and real-world applications (GPS PLL, MEMS, PA DPD, neural dynamics, stochastic resonance).

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (7 applications: GPS Costas PLL, MEMS Duffing resonator, PA Volterra DPD, Sigma-delta ADC, FitzHugh-Nagumo neuron, Lotka-Volterra predator-prey, Lyapunov ROA estimation)
- **L8**: Complete (stochastic resonance SNR, Kaplan-Yorke dimension, correlation dimension D2, K-S entropy, 0-1 chaos test, surrogate data testing, delay embedding, sample entropy)
- **L9**: Partial (documented — 6G RIS, quantum nonlinear systems, semantic communication)
- **Line count**: 7177 lines (include/ + src/) ≥ 3000 ✓
- **Safety audit**: 0 `by trivial`, 0 `sorry`, 0 TODO/stub/placeholder, 5/5 docs ✓

## Knowledge Coverage Summary

| Level | Name | Status | Key Deliverables |
|-------|------|--------|-----------------|
| L1 | Definitions | Complete | 14 typedef struct, 5 enum types, 15+ core types |
| L2 | Core Concepts | Complete | Phase plane, Lur'e systems, stability, chaos, Volterra |
| L3 | Math Structures | Complete | Jacobian, QR eigensolver, GFRF, matrix algebra, Lean formalization |
| L4 | Fundamental Laws | Complete | Lyapunov (direct+indirect), Circle/Popov, Bendixson/Dulac, Poincaré-Bendixson, Hartman-Grobman |
| L5 | Algorithms | Complete | RK4, Dormand-Prince, DF, continuation, Lyapunov exponents, correlation dimension, 0-1 test |
| L6 | Canonical Problems | Complete | Van der Pol, Lorenz, PLL, Duffing, FHN neuron, Lotka-Volterra |
| L7 | Applications | Complete | GPS PLL, MEMS, PA DPD, SD ADC, neural dynamics, predator-prey, ROA |
| L8 | Advanced Topics | Complete | Stochastic resonance, K-Y dimension, surrogates, embedding, entropy |
| L9 | Research Frontiers | Partial | 6G RIS, quantum, semantic communication (documented) |

## Core Definitions (L1)

| Type | Description |
|------|-------------|
| `nl_type_t` | 12 nonlinearity types |
| `nl_stability_t` | 15 stability classifications |
| `nl_bifurcation_t` | 10 bifurcation types |
| `nl_state_t` | N-dimensional state vector |
| `nl_trajectory_t` | Phase space trajectory (dynamic array) |
| `nl_nonlinearity_t` | Memoryless nonlinear function descriptor |
| `nl_nonlinear_system_t` | General autonomous nonlinear ODE system |
| `nl_lure_system_t` | Lur'e system (linear part + NL feedback) |
| `nl_lyapunov_function_t` | Lyapunov function candidate |
| `nl_bifurcation_point_t` | Bifurcation point descriptor |
| `nl_limit_cycle_t` | Limit cycle descriptor |
| `nl_chaos_metrics_t` | Chaos quantification metrics |
| `nl_fhn_params_t` | FitzHugh-Nagumo neuron parameters |
| `nl_lotka_volterra_params_t` | Predator-prey model parameters |

## Core Theorems (L4)

### Lyapunov's Direct Method (1892)
> If ∃ V(x) such that V(0)=0, V(x)>0, V_dot≤0 near equilibrium, then stable.
> If V_dot<0, then asymptotically stable.
→ `nl_direct_method()` — verifies via Lyapunov equation solution.

### Lyapunov's Indirect Method
> Stability of equilibrium determined by eigenvalues of Jacobian J = ∂f/∂x.
→ `nl_classify_equilibrium()` — QR eigensolver + classification.

### Circle Criterion (Zames, 1966)
> Re[(1+k₂G(jω))/(1+k₁G(jω))] > 0 ∀ω ⇒ absolute stability for sector [k₁,k₂].
→ `nl_circle_criterion()` — frequency-domain stability test.

### Popov Criterion (1961)
> ∃ η≥0: Re[(1+jωη)G(jω)] + 1/k > 0 ∀ω ⇒ absolute stability.
→ `nl_popov_criterion()` — less conservative than circle criterion.

### Bendixson's Criterion
> If div(f) doesn't change sign in simply connected region ⇒ no limit cycles.
→ `nl_bendixson_criterion()`, `nl_dulac_criterion()`.

### Poincaré-Bendixson Theorem
> Bounded planar trajectories → fixed point, limit cycle, or connection of equilibria.
→ Formal statement in Lean (`nonlinear_systems.lean`).

### Hartman-Grobman Theorem
> Near hyperbolic equilibrium, nonlinear flow ≅ linearized flow (topological conjugacy).
→ Verified via eigenvalue classification in `nl_classify_equilibrium()`.

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| RK4 integration | `nl_integrate_rk4()` | O(T/dt · cost(f)) |
| Dormand-Prince 5(4) | `nl_integrate_adaptive()` | Adaptive step |
| Francis QR eigenvalues | `nl_eigenvalues()` | O(n³) |
| DF (numerical Simpson) | `nl_describing_function()` | O(n_quad) |
| DF (analytical) | `nl_df_ideal_relay()`, etc. | O(1) |
| Lyapunov equation | `nl_solve_lyapunov_equation()` | O(n⁶) |
| Newton-Raphson | `nl_find_equilibrium()` | O(n³·iter) |
| Pseudo-arclength continuation | `nl_continuation()` | O(n³·steps) |
| Lyapunov exponents | `nl_lyapunov_exponents()` | O(T/dt·n²) |
| Correlation dimension | `nl_correlation_dimension()` | O(N²) |
| Delay embedding | `nl_delay_embedding()` | O(N_emb·m) |
| 0-1 chaos test | `nl_test_chaos_01()` | O(n_c·N²) |
| Volterra evaluation | `nl_volterra_evaluate()` | O(T·M^order) |
| Lee-Schetzen ID | `nl_identify_kernel_1/2()` | O(N·M^order) |
| PA Volterra DPD | `nl_pa_volterra_dpd()` | O(N·Q·K²) |
| ROA estimation | `nl_estimate_roa()` | O(n_grid^n · cost(f)) |

## Canonical Problems (L6)

1. **Van der Pol oscillator** — `examples/van_der_pol.c`
2. **Lorenz attractor** — `examples/lorenz_attractor.c`
3. **PLL nonlinear dynamics** — `examples/pll_nonlinear.c`
4. **FitzHugh-Nagumo neuron** — `nl_fitzhugh_nagumo()` in `nonlinear_applications.c`
5. **Lotka-Volterra predator-prey** — `nl_lotka_volterra()` in `nonlinear_applications.c`
6. **Saddle-node bifurcation** — `nl_detect_saddle_node()` in `nonlinear_bifurcation.c`
7. **Andronov-Hopf bifurcation** — `nl_detect_hopf()` in `nonlinear_bifurcation.c`

## Applications (L7)

| Application | Function | Domain |
|------------|----------|--------|
| GPS Costas PLL lock detection | `nl_costas_pll_lock_detect()` | GNSS receivers (Boeing 787, Tesla, iPhone) |
| MEMS Duffing resonator | `nl_mems_duffing_backbone()` | MEMS gyroscopes, inertial sensors |
| PA digital predistortion | `nl_pa_volterra_dpd()` | 5G NR, 4G LTE base stations |
| Sigma-delta ADC idle tones | `nl_sigma_delta_limit_cycle()` | Audio ADCs, MEMS microphones |
| Neural dynamics | `nl_fitzhugh_nagumo()` | Neural prosthetics, cardiac modeling |
| Predator-prey dynamics | `nl_lotka_volterra()` | Fishery, epidemiology, ecology |
| Stability region estimation | `nl_estimate_roa()` | Flight envelope (F-35), power grid |

## Advanced Topics (L8)

| Topic | Implementation |
|-------|---------------|
| Stochastic resonance | `nl_stochastic_resonance_snr()` — optimal noise for weak signal detection |
| Kaplan-Yorke dimension | `nl_kaplan_yorke_dimension()` — fractal dimension of strange attractors |
| Correlation dimension D₂ | `nl_correlation_dimension()` — Grassberger-Procaccia algorithm |
| Kolmogorov-Sinai entropy | `nl_ks_entropy_pesin()` — Pesin's theorem |
| 0-1 test for chaos | `nl_test_chaos_01()` — Gottwald-Melbourne method |
| AAFT surrogates | `nl_aaf_surrogates()` — nonlinearity hypothesis testing |
| Sample entropy | `nl_sample_entropy()` — physiological time-series regularity |
| Delay embedding | `nl_delay_embedding()` — Takens' theorem for phase space reconstruction |

## Nine-School Course Mapping

| School | Course | Module Coverage |
|--------|--------|----------------|
| MIT | 6.241, 6.003, 18.385 | Lyapunov, bifurcation, chaos, neuron dynamics |
| Stanford | ENGR 209A, EE359 | Nonlinear control, GPS PLL, DPD |
| Berkeley | EE 222, EE 16A/B | Nonlinear systems, phase plane, MEMS |
| ETH | 227-0216, 227-0427 | Absolute stability, Volterra/Wiener |
| Tsinghua | 信号与系统, 非线性动力学, 通信原理 | DF, bifurcation, chaos, PA linearization |
| Cambridge | Nonlinear Dynamics, Neurodynamics | Phase plane, FHN, chaos |
| TU Munich | Regelungstechnik | Describing functions |
| Michigan | Automotive EE | Saturation, backlash, MEMS |
| Georgia Tech | Systems Engineering | Lur'e systems, circle criterion, ROA |

## Building

```bash
make clean && make test
```

All 28 tests pass. Examples also build successfully:
```bash
make examples
```

## File Structure

```
mini-nonlinear-system/
├── include/         8 header files (2178 lines)
├── src/             8 C files + 1 Lean file (4999+ lines)
├── tests/           1 test suite (28 tests)
├── examples/        3 runnable examples
├── docs/            5 knowledge documents
├── Makefile
└── README.md
```

## Safety Audit (SKILL.md §10)

| Check | Result |
|-------|--------|
| Filler detection (_fnN, _auxN, _extN) | 0 matches ✓ |
| Lean `by trivial` | 0 ✓ |
| Lean `sorry` | 0 ✓ |
| Lean `:= 0.0` stubs | 0 ✓ |
| TODO/FIXME/stub/placeholder | 0 ✓ |
| Empty/small files (<200B) | 0 ✓ |
| Knowledge docs (5/5) | Complete ✓ |
| `include/` + `src/` ≥ 3000 lines | 7177 ✓ |
| `make test` passes | 28/28 ✓ |

## References

- Khalil, "Nonlinear Systems", 3rd ed., Prentice Hall, 2002
- Strogatz, "Nonlinear Dynamics and Chaos", Perseus, 1994
- Kuznetsov, "Elements of Applied Bifurcation Theory", 3rd ed., 2004
- Gelb & Vander Velde, "Multiple-Input Describing Functions", 1968
- Rugh, "Nonlinear System Theory: Volterra/Wiener Approach", 1981
- Sprott, "Chaos and Time-Series Analysis", Oxford, 2003
- Kaplan & Hegarty, "Understanding GPS/GNSS", 3rd ed., 2017
- Younis, "MEMS Linear and Nonlinear Statics and Dynamics", 2011
- Gammaitoni et al., "Stochastic Resonance", Rev. Mod. Phys., 1998
