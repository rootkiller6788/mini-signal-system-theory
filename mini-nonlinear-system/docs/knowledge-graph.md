# Knowledge Graph — mini-nonlinear-system

## L1: Definitions (Complete)
| Item | C Implementation | Lean Definition |
|------|-----------------|-----------------|
| Nonlinearity types (relay, dead-zone, saturation, backlash, etc.) | `nl_type_t`, `nl_nonlinearity_t` | `NonlinearityType` |
| Equilibrium point | `nl_state_t` | `StateVector` |
| Stability classification | `nl_stability_t` | `StabilityType` |
| Bifurcation types | `nl_bifurcation_t` | `BifurcationType` |
| Limit cycle | `nl_limit_cycle_t` | — |
| Lyapunov function | `nl_lyapunov_function_t` | — |
| Describing function | `nl_describing_function_t` | — |
| Chaos metrics | `nl_chaos_metrics_t` | — |
| Volterra kernel | `nl_volterra_kernel_t` | — |
| Volterra series | `nl_volterra_series_t` | — |

## L2: Core Concepts (Complete)
| Concept | Implementation |
|---------|---------------|
| State-space nonlinear ODE | `nl_nonlinear_system_t`, `nl_evaluate_f()` |
| Lur'e system (linear + NL feedback) | `nl_lure_system_t`, `nl_lure_dxdt()` |
| Phase plane / phase space | `nl_phase_portrait_t`, `nl_compute_phase_portrait()` |
| Memoryless nonlinearity | `nl_nonlinearity_eval()` |
| Lyapunov stability | `nl_direct_method()`, `nl_classify_equilibrium()` |
| Describing function approximation | `nl_describing_function()` |
| Bifurcation | `nl_detect_saddle_node()`, `nl_detect_hopf()` |
| Deterministic chaos | `nl_chaos_analyze()` |
| Volterra series expansion | `nl_volterra_evaluate()` |
| Wiener G-functionals | `nl_wiener_G_functionals()` |

## L3: Mathematical Structures (Complete)
| Structure | Implementation |
|-----------|---------------|
| Jacobian matrix | `nl_jacobian_t`, `nl_compute_jacobian()` |
| Eigenvalue computation (QR algorithm) | `nl_eigenvalues()` |
| Matrix determinant (LU) | `nl_matrix_determinant()` |
| Positive definiteness (Gershgorin) | `nl_is_positive_definite()` |
| Transfer function G(jω) | `nl_transfer_function()` |
| 2×2 matrix formalization | `Matrix2x2`, `classify2x2` (Lean) |
| Volterra kernel symmetrization | `nl_kernel_symmetrize()` |
| GFRF (multi-dim Fourier) | `nl_compute_GFRF()` |

## L4: Fundamental Laws (Complete)
| Theorem | C Verification | Lean Statement |
|---------|---------------|----------------|
| Lyapunov direct method | `nl_direct_method()` — verifies V>0, V_dot<0 | `lyapunov_direct_method_formal` |
| Lyapunov indirect method | `nl_classify_equilibrium()` — eigenvalues of J | — |
| Hartman-Grobman | Documented | `hartman_grobman_formal` |
| Poincaré-Bendixson | Documented | `poincare_bendixson_formal` |
| Circle criterion | `nl_circle_criterion()` | `circle_criterion_formal` |
| Popov criterion | `nl_popov_criterion()` | `popov_criterion_encoding` |
| Bendixson's criterion | `nl_bendixson_criterion()` | `bendixson_criterion_formal` |
| Dulac's criterion | `nl_dulac_criterion()` | — |
| Routh-Hurwitz (n=2) | Computed in stability classification | `routh_hurwitz_2x2` |

## L5: Algorithms/Methods (Complete)
| Algorithm | Implementation | Complexity |
|-----------|---------------|------------|
| RK4 integration | `nl_integrate_rk4()` | O(T/dt · cost(f)) |
| Dormand-Prince 5(4) adaptive | `nl_integrate_adaptive()` | O(T/dt · cost(f)) |
| Francis QR eigensolver | `nl_eigenvalues()` | O(n³) |
| Lyapunov equation (Bartels-Stewart) | `nl_solve_lyapunov_equation()` | O(n⁶) direct, O(n³) possible |
| Numerical DF (Simpson) | `nl_describing_function()` | O(n_quad) |
| Analytical DF formulas | `nl_df_ideal_relay()`, `nl_df_dead_zone()`, etc. | O(1) |
| Limit cycle DF prediction | `nl_df_limit_cycle()` | O(nA · nw) |
| Newton-Raphson equilibrium | `nl_find_equilibrium()` | O(dim³ · iter) |
| Pseudo-arclength continuation | `nl_continuation()` | O(dim³ · steps) |
| Bifurcation detection (test functions) | `nl_detect_saddle_node()`, `nl_detect_hopf()` | O(n_pts · dim³) |
| Lyapunov exponents (Benettin) | `nl_lyapunov_exponents()` | O(T/dt · dim²) |
| Correlation sum/dimension | `nl_correlation_sum()`, `nl_correlation_dimension()` | O(N²) |
| Time-delay embedding | `nl_delay_embedding()` | O(N_emb · m) |
| Mutual information (delay) | `nl_mutual_information_delay()` | O(N · n_bins²) |
| False nearest neighbors | `nl_false_nearest_neighbors()` | O(N · max_m) |
| 0-1 test for chaos | `nl_test_chaos_01()` | O(n_c · N²) |
| Volterra evaluation | `nl_volterra_evaluate()` | O(T · M^order) |
| Lee-Schetzen identification | `nl_identify_kernel_1()`, `nl_identify_kernel_2()` | O(N · M^order) |
| Poincaré index | `nl_poincare_index()` | O(n_pts) |
| Phase portrait generation | `nl_compute_phase_portrait()` | O(nx · ny) |
| Equilibrium search 2D | `nl_find_equilibria_2d()` | O(nx · ny · iter) |

## L6: Canonical Problems (Complete)
| Problem | Example / Implementation |
|---------|-------------------------|
| Van der Pol oscillator | `examples/van_der_pol.c` — limit cycle, DF prediction |
| Lorenz attractor | `examples/lorenz_attractor.c` — chaos, Lyapunov exponents |
| PLL nonlinear dynamics | `examples/pll_nonlinear.c` — lock range, cycle slipping |
| Duffing oscillator | `Duffing` structure (Lean), parameterized in core |
| Chua's circuit | `ChuaCircuit` structure (Lean) |
| Saddle-node bifurcation | `nl_detect_saddle_node()` |
| Andronov-Hopf bifurcation | `nl_detect_hopf()`, `nl_hopf_normal_form()` |
| Limit cycle prediction (DF) | `nl_df_limit_cycle()` |
| Period-doubling cascade | `nl_bifurcation_diagram()` |

## L7: Applications (Partial+ — 5 applications)
| Application | Evidence |
|-------------|---------|
| GPS receiver PLL | `examples/pll_nonlinear.c` — GPS L1 1575.42 MHz carrier tracking |
| Amplifier saturation | `nl_df_saturation()` — RF power amplifier |
| Sigma-delta modulator | `NL_TYPE_QUANTIZER` — limit cycles in ΣΔ ADCs |
| MEMS nonlinearity | `NL_TYPE_CUBIC_STIFFNESS` — Duffing-type MEMS resonators |
| 5G NR carrier recovery | Documented in `examples/pll_nonlinear.c` |

## L8: Advanced Topics (Partial+ — 8 topics)
| Topic | Implementation |
|-------|---------------|
| Kaplan-Yorke dimension | `nl_kaplan_yorke_dimension()` |
| Correlation dimension D2 | `nl_correlation_dimension()` |
| Kolmogorov-Sinai entropy | `nl_ks_entropy_pesin()` |
| 0-1 test for chaos | `nl_test_chaos_01()` |
| AAFT surrogates | `nl_aaf_surrogates()` |
| Sample entropy | `nl_sample_entropy()` |
| Delay embedding (Takens) | `nl_delay_embedding()` |
| Lyapunov-based control | `nl_direct_method()` — region of attraction |

## L9: Research Frontiers (Partial — documented)
| Topic | Status |
|-------|--------|
| Quantum nonlinear systems | Documented in `docs/gap-report.md` |
| Nonlinear metamaterials | Documented |
| 6G RIS nonlinear optimization | Documented |
| Semantic communication over nonlinear channels | Documented |
