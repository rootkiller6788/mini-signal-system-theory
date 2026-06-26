# Course Alignment — System Identification

## Nine-School Curriculum Mapping

### MIT — 6.435 System Identification (Ljung)
| Topic | Implementation |
|-------|---------------|
| LS estimation | `sysid_ls_solve()` |
| PEM (prediction error) | `sysid_pem_armax/oe/bj()` |
| Model structure selection | `sysid_order_selection_aicc()` |
| Validation and residual analysis | `sysid_ljung_box_q()`, cross-correlation tests |

### Stanford — EE264 Digital Signal Processing / EE392M System ID
| Topic | Implementation |
|-------|---------------|
| Adaptive filtering (LMS, RLS) | `sysid_lms_step()`, `sysid_rls_step()` |
| Recursive identification | `sysid_rpem_armax_step()` |
| Kalman filtering | `sysid_kalman_step()` |

### Berkeley — EECS 225D Audio Signal Processing / EE C220B Experience
| Topic | Implementation |
|-------|---------------|
| Spectral estimation | `sysid_blackman_tukey()` |
| Frequency-domain methods | `sysid_freq_ls_fit()`, `sysid_etfe()` |
| Subspace methods | `sysid_n4sid()`, `sysid_moesp()` |

### Illinois — ECE 544 Topics in Signal Processing
| Topic | Implementation |
|-------|---------------|
| System modeling fundamentals | Core model structures (ARX, ARMAX, etc.) |
| Real-time implementation | RLS, LMS adaptive methods |

### Michigan — EECS 551 Matrix Methods in Signal Processing
| Topic | Implementation |
|-------|---------------|
| SVD for identification | `sysid_svd_power()` in subspace methods |
| QR decomposition | `sysid_qr_decomp()` |
| Cholesky factorization | `sysid_cholesky()`, `sysid_cholesky_solve()` |

### Georgia Tech — ECE 6250 Advanced DSP
| Topic | Implementation |
|-------|---------------|
| Adaptive filter theory | Complete LMS/NLMS/Sign-LMS family |
| System ID applications | DC motor, RC circuit, mechanical system examples |

### TU Munich — System Identification
| Topic | Implementation |
|-------|---------------|
| Industrial process ID | Hammerstein/Wiener models |
| Nonlinear identification | NARX polynomial/sigmoid networks |
| Frequency-domain ID | SK iteration, LPM |

### ETH Zurich — 227-0427 Signal Processing & System Identification
| Topic | Implementation |
|-------|---------------|
| Subspace identification | Complete N4SID/MOESP/CVA |
| Closed-loop identification | SSARX |
| Model validation | Full validation pipeline |

### Tsinghua (清华) — 信号与系统 / 系统辨识
| Topic | Implementation |
|-------|---------------|
| LS estimation | OLS, WLS, TLS, RLS |
| Model structure analysis | Pole computation, stability, DC gain |
| Signal generation for ID | PRBS, chirp, multisine |

## Unique Institutional Perspectives Captured
- MIT: PEM paradigm and model validation rigor
- Stanford: Adaptive/recursive filter focus
- ETH: Numerical linear algebra for subspace methods
- TUM: Industrial nonlinear modeling
- Georgia Tech: Complete DSP+ID integration
- Berkeley: Spectral and frequency-domain methods
