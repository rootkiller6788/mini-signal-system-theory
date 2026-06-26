# Mini Signal & System Theory（迷你信号与系统理论）

**从零开始、零依赖的 C 语言实现**，涵盖大学级信号处理与系统理论。每个模块对应 MIT、Stanford、Berkeley 等顶尖大学的课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-adaptive-filter](mini-adaptive-filter/) | LMS、RLS、Kalman 滤波、自适应 FIR/IIR/Lattice/Subband 结构 | MIT 6.341, Stanford EE264 |
| [mini-filter-theory](mini-filter-theory/) | Butterworth、Chebyshev I/II、Elliptic、Bessel 模拟/数字滤波器设计 | MIT 6.003, Stanford EE102A, Berkeley EE105 |
| [mini-fourier-analysis](mini-fourier-analysis/) | Cooley-Tukey FFT、DFT、卷积、相关、谱估计、窗函数 | MIT 6.003, Stanford EE261 |
| [mini-laplace-z-transform](mini-laplace-z-transform/) | Laplace 变换对、Z 变换、有理函数、s 域到 z 域映射 | MIT 6.003, Stanford EE102A |
| [mini-nonlinear-system](mini-nonlinear-system/) | 分岔分析、混沌、Lyapunov 稳定性、描述函数、MEMS、GPS PLL | MIT 6.450, Stanford EE359, Berkeley EE222 |
| [mini-signal-basis](mini-signal-basis/) | 连续/离散信号类型、基本信号、信号分解、DTMF、ECG QRS 检测 | Stanford EE359, Michigan EECS 455 |
| [mini-system-analysis](mini-system-analysis/) | 连续/离散卷积、频率响应、稳定性分析、冲激/阶跃响应 | MIT 6.003, Stanford EE102A |
| [mini-system-identification](mini-system-identification/) | LS/WLS/PEM/ML/IV 估计、RLS、LMS、NLMS、Kalman 滤波、EKF | MIT 6.435, Stanford EE263 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明
- **实用演示程序** — 自适应噪声消除、滤波器设计器、频谱分析仪、心电检测器、系统辨识器等

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-adaptive-filter
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-signal-system-theory/
├── mini-adaptive-filter/       # 自适应滤波：LMS、RLS、Kalman、格型结构
├── mini-filter-theory/         # 模拟/数字滤波器设计：Butterworth、Chebyshev、Elliptic
├── mini-fourier-analysis/      # 傅里叶分析：FFT、DFT、谱估计、窗函数
├── mini-laplace-z-transform/   # Laplace 与 Z 变换：变换对、s-to-z 映射
├── mini-nonlinear-system/      # 非线性系统：分岔、混沌、Lyapunov 稳定性
├── mini-signal-basis/          # 信号基础：信号类型、分解、实际应用
├── mini-system-analysis/       # LTI 系统分析：卷积、频率响应、稳定性
└── mini-system-identification/ # 系统辨识：LS、PEM、ML、RLS、Kalman 参数估计
```

## 许可证

MIT
