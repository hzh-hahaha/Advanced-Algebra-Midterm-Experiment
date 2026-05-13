# 研究日志

## 2026-05-13

### 已完成

- 确定课题：病态算子构造与 Tikhonov 正则化。
- 完成研究计划：[research_plan.md](research_plan.md)。
- 搭建可复现实验项目结构。
- 使用 `uv` 创建虚拟环境并安装依赖。
- 实现核心模块：
  - `src/operators.py`：病态矩阵、测试信号、相对噪声；
  - `src/solvers.py`：伪逆、Tikhonov、截断 SVD、谱滤波因子；
  - `src/metrics.py`：相对误差、残差、误差放大倍数；
  - `src/plotting.py`：理论图、误差图、信号恢复图、图像网格图。
- 完成并运行 4 个实验：
  - `experiments.exp_filter_factors`
  - `experiments.exp_1d_condition_number`
  - `experiments.exp_1d_lambda_selection`
  - `experiments.exp_2d_image_recovery`

### 当前关键结果

一维实验设置：

- $n=256$
- 相对噪声水平 $\rho=10^{-3}$
- 条件数 $\kappa(A)\in\{10^2,10^4,10^6,10^8\}$
- 每组重复 20 次

一维条件数实验均值：

| $\kappa(A)$ | 伪逆平均相对误差 | Tikhonov 平均相对误差 | 伪逆平均放大倍数 | Tikhonov 平均放大倍数 |
|---:|---:|---:|---:|---:|
| $10^2$ | 0.0115 | 0.0116 | 11.50 | 11.61 |
| $10^4$ | 0.5507 | 0.4460 | 550.71 | 445.97 |
| $10^6$ | 36.7971 | 0.6807 | 36797.09 | 680.65 |
| $10^8$ | 2741.0088 | 0.7791 | 2741008.84 | 779.15 |

解释：

- 当条件数较小时，伪逆和 Tikhonov 差别不大。
- 当条件数增大到 $10^6$ 以上时，伪逆恢复误差灾难性增长。
- Tikhonov 正则化虽然存在偏差，但显著抑制了小奇异值方向上的噪声放大。

一维 $\lambda$ 扫描：

- 条件数：$\kappa(A)=10^8$
- 噪声水平：$\rho=10^{-3}$
- 扫描范围：$\lambda\in[10^{-16},1]$
- oracle 最优 $\lambda\approx2.96\times10^{-8}$
- oracle 相对误差约为 0.7563

二维图像实验：

- 图像尺寸：$64\times64$
- 一维退化矩阵条件数：$\kappa(B)=10^3$
- 向量化二维算子条件数：$\kappa(B\otimes B)=10^6$
- 噪声水平：$\rho=10^{-3}$

二维图像恢复结果：

| 方法 | $\lambda$ | Frobenius 相对误差 | PSNR | SSIM |
|---|---:|---:|---:|---:|
| 伪逆 | - | 5.1718 | 5.5508 | 0.0228 |
| Tikhonov | $10^{-10}$ | 1.0865 | 8.2722 | 0.0807 |
| Tikhonov | $10^{-7}$ | 0.5332 | 11.2899 | 0.2111 |
| Tikhonov | $10^{-5}$ | 0.7526 | 8.0632 | 0.1053 |
| Tikhonov | $10^{-3}$ | 0.9015 | 6.3246 | 0.0472 |

解释：

- 伪逆恢复在图像上表现为严重噪声放大。
- $\lambda=10^{-7}$ 在当前网格中取得最小 Frobenius 相对误差。
- 过小的 $\lambda$ 正则化不足，过大的 $\lambda$ 抑制过强。

### 已生成图表

- `figures/theory/fig_filter_factors.png`
- `figures/signal/fig_1d_pinv_error_vs_kappa.png`
- `figures/signal/fig_1d_tikhonov_error_vs_kappa.png`
- `figures/signal/fig_1d_signal_recovery_kappa_1e8.png`
- `figures/signal/fig_1d_lambda_selection_and_lcurve.png`
- `figures/image/fig_2d_image_recovery_main.png`
- `figures/image/fig_2d_lambda_grid.png`

### 已记录数据

- `results/exp_filter_factors.json`
- `results/exp_1d_condition_number.csv`
- `results/exp_1d_condition_number.json`
- `results/exp_1d_lambda_selection.csv`
- `results/exp_1d_lambda_selection.json`
- `results/exp_2d_image_recovery.csv`
- `results/exp_2d_image_recovery.json`

### 当前判断

目前实验已经形成闭环：理论公式、数值实验、图像展示和参数选择都有对应证据。下一步应开始撰写报告正文，优先完成以下章节：

1. 奇异值分解与病态性；
2. 伪逆误差放大机制；
3. Tikhonov 正则化的谱滤波解释；
4. 一维实验结果分析；
5. 二维图像实验结果分析。
