# 病态线性算子的误差放大与 Tikhonov 正则化研究计划

## 1. 课题定位

本项目选择课程课题 **“病态算子构造与正则化”**，研究病态线性系统中微小扰动如何通过伪逆求解被显著放大，并进一步分析 Tikhonov 正则化如何通过奇异值方向上的谱滤波稳定求解。

建议正式题目：

> **从奇异值分解看病态线性系统的误差放大与 Tikhonov 正则化稳定机制**

本题与高等代数联系紧密，核心内容包括：

- 线性映射与矩阵表示；
- 正交矩阵与正交基变换；
- 奇异值分解；
- 矩阵条件数；
- Moore-Penrose 伪逆；
- 正规方程；
- 最小二乘问题；
- 二次型与正定矩阵；
- 谱分解视角下的滤波机制。

本项目不把图像复原作为单纯应用，而是把图像或一维信号作为向量数据，用它们可视化病态线性算子的数值不稳定性。项目的数学主线是：

> 病态矩阵的小奇异值导致伪逆中的巨大放大因子；Tikhonov 正则化通过改变奇异值方向上的放大系数，抑制噪声在小奇异值方向上的爆炸。

## 2. 联网调研依据与科研规范来源

本计划参考了以下科研计算与可重复研究的权威实践：

1. Wilson 等人的 *Good enough practices in scientific computing* 指出，科研计算项目应当像实验室实验一样组织数据、代码、文档和计算步骤，避免数据丢失、步骤不可追踪和结果难以复现。  
   来源：<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510>

2. Wilkinson 等人的 FAIR 原则强调科研数据和数字对象应尽可能做到 Findable、Accessible、Interoperable、Reusable，即可发现、可访问、可互操作、可复用。  
   来源：<https://www.nature.com/articles/sdata201618>

3. 美国国家科学院报告 *Reproducibility and Replicability in Science* 将计算可重复性定义为在相同输入数据、计算步骤、方法、代码和分析条件下得到一致计算结果。  
   来源：<https://www.ncbi.nlm.nih.gov/books/NBK547535/>

4. NIH 数据管理与共享政策强调数据管理、元数据记录、共享和透明性有助于提升研究严谨性和可重复性。  
   来源：<https://sharing.nih.gov/data-management-and-sharing-policy/sharing-scientific-data>

这些原则在本项目中落实为：

- 所有实验脚本化，图表由脚本自动生成；
- 固定随机种子；
- 保存软件环境和依赖版本；
- 区分原始数据、处理后数据和图表输出；
- 记录每次实验的矩阵规模、条件数、噪声水平、正则化参数；
- 对随机噪声实验重复多次，报告均值和标准差；
- 明确区分理论结论、数值现象和展示性案例；
- 让报告中的每一张图都能被代码复现。

## 3. 核心研究问题

本项目围绕以下问题展开：

1. 如何构造一个条件数极大的病态线性算子？
2. 为什么病态算子会使微小观测扰动造成巨大恢复误差？
3. 伪逆求解中的误差放大如何由奇异值分解精确解释？
4. Tikhonov 正则化为什么能稳定病态问题？
5. 正则化参数 $\lambda$ 如何影响恢复误差？
6. 如何用数值实验验证理论推导？
7. 一维信号实验与二维图像实验是否表现出一致的误差放大规律？

## 4. 数学模型

设真实信号或图像向量为：

$$
x\in\mathbb{R}^n
$$

线性观测模型为：

$$
b = Ax
$$

其中：

$$
A\in\mathbb{R}^{m\times n}
$$

是人为构造的病态线性算子。实际观测带有扰动：

$$
\tilde b = b + \delta b = Ax + \delta b
$$

目标是根据 $\tilde b$ 恢复 $x$。

如果直接使用伪逆：

$$
\hat x_{\mathrm{pinv}} = A^+\tilde b
$$

则恢复误差为：

$$
\hat x_{\mathrm{pinv}} - x = A^+\delta b
$$

如果 $A$ 病态，$A^+\delta b$ 可能远大于 $\delta b$。

Tikhonov 正则化求解：

$$
\hat x_\lambda = \arg\min_x \left(\|Ax-\tilde b\|_2^2+\lambda\|x\|_2^2\right)
$$

其正规方程为：

$$
(A^TA+\lambda I)\hat x_\lambda = A^T\tilde b
$$

当 $\lambda>0$ 时，若 $A^TA+\lambda I$ 正定，则闭式解为：

$$
\hat x_\lambda = (A^TA+\lambda I)^{-1}A^T\tilde b
$$

## 5. 奇异值分解视角

令 $A$ 的奇异值分解为：

$$
A = U\Sigma V^T
$$

其中：

$$
U^TU = I,\quad V^TV = I
$$

$$
\Sigma = \mathrm{diag}(\sigma_1,\sigma_2,\dots,\sigma_r)
$$

且：

$$
\sigma_1\geq\sigma_2\geq\cdots\geq\sigma_r>0
$$

矩阵条件数定义为：

$$
\kappa(A)=\frac{\sigma_1}{\sigma_r}
$$

若 $\kappa(A)$ 很大，则 $A$ 是病态矩阵。

### 5.1 伪逆误差放大

伪逆为：

$$
A^+ = V\Sigma^+U^T
$$

其中：

$$
\Sigma^+ = \mathrm{diag}\left( \frac1{\sigma_1}, \frac1{\sigma_2}, \dots, \frac1{\sigma_r} \right)
$$

因此噪声项满足：

$$
A^+\delta b = \sum_{i=1}^r \frac{u_i^T\delta b}{\sigma_i}v_i
$$

这说明噪声在左奇异向量 $u_i$ 方向上的分量 $u_i^T\delta b$，会被 $\frac1{\sigma_i}$ 放大到右奇异向量 $v_i$ 方向上。

如果某些 $\sigma_i$ 极小，则即使 $u_i^T\delta b$ 很小，也可能造成巨大的恢复误差。

### 5.2 Tikhonov 的谱滤波解释

由：

$$
\hat x_\lambda = (A^TA+\lambda I)^{-1}A^T\tilde b
$$

结合 $A=U\Sigma V^T$，可得：

$$
\hat x_\lambda = \sum_{i=1}^r \frac{\sigma_i}{\sigma_i^2+\lambda} (u_i^T\tilde b)v_i
$$

伪逆解中每个奇异值方向上的放大因子是：

$$
\frac1{\sigma_i}
$$

Tikhonov 解中对应的滤波因子是：

$$
g_\lambda(\sigma_i) = \frac{\sigma_i}{\sigma_i^2+\lambda}
$$

当 $\sigma_i$ 很小时：

$$
g_\lambda(\sigma_i) \approx \frac{\sigma_i}{\lambda}
$$

这避免了 $\frac1{\sigma_i}$ 在小奇异值方向上的无限放大。

因此，Tikhonov 正则化的本质不是简单“让图像变平滑”，而是在奇异向量基下对不同谱方向进行选择性抑制。

## 6. 实验总体设计

本项目包含两个主要实验：

1. 一维信号实验：用于严谨验证误差放大、条件数影响和正则化参数影响。
2. 二维图像实验：用于展示病态算子造成的视觉退化，以及正则化恢复的直观效果。

两类实验共用同一数学框架：

$$
\tilde b = Ax + \delta b
$$

并比较：

- 真实数据 $x$；
- 病态观测 $Ax$；
- 加噪观测 $\tilde b$；
- 伪逆恢复 $\hat x_{\mathrm{pinv}}$；
- Tikhonov 恢复 $\hat x_\lambda$。

## 7. 实验一：一维病态线性系统

### 7.1 数据构造

取：

$$
n=256
$$

构造真实信号：

$$
x(t) = \sin(2\pi f_1t) +0.5\sin(2\pi f_2t) +p(t)
$$

其中：

- $f_1=3$；
- $f_2=17$；
- $p(t)$ 是局部脉冲或分段突变，用来模拟非平滑结构。

这样构造的 $x$ 同时含有低频、较高频和局部结构，能检验算法对不同成分的恢复能力。

### 7.2 病态矩阵构造

构造：

$$
A=U\Sigma V^T
$$

其中 $U,V$ 是随机正交矩阵，可由随机高斯矩阵的 QR 分解得到。

设置奇异值为指数衰减：

$$
\sigma_i =10^{-\alpha\frac{i-1}{n-1}}, \quad i=1,2,\dots,n
$$

则：

$$
\sigma_1=1,\quad \sigma_n=10^{-\alpha}
$$

因此：

$$
\kappa(A)=10^\alpha
$$

建议测试：

$$
\kappa(A)\in \{10^2,10^4,10^6,10^8\}
$$

即：

$$
\alpha\in\{2,4,6,8\}
$$

### 7.3 噪声设置

设置相对噪声水平：

$$
\rho = \frac{\|\delta b\|_2}{\|Ax\|_2}
$$

建议测试：

$$
\rho\in\{10^{-4},10^{-3},10^{-2}\}
$$

噪声生成方式：

1. 生成标准正态噪声 $z\sim N(0,I)$；
2. 标准化：

$$
z\leftarrow \frac{z}{\|z\|_2}
$$

3. 缩放：

$$
\delta b = \rho \|Ax\|_2 z
$$

这样可以精确控制扰动大小。

### 7.4 对比方法

比较以下方法：

1. 伪逆解：

$$
\hat x_{\mathrm{pinv}} = A^+\tilde b
$$

2. Tikhonov 解：

$$
\hat x_\lambda = (A^TA+\lambda I)^{-1}A^T\tilde b
$$

3. 截断 SVD 解，可作为补充对照：

$$
\hat x_k = \sum_{i=1}^k \frac{u_i^T\tilde b}{\sigma_i}v_i
$$

截断 SVD 不是主线，但可以帮助解释“舍弃小奇异值方向”和 Tikhonov 连续滤波之间的关系。

### 7.5 评价指标

相对恢复误差：

$$
E_x = \frac{\|\hat x-x\|_2}{\|x\|_2}
$$

残差：

$$
R = \frac{\|A\hat x-\tilde b\|_2}{\|\tilde b\|_2}
$$

解范数：

$$
N = \|\hat x\|_2
$$

误差放大倍数：

$$
M = \frac{\|\hat x-x\|_2/\|x\|_2}{\|\delta b\|_2/\|b\|_2}
$$

该指标直接反映“输入扰动”到“输出误差”的放大。

### 7.6 重复实验

每组参数重复 20 次随机噪声实验：

$$
\text{repeat}=20
$$

记录：

- 平均相对误差；
- 标准差；
- 最大误差；
- 最小误差。

报告中不只画单次结果，应画均值曲线和误差条。

## 8. 实验二：二维图像病态复原

### 8.1 图像数据

选取一张灰度图像：

$$
X\in\mathbb{R}^{m\times m}
$$

建议：

$$
m=64
$$

如果计算性能足够，可使用：

$$
m=128
$$

图像来源可以使用：

- 自己拍摄的静物照片并裁剪成灰度图；
- `scikit-image` 中的标准测试图；
- 自行构造的几何图案。

为了可复现，建议报告主实验使用公开测试图或项目中保存的固定图片；课堂展示可补充自己拍摄图片。

### 8.2 二维线性退化模型

直接构造 $m^2\times m^2$ 的大矩阵计算成本较高。建议使用可分离线性算子：

$$
Y = BXB^T
$$

其中：

$$
B\in\mathbb{R}^{m\times m}
$$

是病态矩阵。

向量化后：

$$
\mathrm{vec}(Y) = (B\otimes B)\mathrm{vec}(X)
$$

因此整体线性算子是：

$$
A = B\otimes B
$$

利用 Kronecker 积性质：

$$
\kappa(A) = \kappa(B\otimes B) = \kappa(B)^2
$$

这是一条非常适合在课堂展示的高代结论：二维图像退化会把一维病态性进一步平方放大。

### 8.3 加噪模型

实际观测：

$$
\tilde Y = BXB^T + E
$$

其中：

$$
\frac{\|E\|_F}{\|BXB^T\|_F} = \rho
$$

建议测试：

$$
\rho\in\{10^{-4},10^{-3},10^{-2}\}
$$

### 8.4 图像恢复方法

伪逆恢复：

$$
\hat X_{\mathrm{pinv}} = B^+\tilde Y(B^+)^T
$$

Tikhonov 恢复可通过向量化表达：

$$
\hat x_\lambda = (A^TA+\lambda I)^{-1}A^T\tilde y
$$

其中：

$$
A=B\otimes B, \quad \tilde y=\mathrm{vec}(\tilde Y)
$$

如果 $m=64$，则 $A$ 为 $4096\times4096$，仍可在普通电脑上计算，但需要注意内存。

为了展示高代结构，也可以利用 $B$ 的 SVD 来实现更高效的二维谱滤波。若：

$$
B=P\Sigma Q^T
$$

则：

$$
B\otimes B = (P\otimes P)(\Sigma\otimes\Sigma)(Q\otimes Q)^T
$$

其奇异值为：

$$
\sigma_i(B)\sigma_j(B)
$$

Tikhonov 滤波因子变为：

$$
\frac{\sigma_i\sigma_j}{(\sigma_i\sigma_j)^2+\lambda}
$$

这可以作为报告中的高级亮点。

### 8.5 图像评价指标

使用 Frobenius 相对误差：

$$
E_X = \frac{\|\hat X-X\|_F}{\|X\|_F}
$$

同时使用图像指标：

- PSNR；
- SSIM。

但报告中应明确：PSNR 和 SSIM 是图像质量辅助指标，主线仍然是矩阵范数误差和奇异值谱分析。

## 9. 正则化参数 $\lambda$ 的选择

$\lambda$ 是本项目的关键变量。不能只说“调参”，必须给出合理选择策略。

### 9.1 Oracle 最优参数

实验中已知真实 $x$，因此可以定义：

$$
\lambda_{\mathrm{oracle}} = \arg\min_\lambda \frac{\|x_\lambda-x\|_2}{\|x\|_2}
$$

该参数用于评估正则化方法的最优可能效果。

注意：真实应用中通常不知道 $x$，所以 oracle 参数不能作为实际方法，只能作为理论对照。

### 9.2 L-curve 方法

对不同 $\lambda$ 计算：

$$
\|A x_\lambda-\tilde b\|_2
$$

和：

$$
\|x_\lambda\|_2
$$

画出：

$$
\log\|A x_\lambda-\tilde b\|_2 \quad \text{vs} \quad \log\|x_\lambda\|_2
$$

L-curve 的拐角处对应“残差较小”和“解不过分震荡”的折中。

### 9.3 残差原则

若噪声水平 $\|\delta b\|_2$ 已知或可估计，则选择满足：

$$
\|A x_\lambda-\tilde b\|_2 \approx \|\delta b\|_2
$$

的 $\lambda$。

该原则的直观含义是：不应把噪声也过度拟合进去。

### 9.4 参数扫描范围

建议使用对数网格：

$$
\lambda\in \{10^{-16},10^{-15},\dots,10^0\}
$$

或更细：

$$
\lambda_j=10^{a_j}, \quad a_j\in[-16,0]
$$

取 100 个点。

## 10. 需要生成的图表

最终报告和展示至少应包含以下图表。

### 10.1 理论解释类图表

1. 奇异值谱图  
   横轴为 $i$，纵轴为 $\sigma_i$，纵轴使用对数坐标。

2. 条件数变化图  
   展示不同 $\alpha$ 下 $\kappa(A)=10^\alpha$ 的变化。

3. 伪逆放大因子图  
   画：

$$
\frac1{\sigma_i}
$$

   随 $i$ 的变化。

4. Tikhonov 滤波因子图  
   对不同 $\lambda$ 画：

$$
\frac{\sigma_i}{\sigma_i^2+\lambda}
$$

5. 伪逆因子与 Tikhonov 因子对比图  
   这是展示中最重要的图之一。

### 10.2 一维实验图表

1. 原始信号、观测信号、伪逆恢复、Tikhonov 恢复对比图。
2. 相对误差随条件数变化的曲线。
3. 相对误差随噪声水平变化的曲线。
4. 相对误差随 $\lambda$ 变化的 U 型曲线。
5. L-curve 图。
6. 20 次重复实验的均值和标准差误差条。

### 10.3 图像实验图表

1. 原图 $X$。
2. 退化图 $BXB^T$。
3. 加噪观测图 $\tilde Y$。
4. 伪逆恢复图。
5. Tikhonov 恢复图。
6. 不同 $\lambda$ 下恢复图像对比。
7. 图像恢复误差、PSNR、SSIM 表格。

## 11. 预期结果与理论解释

### 11.1 预期现象

1. 当 $\kappa(A)$ 较小时，伪逆恢复可以接近真实解。
2. 当 $\kappa(A)$ 增大到 $10^6$ 或 $10^8$ 时，即使噪声水平只有 $10^{-4}$，伪逆恢复也可能严重失真。
3. Tikhonov 正则化可以显著降低恢复误差。
4. $\lambda$ 太小时，正则化不足，噪声仍被放大。
5. $\lambda$ 太大时，真实信号也被过度抑制，出现过平滑或幅值偏小。
6. 最优 $\lambda$ 位于二者之间，体现偏差-方差折中。

### 11.2 理论解释

伪逆解：

$$
\hat x_{\mathrm{pinv}} = \sum_i \frac{u_i^T\tilde b}{\sigma_i}v_i
$$

噪声项：

$$
\sum_i \frac{u_i^T\delta b}{\sigma_i}v_i
$$

当 $\sigma_i$ 很小时，噪声被放大。

Tikhonov 解：

$$
\hat x_\lambda = \sum_i \frac{\sigma_i}{\sigma_i^2+\lambda} (u_i^T\tilde b)v_i
$$

其中：

$$
\frac{\sigma_i}{\sigma_i^2+\lambda}
$$

对小奇异值方向起抑制作用。

因此，实验中看到的“伪逆恢复崩溃”和“Tikhonov 恢复稳定”不是偶然现象，而是由奇异值谱结构决定的。

## 12. 项目代码与文件组织

建议项目结构如下：

```text
.
├── docs/
│   └── research_plan.md
├── data/
│   ├── raw/
│   └── processed/
├── src/
│   ├── operators.py
│   ├── solvers.py
│   ├── metrics.py
│   ├── plotting.py
│   └── config.py
├── experiments/
│   ├── exp_1d_condition_number.py
│   ├── exp_1d_lambda_selection.py
│   ├── exp_2d_image_recovery.py
│   └── exp_filter_factors.py
├── figures/
│   ├── theory/
│   ├── signal/
│   └── image/
├── report/
│   ├── report.md
│   └── references.bib
├── slides/
└── pyproject.toml
```

### 12.1 环境管理

建议使用 `uv` 管理 Python 环境：

```powershell
uv init
uv add numpy scipy matplotlib scikit-image pillow pandas jupyter
```

如果后续制作网页 PPT，可使用 `pnpm` 管理前端依赖。

`cargo` 暂不需要，除非后续决定实现高性能矩阵实验或命令行工具。

### 12.2 随机性控制

所有随机实验统一使用固定随机种子：

```python

rng = np.random.default_rng(20260513)

```

每个重复实验可以使用：

```python

seed = 20260513 + trial_id

```

报告中记录主随机种子。

### 12.3 结果复现规范

所有图表必须由脚本生成，不手工改图。

建议命令：

```powershell

uv run python experiments/exp_filter_factors.py
uv run python experiments/exp_1d_condition_number.py
uv run python experiments/exp_1d_lambda_selection.py
uv run python experiments/exp_2d_image_recovery.py

```

每个脚本输出：

- `.png` 或 `.pdf` 图表；
- `.csv` 数值结果；
- 一份 `.json` 参数记录。

## 13. 科研最佳实践落实清单

### 13.1 数据管理

- `data/raw/` 只保存原始数据，不覆盖、不修改。
- `data/processed/` 保存处理后的数据。
- 每个数据文件应有说明：来源、日期、处理方式、尺寸、格式。
- 如果使用自己拍摄图片，应记录拍摄设备、裁剪方式和灰度化方法。

### 13.2 代码管理

- 每个实验脚本只负责一个明确实验问题。
- 公共函数放入 `src/`，避免在多个实验脚本中复制粘贴。
- 函数命名应反映数学含义，例如 `make_ill_conditioned_matrix`、`solve_tikhonov`、`relative_error`。
- 不在代码中硬编码输出路径，应统一放在配置中。

### 13.3 参数记录

每次实验记录：

- 矩阵规模 $n$ 或图像尺寸 $m$；
- 条件数 $\kappa(A)$；
- 奇异值衰减参数 $\alpha$；
- 噪声水平 $\rho$；
- 正则化参数 $\lambda$；
- 重复次数；
- 随机种子；
- 软件版本。

### 13.4 图表规范

每张图必须包含：

- 标题；
- 横轴名称；
- 纵轴名称；
- 图例；
- 单位或量纲说明；
- 对数坐标说明；
- 参数说明。

所有图表文件命名应可读，例如：

```text

fig_singular_values_kappa_1e8.png
fig_filter_factors_lambda_compare.png
fig_1d_error_vs_kappa.png
fig_image_recovery_lambda_grid.png

```

### 13.5 结果解释规范

报告中每个实验图都要回答三个问题：

1. 图中现象是什么？
2. 该现象与理论推导如何对应？
3. 该现象说明了什么结论？

避免使用没有量化依据的表述，例如：

- “效果很好”；
- “误差比较小”；
- “恢复比较明显”。

应改成：

- “相对误差从 12.4 降至 0.18”；
- “PSNR 从 8.1 dB 提升至 24.6 dB”；
- “当 $\kappa(A)=10^8$ 时，伪逆误差均值约为 Tikhonov 误差的 70 倍”。

## 14. 组内分工

### 成员 A：理论推导与报告主线

负责：

- SVD 与伪逆推导；
- 条件数与病态性解释；
- Tikhonov 正则化推导；
- 谱滤波机制说明；
- 报告第 1 至第 5 节初稿；
- 展示中数学理论部分。

核心交付：

- 理论推导文档；
- 关键公式解释；
- 课堂展示的数学主线。

### 成员 B：算法实现与数值实验

负责：

- 病态矩阵构造；
- 一维信号实验；
- Tikhonov 求解器；
- $\lambda$ 扫描；
- 重复实验统计；
- 保存 `.csv` 和参数文件。

核心交付：

- `src/operators.py`；
- `src/solvers.py`；
- `experiments/exp_1d_condition_number.py`；
- `experiments/exp_1d_lambda_selection.py`。

### 成员 C：图像实验、可视化与展示

负责：

- 图像数据准备；
- 二维图像复原实验；
- 图表美化；
- PPT 或网页展示；
- 展示讲稿；
- 最终报告排版检查。

核心交付：

- `experiments/exp_2d_image_recovery.py`；
- `figures/image/`；
- 展示文件；
- 结果图注和展示讲稿。

三名成员都必须理解以下核心公式：

$$
A^+\delta b = \sum_i \frac{u_i^T\delta b}{\sigma_i}v_i
$$

和：

$$
\hat x_\lambda = \sum_i \frac{\sigma_i}{\sigma_i^2+\lambda} (u_i^T\tilde b)v_i
$$

这两条公式是答辩时最可能被追问的内容。

## 15. 时间安排

当前日期为 2026 年 5 月 13 日，目标是在 2026 年 5 月 31 日前完成报告与展示材料。

| 日期 | 工作内容 | 产出 |
|---|---|---|
| 5.13 | 确定课题、研究问题、项目结构、分工 | 研究计划 |
| 5.14 | 完成 SVD、伪逆、条件数理论推导 | 理论初稿 |
| 5.15 | 完成 Tikhonov 正则化与谱滤波推导 | 理论公式图 |
| 5.16 | 实现病态矩阵构造和一维信号实验 | 初版代码 |
| 5.17 | 完成条件数与噪声水平实验 | 误差曲线 |
| 5.18 | 实现 $\lambda$ 扫描、oracle 参数、L-curve | 参数选择图 |
| 5.19 | 完成二维图像退化与伪逆恢复 | 图像对比图 |
| 5.20 | 完成二维 Tikhonov 恢复 | 图像恢复结果 |
| 5.21 | 重复实验 20 次，统计均值和标准差 | 统计表格 |
| 5.22 | 整理所有图表，检查坐标轴和图注 | 图表定稿 |
| 5.23 | 撰写报告主体 | 报告初稿 |
| 5.24 | 补充参考文献、实验设置和局限性 | 报告二稿 |
| 5.25 | 组内交叉检查公式、代码、图表一致性 | 修改清单 |
| 5.26 | 制作课堂展示 | 展示初稿 |
| 5.27 | 第一次试讲，记录问题 | 试讲反馈 |
| 5.28 | 修改展示节奏和重点图表 | 展示二稿 |
| 5.29 | 第二次试讲，压缩讲述时间 | 最终展示 |
| 5.30 | 最终检查报告、代码、图表、引用 | 提交包 |
| 5.31 | 预留缓冲时间 | 最终版本 |

## 16. 报告结构

建议最终报告结构如下：

1. 摘要
2. 引言
3. 问题建模
4. 奇异值分解与病态性
5. 伪逆求解的误差放大机制
6. Tikhonov 正则化与谱滤波解释
7. 一维信号实验设计
8. 一维实验结果与分析
9. 二维图像实验设计
10. 二维实验结果与分析
11. 正则化参数选择
12. 讨论：理论与实验的一致性
13. 局限性与改进方向
14. 结论
15. 参考文献
16. 附录：代码结构、参数表和复现说明

## 17. 课堂展示结构

建议展示控制在 10 至 12 分钟。

| 时间 | 内容 | 目标 |
|---|---|---|
| 1 分钟 | 问题引入：微小噪声为什么会造成巨大错误 | 建立直觉 |
| 2 分钟 | SVD、条件数和伪逆误差放大 | 展示高代核心 |
| 2 分钟 | Tikhonov 正则化的谱滤波解释 | 说明稳定机制 |
| 3 分钟 | 一维实验：误差曲线、滤波因子、L-curve | 理论验证 |
| 2 分钟 | 图像实验：退化、伪逆失败、正则化恢复 | 直观展示 |
| 1 分钟 | 结论与局限性 | 收束观点 |

展示中最重要的视觉逻辑：

1. 先展示奇异值快速衰减；
2. 再展示 $\frac1{\sigma_i}$ 爆炸；
3. 再展示 $\frac{\sigma_i}{\sigma_i^2+\lambda}$ 被压住；
4. 最后展示伪逆恢复失败和正则化恢复改善。

## 18. 答辩可能问题与准备

### 问题 1：为什么不用普通逆矩阵？

因为病态问题中即使 $A$ 可逆，$A^{-1}$ 也可能对扰动极端敏感。伪逆或逆矩阵中的 $\frac1{\sigma_i}$ 会放大小奇异值方向上的噪声。

### 问题 2：Tikhonov 为什么能改善结果？

它将伪逆中的放大因子 $\frac1{\sigma_i}$ 替换为 $\frac{\sigma_i}{\sigma_i^2+\lambda}$，从而抑制小奇异值方向上的噪声放大。

### 问题 3：$\lambda$ 越大越好吗？

不是。过大的 $\lambda$ 会过度压制真实信号，导致偏差增大。最优 $\lambda$ 是稳定性和准确性之间的折中。

### 问题 4：为什么用随机正交矩阵构造 $A$？

随机正交矩阵保证 $U,V$ 不改变向量范数，病态性完全由奇异值谱控制。这使实验能准确研究奇异值衰减对误差放大的影响。

### 问题 5：图像实验和一维实验有什么数学联系？

图像矩阵可以向量化为高维向量。二维退化模型 $Y=BXB^T$ 向量化后变成 $\mathrm{vec}(Y)=(B\otimes B)\mathrm{vec}(X)$，仍是线性系统。

### 问题 6：为什么二维问题更病态？

因为：

$$
\kappa(B\otimes B)=\kappa(B)^2
$$

如果 $B$ 已经病态，Kronecker 积构造的二维算子会更加病态。

## 19. 风险与控制

### 风险 1：矩阵太大导致计算慢

控制方法：

- 一维实验先取 $n=256$；
- 图像实验先取 $m=64$；
- 只在最终展示中尝试 $m=128$；
- 优先使用 SVD 结构计算，不盲目构造超大矩阵。

### 风险 2：伪逆结果数值溢出或图像全黑全白

控制方法：

- 图像展示时对结果单独归一化用于显示；
- 误差计算必须使用原始数值，不用显示归一化后的数据；
- 报告中说明显示归一化和数值评价的区别。

### 风险 3：正则化参数选得不好

控制方法：

- 使用对数网格扫描；
- 同时报告 oracle、L-curve 和残差原则；
- 展示 $\lambda$ 过小、合适、过大的三种结果。

### 风险 4：报告变成图像处理项目

控制方法：

- 理论部分以 SVD 和条件数为主；
- 图像只作为线性代数现象的可视化案例；
- 每张图都回扣奇异值谱、伪逆放大或谱滤波。

## 20. 最终交付物

最终应提交：

1. 完整研究报告；
2. 课堂展示 PPT 或网页 PPT；
3. 可复现实验代码；
4. 所有生成图表；
5. 实验参数记录；
6. 数据说明；
7. 小组分工说明；
8. 参考文献；
9. 复现说明。

## 21. 最终结论预期

本项目预期得出如下结论：

1. 病态线性系统的核心困难来自奇异值谱的快速衰减。
2. 伪逆求解会沿小奇异值方向显著放大观测噪声。
3. 条件数越大，系统对扰动越敏感，恢复误差通常越大。
4. Tikhonov 正则化通过谱滤波抑制小奇异值方向上的噪声放大。
5. 正则化参数 $\lambda$ 控制稳定性和准确性的折中。
6. 一维信号与二维图像实验都能验证同一套高等代数理论。

项目的核心价值在于：用高等代数中的奇异值分解，把“数值计算为什么会失败”和“正则化为什么有效”解释清楚，并通过可复现实验给出量化验证。
