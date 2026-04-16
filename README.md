# OFDM PAPR 降低技术仿真

基于 MATLAB 的 OFDM 峰均功率比（PAPR）降低技术仿真平台，涵盖经典信号处理方法与深度学习方法的全面对比分析。系统包括接收机仿真、发射机仿真、循环前缀仿真、多径信道仿真、星座映射与解映射仿真、高斯信道仿真、功率谱密度分析、计算复杂度分析、参数敏感性扫描以及BER vs SNR 对比和CCDF 曲线对比。此外，系统还支持基于信道状态信息的逐子载波自适应调制（BPSK/QPSK/16QAM/64QAM/256QAM），可根据各子载波的有效 SNR 动态分配调制方式以最大化频谱效率。

## 项目概述

峰均功率比（PAPR）是 OFDM 系统中的核心问题之一。过高的 PAPR 会导致功率放大器工作在非线性区，引起信号失真和效率下降。本项目实现并对比了 7 种主流 PAPR 降低算法，包括基于 Davis & Jedwab (1999) 论文的真正 Golay 互补序列编码（PMEPR 严格 ≤ 3 dB），以及一种基于深度神经网络的"限幅+逐子载波 DNN 补偿"方法。

## 目录结构

```
OFDM/
├── main_run_all.m              # 一键运行全部 6 个实验
├── sim_01_ccdf_comparison.m    # 实验 1：CCDF 曲线对比
├── sim_02_ber_comparison.m     # 实验 2：BER vs SNR 对比
├── sim_03_psd_analysis.m       # 实验 3：功率谱密度分析
├── sim_04_complexity.m         # 实验 4：计算复杂度分析
├── sim_05_param_sweep.m        # 实验 5：参数敏感性扫描
├── sim_06_dnn_papr.m           # 实验 6：DNN 限幅+逐子载波补偿方法
│
├── core/                       # OFDM 系统核心模块
│   ├── get_default_params.m    # 默认系统参数
│   ├── ofdm_mod_adaptive.m     # 逐子载波自适应调制映射
│   ├── ofdm_demod_adaptive.m   # 逐子载波自适应解调映射
│   ├── adaptive_modulation.m   # 信道自适应调制阶数分配
│   ├── ofdm_transmitter.m      # OFDM 发射机（过采样 IFFT + CP）
│   ├── ofdm_receiver.m         # OFDM 接收机（去 CP + FFT + 均衡）
│   ├── cp_add.m                # 循环前缀添加
│   ├── channel_awgn.m          # AWGN 信道
│   ├── channel_awgn_fixed.m    # 固定参考功率 AWGN 信道（公平比较用）
│   └── channel_multipath.m     # ITU 多径衰落信道（EPA/ETU）
│
├── papr_reduction/             # PAPR 降低算法
│   ├── clipping_filtering.m    # 限幅滤波法
│   ├── slm.m                   # 选择映射法（SLM）
│   ├── pts.m                   # 部分传输序列法（PTS）
│   ├── tone_reservation.m      # 预留音调法（TR）
│   ├── companding_mu.m         # μ 律压扩法
│   ├── golay_coding.m          # Golay Davis-Jedwab 编码（QPSK，PMEPR ≤ 3dB）
│   ├── golay_decode.m          # Golay 解码器（FHT 算法）
│   ├── golay_build_perm_table.m # Golay coset 排列表生成
│   └── dnn_papr_reduction.m    # DNN 限幅+逐子载波补偿法（含条件特征）
│
├── analysis/                   # 分析与可视化工具
│   ├── compute_papr.m          # 计算 PAPR（dB）
│   ├── compute_ccdf.m          # 计算 CCDF
│   ├── plot_config.m           # 统一绘图样式配置
│   └── save_figure.m           # 保存图形为 .fig 和 .png
│
├── results/                    # 仿真结果保存
│   ├── figures/                # 输出图像（.fig + .png）
│   └── data/                   # 仿真数据（.mat）
│
└── LIMITATIONS.md              # 仿真系统局限性说明
```

## 系统参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| N_fft | 256 | FFT 点数&子载波数 |
| N_cp | 64 | 循环前缀长度16 |
| 调制方式 | 逐子载波自适应 | 根据各子载波信道 SNR 动态选择 BPSK/QPSK/16QAM/64QAM/256QAM |
| 自适应调制 | 开启（固定装载） | 在工作 SNR=15dB 下一次性分配，所有实验共用（非动态 AMC） |
| 信道模型 | ETU | ITU ETU 多径信道，可选 EPA |
| SNR 门限 | [6,12,18,24] dB | 自适应调制切换门限 |
| 过采样倍数 L | 4 | 精确捕捉连续时间峰值 |
| 仿真符号数 | 10000 | 保证 CCDF 可靠至 1e-3 |
| SNR 范围 | 0~30 dB | 步进 2 dB |

## 已实现算法

| 算法 | 类型 | 关键参数 | 复杂度 |
|------|------|---------|--------|
| 限幅滤波（CF） | 有失真 | 限幅比 CR | O(NL) |
| 选择映射（SLM） | 无失真 | 候选数 U | O(U·NL log NL) |
| 部分传输序列（PTS） | 无失真 | 子块数 V，相位因子集 W | O(\|W\|^(V-1)·NL) |
| 预留音调（TR） | 无失真 | 预留比例，迭代次数 | O(I·NL log NL) |
| μ 律压扩 | 有失真 | μ 参数 | O(NL) |
| Golay Davis-Jedwab 编码 | 编码型 | QPSK，PMEPR ≤ 3dB | O(N·m) |
| DNN 限幅+逐子载波补偿 | 学习型 | 限幅比 CR | 离线训练 + O(N) 推理 |

## 实验说明

### 实验 1：CCDF 对比（`sim_01_ccdf_comparison.m`）
对比 7 种方案的互补累积分布函数，以 PAPR 超越概率 10^-3 作为性能指标。输出：`fig01_ccdf_comparison.png`。

### 实验 2：BER 对比（`sim_02_ber_comparison.m`）
在多径信道下测量 7 种算法的误比特率，公平比较原则：所有算法使用同一组输入数据；噪声功率统一基于原始信号功率计算（`channel_awgn_fixed`）；TR 按有效信息比特统计 BER；Golay-DJ 独立生成 QPSK Golay 序列，以自身信号功率为参考，按 32 info bits/symbol 统计 BER。输出：`fig02_ber_comparison.png`。

### 实验 3：PSD 分析（`sim_03_psd_analysis.m`）
使用 Welch 法估计功率谱密度，重点对比限幅与限幅+滤波在带外频谱再生方面的差异。输出：`fig04_psd_comparison.png`。

### 实验 4：复杂度分析（`sim_04_complexity.m`）
实测各算法的 MATLAB 运行时间，并与理论大 O 复杂度对照，绘制 PAPR 降低量与计算代价的权衡散点图。输出：`fig05_complexity_bar.png`、`fig06_complexity_tradeoff.png`。

### 实验 5：参数敏感性（`sim_05_param_sweep.m`）
分别扫描各算法的关键参数（CF 的 CR、μ 律的 μ、SLM 的 U、PTS 的 V 和 W、TR 的预留比例），输出参数-性能曲线。输出：`fig07~fig11`。

### 实验 6：DNN 限幅+逐子载波补偿（`sim_06_dnn_papr.m`）

采用"**限幅 + DNN 逐子载波补偿**"的两阶段方案，结合传统限幅的 PAPR 降低能力与深度学习的失真补偿能力。

#### DNN 限幅+逐子载波补偿方法原理

**核心思想**：发射端用硬限幅将 PAPR 强制降低，接收端用 DNN 学习补偿限幅引入的非线性失真，从而在获得 PAPR 增益的同时最大程度恢复 BER 性能。

**1. 发射端——硬限幅**

对过采样时域信号 $x[n]$ 做幅度限幅：

$$x_{\text{clip}}[n] = \min(|x[n]|, A) \cdot e^{j\angle x[n]}$$

其中限幅门限 $A = \text{CR} \cdot \text{rms}(x)$，CR 为限幅比（默认 1.2）。限幅后信号的 PAPR 被强制限制在较低水平（实测均值约 2.76 dB，相比原始 8.38 dB 降低 5.6 dB）。

**2. 限幅失真分析**

限幅是时域非线性操作，会在频域引入子载波间的"限幅噪声"：

$$X_{\text{clip}}[k] = X[k] + D[k]$$

其中 $D[k]$ 是限幅失真分量，分布在所有子载波上。该失真是确定性的（取决于输入信号），非高斯噪声，因此可用 DNN 学习补偿。

**3. 接收端——DNN 逐子载波补偿**

关键设计：**将 256 维的全局补偿问题分解为 256 个独立的 4→2 映射**。

- **输入**：单个子载波的接收 IQ 值 $(I_{\text{rx}}, Q_{\text{rx}})$、归一化调制阶数、归一化信道增益
- **输出**：补偿后的 IQ 值 $(\hat{I}, \hat{Q})$
- **网络结构**：3 层全连接网络（4→128→64→2），ReLU 激活
- **所有子载波共享同一网络**

这种逐子载波设计的优势：
- **训练数据量大**：$N_{\text{train}} \times N_{\text{fft}}$ = 20000 × 256 = 512 万逐子载波样本
- **网络极小**（仅 ~9000 参数），训练快速收敛
- **完美泛化**：限幅失真在每个子载波上的统计特性相似，网络在测试集上直接生效

**4. 训练过程**

- 目标函数：MSE 损失 $\mathcal{L} = \frac{1}{B}\sum_{i=1}^{B}\|(\hat{I}_i, \hat{Q}_i) - (I_i, Q_i)\|^2$
- 优化器：Adam（lr=1e-3）
- 训练数据包含信道效应（多径 + 多 SNR 噪声），提升推理泛化能力
- 30 个 epoch 即可收敛

**5. 实验结果**

| 指标 | Original | 纯限幅 (CR=1.2) | 限幅+DNN |
|------|----------|-----------------|----------|
| PAPR 均值 | 8.38 dB | **2.76 dB** | **2.76 dB** |
| BER@10dB | 1.78e-3 | 2.77e-2 | **1.58e-2** |
| BER@20dB | 0 | 8.04e-3 | **1.93e-3** |
| BER@30dB | 0 | 6.36e-3 | **1.31e-3** |

DNN 补偿将限幅引入的 BER 降低了 2~5 倍，同时保持了与纯限幅完全相同的 PAPR 降低效果。

输出：`fig12_dnn_training`（训练损失）、`fig13_dnn_ccdf`（CCDF）、`fig14_dnn_ber`（BER）、`fig15_dnn_constellation`（星座图）。

> **注意**：实验 6 需要安装 MATLAB Deep Learning Toolbox。DNN 训练和测试基于单个固定的信道实现，结果反映的是对特定信道实例的补偿能力，非通用补偿器。

## 运行方法

### 运行全部实验
在 MATLAB 命令窗口中，将工作目录切换至项目根目录后执行：
```matlab
main_run_all
```

### 单独运行某个实验
```matlab
addpath('core', 'analysis', 'papr_reduction');
sim_01_ccdf_comparison   % 或其他实验脚本
```

### 修改仿真参数
编辑 `core/get_default_params.m` 或在脚本中直接覆盖：
```matlab
params = get_default_params();
params.mod_type = 'QPSK';   % 修改调制方式
params.N_sym = 5000;        % 减少仿真符号数以加快速度
```

## 输出文件

运行后结果自动保存至：
- `results/figures/`：所有图形（`.fig` 和 300 dpi `.png`）
- `results/data/`：仿真原始数据（`.mat`）

| 文件名 | 内容 |
|--------|------|
| `fig01_ccdf_comparison` | 各算法 CCDF 对比 |
| `fig02_ber_comparison` | 多径信道 BER 对比 |
| `fig04_psd_comparison` | 功率谱密度对比 |
| `fig05_complexity_bar` | 运行时间柱状图 |
| `fig06_complexity_tradeoff` | PAPR 降低 vs 复杂度散点图 |
| `fig07~fig11` | 各算法参数敏感性曲线 |
| `fig12_dnn_training` | DNN 逐子载波补偿器训练损失曲线 |
| `fig13_dnn_ccdf` | DNN 方案 vs 经典方法 CCDF |
| `fig14_dnn_ber` | Original vs 纯限幅 vs 限幅+DNN BER |
| `fig15_dnn_constellation` | 原始/限幅/DNN补偿 星座图对比 |

## 软件要求

- MATLAB R2021a 或更高版本
- Signal Processing Toolbox
- Deep Learning Toolbox（仅实验 6 需要）

## 局限性说明

### 1. 信道模型局限

#### 1.1 单信道实现

所有仿真使用**一次随机生成的固定信道实现**（`channel_multipath` 只调用一次）。这意味着：
- 所有 SNR 扫描点共用同一信道频率响应
- 结果反映的是特定信道实例下的性能，非多信道统计平均
- 实际系统中信道时变，性能会有波动

#### 1.2 完美信道估计

接收端使用发射端已知的理想信道频率响应 `H_ch` 做均衡，不考虑信道估计误差。实际系统中导频开销和估计误差会进一步降低性能。

#### 1.3 ZF 均衡

使用零迫（ZF）均衡 `X_hat = Y ./ H`，在信道深衰落子载波上会严重放大噪声。实际系统通常使用 MMSE 均衡以获得更好的噪声-干扰权衡。

### 2. 公平比较局限

#### 2.1 SNR 定义

BER 仿真中 SNR 定义为 **Es/N0**（每符号能量与噪声功率谱密度之比），而非 Eb/N0（每比特能量）。对于编码率不同的算法：
- **Golay 编码**（Rate-1/2）：相同 Es/N0 下 Eb/N0 比其他算法高约 **3 dB**，BER 曲线相对占优
- **Tone Reservation**（~10% 预留）：Eb/N0 高约 **0.46 dB**

图中 x 轴标注的 SNR 为 Es/N0，阅读 BER 曲线时需注意此差异。

#### 2.2 Side Information

SLM 和 PTS 假设接收端**完美已知 side information**（所选相位序列索引或最优相位因子）。实际系统中 side information 需要额外比特传输，存在传输错误风险，且会降低有效数据率。Golay 编码的序列索引（候选编号）同样假设完美已知。

### 3. 算法实现局限

#### 3.1 μ-law 压扩

- **解压扩参考幅度不匹配**：发射端压扩使用发射信号的 `v_max`，接收端解压扩使用接收信号的 `v_max`。经过信道后两者不同，导致压扩/解压扩不完全互逆，引入额外失真。
- **多径信道下噪声放大**：接收链路为 FFT→ZF均衡→IFFT→逆压扩→FFT→解调。ZF 均衡在深衰落子载波上放大的噪声，会被逆压扩的指数函数进一步放大，导致 mu-law 在多径信道下 BER 性能劣化明显。这是 mu-law 压扩方法在频率选择性信道中的固有缺陷。

#### 3.2 Tone Reservation

- 目标限幅门限硬编码为 `3.0 × RMS`，未作为可调参数暴露
- 步长因子 `mu_step = 0.5` 固定，可能非最优
- 迭代次数 30 次对部分符号可能不足以收敛

#### 3.3 Golay Davis-Jedwab 编码

- 基于 Davis & Jedwab (1999) 论文实现，每个 OFDM 码字是 Z_4 上长度 N 的 Golay 序列，**PMEPR 严格 ≤ 2（即 ≤ 3.01 dB）**
- 调制方式固定为 **QPSK**（h=4），不兼容自适应调制或更高阶调制
- **编码率极低**：N=256 (m=8) 下，每 OFDM 符号仅传输 32 信息比特（14 bits coset 选择 + 18 bits 线性系数），编码率 = 32/512 ≈ 1/16。这是论文方案在大 N 下的固有特性（论文原文主要针对 N=16~32）
- 接收端假设 **coset 索引（排列 π）完美已知**（side information），实际系统中需额外比特传输
- 解码使用 **硬判决 + 两次 Fast Hadamard Transform**（Algorithm 14, h=4），未实现软判决
- 排列表的去重方式（去除逆序重复）可能未覆盖所有 m!/2 个等价类，但不影响 PMEPR 保证

#### 3.4 DNN 限幅+逐子载波补偿

- **单信道泛化**：DNN 在单个固定信道实现上训练和测试，不具备跨信道泛化能力
- **无验证集/早停**：训练 30 epochs 无过拟合监控，对于当前小网络（~9000 参数）和大数据量（512 万样本）过拟合风险较低，但不是最佳实践
- **推理延迟未考虑**：复杂度分析中 DNN 的推理时间与 MATLAB 的 `predict` 函数实现相关，不代表实际硬件部署性能

### 4. 仿真框架局限

#### 4.1 BER 统计精度

- 最大仿真比特数 `max_bits = 1e6`，在高 SNR（>25 dB）区间 BER 低于 1e-5 时统计样本不足
- Golay 编码每 OFDM 符号的有效信息比特远少于其他算法，高 SNR 下误比特积累更慢
- 提前停止条件要求所有 7 种算法均达到 100 个错误比特，低 BER 算法可能拖慢整体进度

#### 4.2 自适应调制简化

- 使用**固定装载**方案：在工作 SNR = 15 dB 下一次性分配 `mod_map`，所有 SNR 扫描点共用
- 非真正的动态 AMC（Adaptive Modulation and Coding），不会随 SNR 变化重新分配
- 实际系统中调制方式会随信道和 SNR 变化动态调整

#### 4.3 无同步误差

系统假设完美的时间和频率同步，不考虑：
- 符号定时偏移（STO）
- 载波频率偏移（CFO）
- 采样时钟偏移（SCO）

这些在实际 OFDM 系统中会显著影响性能。

#### 4.4 无功率放大器模型

PAPR 降低的实际目的是使功率放大器（PA）工作在线性区。当前仿真没有引入 PA 非线性模型（如 Rapp 模型或 Saleh 模型），因此无法直接展示 PAPR 降低对 PA 效率和信号质量的实际改善。

### 5. 代码实现注意事项

| 项目 | 说明 |
|------|------|
| `clipping_filtering` 返回的 `X_clipped` | 未做 `sqrt(N_os)` 归一化，比 `ofdm_receiver` 输出大 `sqrt(N_os)` 倍。当前仿真脚本未使用此返回值做解调，不影响结果 |
| `compute_papr` | 输入全零信号时会产生 `Inf`/`NaN`（除零） |
| PTS V=8 参数扫描 | W={±1, ±j} 时搜索 4^7=16384 种组合，运行时间极长 |
| `save_figure` | 临时将文本 Interpreter 设为 'none'，保存后不恢复原设置 |

## 许可证

本项目用于孙剑宇的西北工业大学MIMO-OFDM无线通信技术课程设计
