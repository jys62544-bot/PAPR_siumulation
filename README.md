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

详见 [LIMITATIONS.md](LIMITATIONS.md)，包括信道模型简化、SNR 定义（Es/N0 vs Eb/N0）、各算法实现局限、以及仿真框架假设等。

## 许可证

本项目用于孙剑宇的西北工业大学MIMO-OFDM无线通信技术课程设计
