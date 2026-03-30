# OFDM PAPR 降低仿真套件

基于 MATLAB 的 OFDM 峰均功率比（PAPR）降低算法仿真平台，涵盖经典信号处理方法与深度学习方法的全面对比分析。

## 项目概述

峰均功率比（PAPR）是 OFDM 系统中的核心问题之一。过高的 PAPR 会导致功率放大器工作在非线性区，引起信号失真和效率下降。本项目实现并对比了 6 种主流 PAPR 降低算法，并新增了一种基于深度神经网络自编码器的方法。

## 目录结构

```
OFDM/
├── main_run_all.m              # 一键运行全部 6 个实验
├── sim_01_ccdf_comparison.m    # 实验 1：CCDF 曲线对比
├── sim_02_ber_comparison.m     # 实验 2：BER vs SNR 对比
├── sim_03_psd_analysis.m       # 实验 3：功率谱密度分析
├── sim_04_complexity.m         # 实验 4：计算复杂度分析
├── sim_05_param_sweep.m        # 实验 5：参数敏感性扫描
├── sim_06_dnn_papr.m           # 实验 6：DNN 自编码器方法
│
├── core/                       # OFDM 系统核心模块
│   ├── get_default_params.m    # 默认系统参数
│   ├── ofdm_mod.m              # 星座映射（QPSK/16QAM/64QAM）
│   ├── ofdm_demod.m            # 星座解映射（硬判决）
│   ├── ofdm_transmitter.m      # OFDM 发射机（过采样 IFFT）
│   ├── ofdm_receiver.m         # OFDM 接收机（FFT + 均衡）
│   ├── channel_awgn.m          # AWGN 信道
│   └── channel_multipath.m     # ITU 多径衰落信道（EPA/ETU）
│
├── papr_reduction/             # PAPR 降低算法
│   ├── clipping_filtering.m    # 限幅滤波法
│   ├── companding_mu.m         # μ 律压扩法
│   ├── companding_a.m          # A 律压扩法
│   ├── slm.m                   # 选择映射法（SLM）
│   ├── pts.m                   # 部分传输序列法（PTS）
│   ├── tone_reservation.m      # 预留音调法（TR）
│   └── dnn_papr_reduction.m    # DNN 自编码器法
│
├── analysis/                   # 分析与可视化工具
│   ├── compute_papr.m          # 计算 PAPR（dB）
│   ├── compute_ccdf.m          # 计算 CCDF
│   ├── plot_config.m           # 统一绘图样式配置
│   └── save_figure.m           # 保存图形为 .fig 和 .png
│
└── results/                    # 仿真结果（自动生成）
    ├── figures/                # 输出图像（.fig + .png）
    └── data/                   # 仿真数据（.mat）
```

## 系统参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| N_fft | 256 | FFT 点数（子载波数） |
| N_cp | 64 | 循环前缀长度（N_fft/4） |
| 调制方式 | 16QAM | 支持 QPSK / 16QAM / 64QAM |
| 过采样倍数 L | 4 | 精确捕捉连续时间峰值 |
| 仿真符号数 | 10000 | 保证 CCDF 可靠至 1e-3 |
| SNR 范围 | 0~30 dB | 步进 2 dB |
| 信道模型 | AWGN | 可选 ITU EPA/ETU 多径 |

## 已实现算法

| 算法 | 类型 | 关键参数 | 复杂度 |
|------|------|---------|--------|
| 限幅滤波（CF） | 有失真 | 限幅比 CR | O(NL) |
| μ 律压扩 | 有失真 | μ 参数 | O(NL) |
| A 律压扩 | 有失真 | A 参数 | O(NL) |
| 选择映射（SLM） | 无失真 | 候选数 U | O(U·NL log NL) |
| 部分传输序列（PTS） | 无失真 | 子块数 V，相位因子集 W | O(\|W\|^(V-1)·NL) |
| 预留音调（TR） | 无失真 | 预留比例，迭代次数 | O(I·NL log NL) |
| DNN 自编码器 | 学习型 | 扰动预算 α，损失权重 | 离线训练 |

## 实验说明

### 实验 1：CCDF 对比（`sim_01_ccdf_comparison.m`）
对比 7 种方案（含原始信号）的互补累积分布函数，以 PAPR 超越概率 10^-3 作为性能指标。输出：`fig01_ccdf_comparison.png`。

### 实验 2：BER 对比（`sim_02_ber_comparison.m`）
在 AWGN 信道下测量各算法的误比特率，展示 PAPR 降低方法对通信质量的影响。输出：`fig02_ber_awgn.png`。

### 实验 3：PSD 分析（`sim_03_psd_analysis.m`）
使用 Welch 法估计功率谱密度，重点对比限幅与限幅+滤波在带外频谱再生方面的差异。输出：`fig04_psd_comparison.png`。

### 实验 4：复杂度分析（`sim_04_complexity.m`）
实测各算法的 MATLAB 运行时间，并与理论大 O 复杂度对照，绘制 PAPR 降低量与计算代价的权衡散点图。输出：`fig05_complexity_bar.png`、`fig06_complexity_tradeoff.png`。

### 实验 5：参数敏感性（`sim_05_param_sweep.m`）
分别扫描各算法的关键参数（CF 的 CR、μ 律的 μ、SLM 的 U、PTS 的 V 和 W、TR 的预留比例），输出参数-性能曲线。输出：`fig07~fig11`。

### 实验 6：DNN 自编码器（`sim_06_dnn_papr.m`）
训练编码器-解码器结构的 DNN，使用可微 PAPR 代理损失函数（p 范数近似）联合优化重建质量与 PAPR。与 SLM、PTS 在 CCDF 和 BER 上进行对比。输出：`fig12~fig15`。

> **注意**：实验 6 需要安装 MATLAB Deep Learning Toolbox。

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
| `fig02_ber_awgn` | AWGN 信道 BER 对比 |
| `fig04_psd_comparison` | 功率谱密度对比 |
| `fig05_complexity_bar` | 运行时间柱状图 |
| `fig06_complexity_tradeoff` | PAPR 降低 vs 复杂度散点图 |
| `fig07~fig11` | 各算法参数敏感性曲线 |
| `fig12_dnn_training` | DNN 训练损失曲线 |
| `fig13_dnn_ccdf` | DNN vs 经典方法 CCDF |
| `fig14_dnn_ber` | DNN vs 原始信号 BER |
| `fig15_dnn_constellation` | DNN 修改前后星座图 |

## 软件要求

- MATLAB R2021a 或更高版本
- Signal Processing Toolbox（`pwelch` 等函数）
- Deep Learning Toolbox（仅实验 6 需要）

## 许可证

本项目用于学术用途。
