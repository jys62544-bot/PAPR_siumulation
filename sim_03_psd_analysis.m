% SIM_03_PSD_ANALYSIS
% 对比 PAPR 降低后各信号的功率谱密度
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);

N = params.N_fft;
N_os = params.N_os;
N_sym_psd = 200;  % 用于 PSD 平均估计的符号数

fprintf('实验 3：PSD 分析\n');

% 拼接多个符号用于频谱估计
x_concat = struct();
algo_names = {'Original', 'Clipping (CR=0.8)', 'Clipping+Filter (CR=0.8)', ...
              'mu-law (mu=10)', 'SLM (U=8)', 'PTS (V=4)'};
n_algo = length(algo_names);

% 初始化信号缓冲区
for a = 1:n_algo
    x_concat(a).data = [];
end

P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N,1);
W_pts = [1, -1, 1j, -1j];

for i = 1:N_sym_psd
    bits = randi([0 1], N * params.bps, 1);
    X = ofdm_mod(bits, params.mod_type);
    [x_orig, ~] = ofdm_transmitter(X, params);
    % 1) 原始信号
    x_concat(1).data = [x_concat(1).data; x_orig];
    % 2) 仅限幅
    A = 0.8 * rms(x_orig);
    mag = abs(x_orig);
    x_clip_only = min(mag, A) .* exp(1j * angle(x_orig));
    x_concat(2).data = [x_concat(2).data; x_clip_only];
    % 3) 限幅 + 滤波
    [x_cf, ~] = clipping_filtering(x_orig, params, 0.8);
    x_concat(3).data = [x_concat(3).data; x_cf];
    % 4) μ 律压扩
    [x_mu, ~] = companding_mu(x_orig, 10);
    x_concat(4).data = [x_concat(4).data; x_mu];
    % 5) SLM
    [x_slm, ~, ~] = slm(X, params, 8, P_slm);
    x_concat(5).data = [x_concat(5).data; x_slm];
    % 6) PTS
    [x_pts, ~, ~] = pts(X, params, 4, 'interleaved', W_pts);
    x_concat(6).data = [x_concat(6).data; x_pts];
end
% 使用 Welch 方法计算并绘制 PSD
fig = figure('Position', [100, 100, 900, 600]);
hold on;
Fs = 1;% 归一化频率
nfft_psd = 1024;

for a = 1:n_algo
    [pxx, f] = pwelch(x_concat(a).data, hamming(nfft_psd), nfft_psd/2, nfft_psd, Fs, 'centered');
    pxx_dB = 10*log10(pxx / max(pxx));  % 归一化到 0 dB 峰值
    plot(f, pxx_dB, 'Color', colors(a,:), 'LineStyle', lstyles{a}, ...
         'LineWidth', 1.5, 'DisplayName', algo_names{a});
end
hold off;
xlabel('归一化频率');
ylabel('功率谱密度（dB）');
title('PSD 对比：频谱再生分析');
legend('Location', 'south');
ylim([-60 5]);
grid on;
save_figure(fig, 'fig04_psd_comparison');
% 保存数据
save('results/data/psd_results.mat', 'algo_names', 'params');
fprintf('实验 3 完成。\n');
