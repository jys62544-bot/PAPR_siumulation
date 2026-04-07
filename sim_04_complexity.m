% SIM_04_COMPLEXITY
% 测量各 PAPR 降低算法的计算复杂度。
% 支持逐子载波自适应调制。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);

N = params.N_fft;
N_trials = 500;

fprintf('实验 4：复杂度分析\n');

%% 自适应调制设置
if params.adaptive
    dummy_x = zeros(params.N_os, 1);
    [~, H_ch, ~] = channel_multipath(dummy_x, params, params.channel_model);
    snr_work = 15;
    [mod_map, ~] = adaptive_modulation(H_ch, snr_work, params.snr_thresholds);
    fprintf('自适应调制，%s 信道，平均 bps=%.2f\n', params.channel_model, mean(mod_map));
end

% 预生成测试数据
if params.adaptive
    [X, ~, ~] = ofdm_mod_adaptive(mod_map);
else
    bits = randi([0 1], N * params.bps, 1);
    X = ofdm_mod(bits, params.mod_type);
end
[x_orig, ~] = ofdm_transmitter(X, params);

algo_names = {'Clipping', 'SLM(U=8)', ...
              'PTS(V=4)', 'Tone Res.', 'mu-law', 'Golay(N=16)'};
n_algo = length(algo_names);

theory_labels = {
    'O(NL)';
    'O(U \cdot NL\log NL)';
    'O(|W|^{V-1} \cdot NL)';
    'O(I \cdot NL\log NL)';
    'O(NL)';
    'O(C \cdot NL\log NL)';
};

runtimes = zeros(n_algo, 1);

P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N, 1);
W_pts = [1, -1, 1j, -1j];
% 预热
clipping_filtering(x_orig, params, 1.2);
companding_mu(x_orig, 10);
golay_coding(X, params, 16);

for trial = 1:N_trials
    if params.adaptive
        [X_t, ~, ~] = ofdm_mod_adaptive(mod_map);
    else
        bits_t = randi([0 1], N * params.bps, 1);
        X_t = ofdm_mod(bits_t, params.mod_type);
    end
    [x_t, ~] = ofdm_transmitter(X_t, params);

    tic; clipping_filtering(x_t, params, 1.2); runtimes(1) = runtimes(1) + toc;
    tic; slm(X_t, params, 8, P_slm); runtimes(2) = runtimes(2) + toc;
    tic; pts(X_t, params, 4, 'interleaved', W_pts); runtimes(3) = runtimes(3) + toc;
    tic; tone_reservation(X_t, params, 0.10, 30); runtimes(4) = runtimes(4) + toc;
    tic; companding_mu(x_t, 10); runtimes(5) = runtimes(5) + toc;
    tic; golay_coding(X_t, params, 16); runtimes(6) = runtimes(6) + toc;
end
runtimes = runtimes / N_trials * 1000;

fprintf('\n运行时间结果（毫秒/符号）：\n');
fprintf('%-20s %-15s %-10s\n', '算法', '复杂度', '运行时间(ms)');
fprintf('%s\n', repmat('-', 1, 45));
for a = 1:n_algo
    fprintf('%-20s %-15s %.3f\n', algo_names{a}, theory_labels{a}, runtimes(a));
end

fig1 = figure('Position', [100, 100, 900, 500]);
ax1 = axes(fig1);
b = bar(ax1, runtimes, 'FaceColor', 'flat');
for a = 1:n_algo
    b.CData(a,:) = colors(a+1, :);
end
set(ax1, 'XTickLabel', algo_names, 'TickLabelInterpreter', 'latex');
ylabel(ax1, '每 OFDM 符号运行时间（ms）');
title(ax1, '计算复杂度对比');
grid(ax1, 'on');
save_figure(fig1, 'fig05_complexity_bar');

if exist('results/data/ccdf_results.mat', 'file')
    ccdf_info = dir('results/data/ccdf_results.mat');
    if (now - ccdf_info.datenum) > 1  % 超过 1 天
        fprintf('警告：ccdf_results.mat 可能已过期（%.1f 天前），建议先运行 sim_01。\n', ...
            now - ccdf_info.datenum);
    end
    loaded = load('results/data/ccdf_results.mat');
    papr_at_1e3 = prctile(loaded.papr_results, 99.9);
    % papr_results 列: Original, Clipping, SLM, PTS, TR, mu-law, Golay（共 7 列）
    papr_reduction = papr_at_1e3(1) - papr_at_1e3(2:end);  % 6 个算法

    fig2 = figure('Position', [100, 100, 900, 600]);
    ax2 = axes(fig2);
    hold(ax2, 'on');
    for a = 1:n_algo
        scatter(ax2, runtimes(a), papr_reduction(a), 100, colors(a+1,:), 'filled', ...
            'Marker', markers{a+1}, 'DisplayName', algo_names{a});
        text(ax2, runtimes(a)*1.05, papr_reduction(a), algo_names{a}, ...
            'Interpreter', 'latex', 'FontSize', 10);
    end
    hold(ax2, 'off');
    xlabel(ax2, '每符号运行时间（ms）');
    ylabel(ax2, 'CCDF=10^{-3} 处的 PAPR 降低量（dB）');
    title(ax2, 'PAPR 降低量 vs 复杂度权衡');
    grid(ax2, 'on');
    lg = legend(ax2, 'Location', 'eastoutside');
    set(lg, 'FontSize', 9);
    save_figure(fig2, 'fig06_complexity_tradeoff');
end

save('results/data/complexity_results.mat', 'algo_names', 'runtimes', ...
     'theory_labels', 'params');
fprintf('实验 4 完成。\n');
