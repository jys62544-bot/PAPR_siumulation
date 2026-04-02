% SIM_01_CCDF_COMPARISON
% 对比所有 PAPR 降低算法的 CCDF 曲线。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);

N_sym = params.N_sym;
N = params.N_fft;

fprintf('实验 1：CCDF 对比\n');
fprintf('N_fft=%d, L=%d, %s, %d 符号\n', N, params.L, params.mod_type, N_sym);

% 预生成所有 OFDM 符号（频域）
X_all = zeros(N, N_sym);
for i = 1:N_sym
    bits = randi([0 1], N * params.bps, 1);
    X_all(:, i) = ofdm_mod(bits, params.mod_type);
end

% 预生成时域原始信号
x_orig_all = zeros(params.N_os, N_sym);
for i = 1:N_sym
    [x_orig_all(:,i), ~] = ofdm_transmitter(X_all(:,i), params);
end

% 待对比的算法列表
algo_names = {'Original', 'Clipping (CR=1.2)', 'SLM (U=8)', ...
              'PTS (V=4)', 'Tone Res. (10%)', 'mu-law (mu=10)', ...
              'Golay Coding'};
n_algo = length(algo_names);
papr_results = zeros(N_sym, n_algo);

% 1) 原始信号
fprintf('计算：原始信号...\n');
for i = 1:N_sym
    papr_results(i, 1) = compute_papr(x_orig_all(:, i));
end

% 2) 限幅滤波
fprintf('计算：限幅滤波...\n');
for i = 1:N_sym
    [x_cl, ~] = clipping_filtering(x_orig_all(:,i), params, 1.2);
    papr_results(i, 2) = compute_papr(x_cl);
end

% 3) SLM
fprintf('计算：SLM (U=8)...\n');
P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N, 1);
for i = 1:N_sym
    [~, ~, papr_results(i, 3)] = slm(X_all(:,i), params, 8, P_slm);
end

% 4) PTS
fprintf('计算：PTS (V=4)...\n');
W_pts = [1, -1, 1j, -1j];
for i = 1:N_sym
    [~, ~, papr_results(i, 4)] = pts(X_all(:,i), params, 4, 'interleaved', W_pts);
end

% 5) 预留音调
fprintf('计算：预留音调...\n');
for i = 1:N_sym
    [x_tr, ~, ~] = tone_reservation(X_all(:,i), params, 0.10, 10);
    papr_results(i, 5) = compute_papr(x_tr);
end

% 6) μ 律压扩
fprintf('计算：μ 律压扩...\n');
for i = 1:N_sym
    [x_mu, ~] = companding_mu(x_orig_all(:,i), 10);
    papr_results(i, 6) = compute_papr(x_mu);
end

% 7) Golay 互补序列编码
fprintf('计算：Golay 编码...\n');
for i = 1:N_sym
    [x_gl, ~, gl_info] = golay_coding(N, params, 16);
    papr_results(i, 7) = gl_info.papr;
end

% 保存数据
save('results/data/ccdf_results.mat', 'papr_results', 'algo_names', 'params');

% 绘制 CCDF 曲线
fig = figure('Position', [100, 100, 900, 600]);
ax = axes(fig);
hold(ax, 'on');
for a = 1:n_algo
    [pa, cc] = compute_ccdf(papr_results(:, a));
    mk_idx = 1:max(1,floor(N_sym/15)):N_sym;
    plot(ax, pa, cc, 'Color', colors(a,:), 'LineStyle', lstyles{a}, ...
         'LineWidth', 1.5, 'DisplayName', algo_names{a});
    plot(ax, pa(mk_idx), cc(mk_idx), markers{a}, 'Color', colors(a,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(a,:), ...
         'HandleVisibility', 'off');
end
hold(ax, 'off');

set(ax, 'YScale', 'log');
xlabel(ax, 'PAPR_0 (dB)');
ylabel(ax, 'Pr(PAPR > PAPR_0)');
title(ax, '各 PAPR 降低算法的 CCDF 对比');
xlim(ax, [2 14]);
ylim(ax, [1e-3 1]);
grid(ax, 'on');

% 图例放在绘图区域外右侧，避免遮挡
lg = legend(ax, 'Location', 'eastoutside');
set(lg, 'FontSize', 9);

save_figure(fig, 'fig01_ccdf_comparison');
fprintf('实验 1 完成。\n');
