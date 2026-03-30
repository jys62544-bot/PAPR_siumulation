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
algo_names = {'Original', 'Clipping (CR=1.2)', 'mu-law (mu=10)', ...
              'A-law (A=87.6)', 'SLM (U=8)', 'PTS (V=4)', 'Tone Res. (10%)'};
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

% 3) μ 律压扩
fprintf('计算：μ 律压扩...\n');
for i = 1:N_sym
    [x_mu, ~] = companding_mu(x_orig_all(:,i), 10);
    papr_results(i, 3) = compute_papr(x_mu);
end

% 4) A 律压扩
fprintf('计算：A 律压扩...\n');
for i = 1:N_sym
    [x_a, ~] = companding_a(x_orig_all(:,i), 87.6);
    papr_results(i, 4) = compute_papr(x_a);
end

% 5) SLM
fprintf('计算：SLM (U=8)...\n');
P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N, 1);
for i = 1:N_sym
    [~, ~, papr_results(i, 5)] = slm(X_all(:,i), params, 8, P_slm);
end

% 6) PTS
fprintf('计算：PTS (V=4)...\n');
W_pts = [1, -1, 1j, -1j];
for i = 1:N_sym
    [~, ~, papr_results(i, 6)] = pts(X_all(:,i), params, 4, 'interleaved', W_pts);
end

% 7) 预留音调
fprintf('计算：预留音调...\n');
for i = 1:N_sym
    [x_tr, ~, ~] = tone_reservation(X_all(:,i), params, 0.10, 10);
    papr_results(i, 7) = compute_papr(x_tr);
end

% 保存数据
save('results/data/ccdf_results.mat', 'papr_results', 'algo_names', 'params');

% 绘制 CCDF 曲线
fig = figure('Position', [100, 100, 800, 600]);
hold on;
for a = 1:n_algo
    [pa, cc] = compute_ccdf(papr_results(:, a));
    mk_idx = 1:max(1,floor(N_sym/15)):N_sym;
    plot(pa, cc, 'Color', colors(a,:), 'LineStyle', lstyles{a}, ...
         'LineWidth', 1.5, 'DisplayName', algo_names{a});
    plot(pa(mk_idx), cc(mk_idx), markers{a}, 'Color', colors(a,:), ...
         'MarkerSize', 6, 'MarkerFaceColor', colors(a,:), ...
         'HandleVisibility', 'off');
end
hold off;

set(gca, 'YScale', 'log');
xlabel('PAPR_0 (dB)');
ylabel('Pr(PAPR > PAPR_0)');
title('各 PAPR 降低算法的 CCDF 对比');
legend('Location', 'southwest');
xlim([4 14]);
ylim([1e-3 1]);
grid on;

save_figure(fig, 'fig01_ccdf_comparison');
fprintf('实验 1 完成。\n');
