% SIM_05_PARAM_SWEEP
% 扫描各算法关键参数，测量 CCDF=10^-3 处的 PAPR，为各方法提供参数选择依据。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
params.N_sym = 5000;
rng(42);

N = params.N_fft;
N_sym = params.N_sym;

fprintf('实验 5：参数敏感性分析\n');

% 预生成 OFDM 符号
X_all = zeros(N, N_sym);
x_orig_all = zeros(params.N_os, N_sym);
for i = 1:N_sym
    bits = randi([0 1], N * params.bps, 1);
    X_all(:, i) = ofdm_mod(bits, params.mod_type);
    [x_orig_all(:, i), ~] = ofdm_transmitter(X_all(:, i), params);
end

% 计算 CCDF=10^-3 处的 PAPR的第 99.9 百分位数
papr_metric = @(papr_vec) prctile(papr_vec, 99.9);

% 基准 PAPR
papr_orig = zeros(N_sym, 1);
for i = 1:N_sym
    papr_orig(i) = compute_papr(x_orig_all(:, i));
end
baseline = papr_metric(papr_orig);
fprintf('基准 PAPR@1e-3：%.2f dB\n', baseline);

%% 限幅滤波：扫描限幅比 CR
fprintf('扫描限幅比 CR...\n');
CR_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6];
papr_clip = zeros(length(CR_vals), 1);
for c = 1:length(CR_vals)
    pvals = zeros(N_sym, 1);
    for i = 1:N_sym
        [x_cl, ~] = clipping_filtering(x_orig_all(:,i), params, CR_vals(c));
        pvals(i) = compute_papr(x_cl);
    end
    papr_clip(c) = papr_metric(pvals);
end

fig7 = figure('Position', [100, 100, 600, 400]);
plot(CR_vals, papr_clip, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerFaceColor', colors(2,:));
hold on; yline(baseline, '--k', 'Original'); hold off;
xlabel('限幅比（CR）'); ylabel('PAPR at CCDF=10^{-3}（dB）');
title('限幅滤波：参数敏感性'); grid on;
save_figure(fig7, 'fig07_param_clipping');

%%  μ 律压扩：扫描参数 μ
fprintf('扫描 μ 律参数 mu...\n');
mu_vals = [1, 2, 5, 10, 50, 255];
papr_mu = zeros(length(mu_vals), 1);
for m = 1:length(mu_vals)
    pvals = zeros(N_sym, 1);
    for i = 1:N_sym
        [x_mu, ~] = companding_mu(x_orig_all(:,i), mu_vals(m));
        pvals(i) = compute_papr(x_mu);
    end
    papr_mu(m) = papr_metric(pvals);
end

fig8 = figure('Position', [100, 100, 600, 400]);
semilogx(mu_vals, papr_mu, '-s', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerFaceColor', colors(3,:));
hold on; yline(baseline, '--k', 'Original'); hold off;
xlabel('mu 参数'); ylabel('PAPR at CCDF=10^{-3}（dB）');
title('mu 律压扩：参数敏感性'); grid on;
save_figure(fig8, 'fig08_param_companding');

%% SLM：扫描候选数 U
fprintf('扫描 SLM 候选数 U...\n');
U_vals = [2, 4, 8, 16, 32];
papr_slm = zeros(length(U_vals), 1);
for u = 1:length(U_vals)
    U = U_vals(u);
    P_slm = exp(1j * 2 * pi * rand(N, U));
    P_slm(:,1) = ones(N, 1);
    pvals = zeros(N_sym, 1);
    for i = 1:N_sym
        [~, ~, pvals(i)] = slm(X_all(:,i), params, U, P_slm);
    end
    papr_slm(u) = papr_metric(pvals);
end

fig9 = figure('Position', [100, 100, 600, 400]);
plot(U_vals, papr_slm, '-^', 'Color', colors(5,:), 'LineWidth', 2, 'MarkerFaceColor', colors(5,:));
hold on; yline(baseline, '--k', 'Original'); hold off;
xlabel('候选相位序列数（U）'); ylabel('PAPR at CCDF=10^{-3}（dB）');
title('SLM：参数敏感性'); grid on;
save_figure(fig9, 'fig09_param_slm');

%% PTS：扫描子块数 V 和相位因子集 W
fprintf('扫描 PTS 子块数 V...\n');
V_vals = [2, 4, 8];
W_sets = {[1, -1], [1, -1, 1j, -1j]};
W_labels = {'W={+/-1}', 'W={+/-1, +/-j}'};

fig10 = figure('Position', [100, 100, 600, 400]);
hold on;
for w = 1:length(W_sets)
    papr_pts_v = zeros(length(V_vals), 1);
    for v = 1:length(V_vals)
        pvals = zeros(N_sym, 1);
        for i = 1:N_sym
            [~, ~, pvals(i)] = pts(X_all(:,i), params, V_vals(v), 'interleaved', W_sets{w});
        end
        papr_pts_v(v) = papr_metric(pvals);
    end
    plot(V_vals, papr_pts_v, ['-' markers{w+5}], 'Color', colors(5+w,:), ...
         'LineWidth', 2, 'MarkerFaceColor', colors(5+w,:), 'DisplayName', W_labels{w});
end
yline(baseline, '--k', 'Original');
hold off;
xlabel('子块数（V）'); ylabel('PAPR at CCDF=10^{-3}（dB）');
title('PTS：参数敏感性'); legend('Location', 'northeast'); grid on;
save_figure(fig10, 'fig10_param_pts');

%% 预留音调：扫描预留比例
fprintf('扫描预留音调比例...\n');
ratio_vals = [0.02, 0.05, 0.10, 0.15, 0.20];
papr_tr = zeros(length(ratio_vals), 1);
for r = 1:length(ratio_vals)
    pvals = zeros(N_sym, 1);
    for i = 1:N_sym
        [x_tr, ~, ~] = tone_reservation(X_all(:,i), params, ratio_vals(r), 10);
        pvals(i) = compute_papr(x_tr);
    end
    papr_tr(r) = papr_metric(pvals);
end

fig11 = figure('Position', [100, 100, 600, 400]);
plot(ratio_vals*100, papr_tr, '-d', 'Color', colors(7,:), 'LineWidth', 2, 'MarkerFaceColor', colors(7,:));
hold on; yline(baseline, '--k', 'Original'); hold off;
xlabel('预留子载波比例（%）'); ylabel('PAPR at CCDF=10^{-3}（dB）');
title('预留音调：参数敏感性'); grid on;
save_figure(fig11, 'fig11_param_tr');

save('results/data/param_sweep_results.mat', 'CR_vals', 'papr_clip', ...
     'mu_vals', 'papr_mu', 'U_vals', 'papr_slm', 'V_vals', ...
     'ratio_vals', 'papr_tr', 'baseline', 'params');
fprintf('实验 5 完成。\n');
