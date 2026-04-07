% SIM_06_DNN_PAPR
% DNN 逐子载波补偿方案：发射端限幅降 PAPR，接收端 DNN 逐子载波补偿失真。
% 支持逐子载波自适应调制。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);
N = params.N_fft;
L = params.L;
N_os = N * L;
CR_dnn = 1.2;

N_train = 20000;
N_test  = 5000;
fprintf('实验 6：基于 DNN 的 PAPR 降低（限幅+逐子载波补偿方案）\n');

%% 自适应调制设置
if params.adaptive
    dummy_x = zeros(N_os, 1);
    [~, H_ch, h_ch] = channel_multipath(dummy_x, params, params.channel_model);
    snr_work = 15;
    [mod_map, ~] = adaptive_modulation(H_ch, snr_work, params.snr_thresholds);
    fprintf('  自适应调制，%s 信道，平均 bps=%.2f\n', params.channel_model, mean(mod_map));
    fprintf('  N_fft=%d, CR=%.1f\n', N, CR_dnn);
else
    fprintf('  N_fft=%d, CR=%.1f, %s\n', N, CR_dnn, params.mod_type);
    mod_map = [];
    H_ch = [];
    h_ch = [];
end

%% 生成训练集和测试集
fprintf('正在生成 %d 个训练符号 + %d 个测试符号...\n', N_train, N_test);
X_train = zeros(N, N_train);
for i = 1:N_train
    if params.adaptive
        [X_train(:, i), ~, ~] = ofdm_mod_adaptive(mod_map);
    else
        bits = randi([0 1], N * params.bps, 1);
        X_train(:, i) = ofdm_mod(bits, params.mod_type);
    end
end

X_test = zeros(N, N_test);
bits_test = cell(N_test, 1);
for i = 1:N_test
    if params.adaptive
        [X_test(:, i), bits_test{i}, ~] = ofdm_mod_adaptive(mod_map);
    else
        bits_test{i} = randi([0 1], N * params.bps, 1);
        X_test(:, i) = ofdm_mod(bits_test{i}, params.mod_type);
    end
end

%% 训练 DNN 逐子载波补偿器
fprintf('训练 DNN 补偿器 (CR=%.1f)...\n', CR_dnn);
opts.CR = CR_dnn;
opts.epochs = 30;
opts.batch_size = 1024;
opts.lr = 1e-3;
[net_trained, train_info] = dnn_papr_reduction(X_train, params, opts, mod_map, H_ch);

%% 图 1：训练损失曲线
fig12 = figure('Position', [100, 100, 600, 400]);
plot(1:opts.epochs, train_info.loss_total, '-b', 'LineWidth', 1.5);
xlabel('轮次'); ylabel('MSE 损失');
title('DNN 逐子载波补偿器训练损失'); grid on;
save_figure(fig12, 'fig12_dnn_training');

%% 预计算限幅信号
fprintf('预计算限幅信号...\n');
x_clip_test = zeros(N_os, N_test);
X_clip_freq = zeros(N, N_test);
papr_orig = zeros(N_test, 1);
papr_clip = zeros(N_test, 1);

for i = 1:N_test
    [x_orig, ~] = ofdm_transmitter(X_test(:,i), params);
    papr_orig(i) = compute_papr(x_orig);
    A = CR_dnn * rms(x_orig);
    x_cl = min(abs(x_orig), A) .* exp(1j * angle(x_orig));
    x_clip_test(:, i) = x_cl;
    papr_clip(i) = compute_papr(x_cl);
    Y_full = fft(x_cl, N_os) / sqrt(N_os);
    X_clip_freq(:, i) = [Y_full(1:N/2); Y_full(N_os - N/2 + 1:end)];
end

%% CCDF 对比
fprintf('计算 CCDF 对比...\n');
P_slm = exp(1j * 2 * pi * rand(N, 16));
P_slm(:,1) = ones(N, 1);
W_pts = [1, -1, 1j, -1j];

papr_slm   = zeros(N_test, 1);
papr_pts   = zeros(N_test, 1);
papr_golay = zeros(N_test, 1);
for i = 1:N_test
    [~, ~, papr_slm(i)] = slm(X_test(:,i), params, 16, P_slm);
    [~, ~, papr_pts(i)] = pts(X_test(:,i), params, 4, 'interleaved', W_pts);
    [~, ~, gl_info6] = golay_coding(X_test(:,i), params, 16);
    papr_golay(i) = gl_info6.papr;
end

fprintf('PAPR 均值(dB): Original=%.2f, Clip+DNN=%.2f, SLM=%.2f, PTS=%.2f, Golay=%.2f\n', ...
    mean(papr_orig), mean(papr_clip), mean(papr_slm), mean(papr_pts), mean(papr_golay));

fig13 = figure('Position', [100, 100, 900, 600]);
ax13 = axes(fig13);
hold(ax13, 'on');
[pa, cc] = compute_ccdf(papr_orig);
plot(ax13, pa, cc, '-k', 'LineWidth', 1.5, 'DisplayName', 'Original');
[pa, cc] = compute_ccdf(papr_slm);
plot(ax13, pa, cc, '--', 'Color', colors(3,:), 'LineWidth', 1.5, 'DisplayName', 'SLM (U=16)');
[pa, cc] = compute_ccdf(papr_pts);
plot(ax13, pa, cc, '-.', 'Color', colors(4,:), 'LineWidth', 1.5, 'DisplayName', 'PTS (V=4)');
[pa, cc] = compute_ccdf(papr_golay);
plot(ax13, pa, cc, '-', 'Color', colors(7,:), 'LineWidth', 1.5, 'DisplayName', 'Golay (N=16)');
[pa, cc] = compute_ccdf(papr_clip);
plot(ax13, pa, cc, '-', 'Color', colors(8,:), 'LineWidth', 2, 'DisplayName', ...
    sprintf('Clip(CR=%.1f)+DNN', CR_dnn));
hold(ax13, 'off');
set(ax13, 'YScale', 'log');
xlabel(ax13, 'PAPR_0 (dB)'); ylabel(ax13, 'Pr(PAPR > PAPR_0)');
if params.adaptive
    title(ax13, sprintf('CCDF：DNN 与经典方法对比（自适应调制，%s）', params.channel_model));
else
    title(ax13, 'CCDF：DNN 与经典方法对比');
end
xlim(ax13, [2 14]); ylim(ax13, [1e-3 1]); grid(ax13, 'on');
lg13 = legend(ax13, 'Location', 'eastoutside');
set(lg13, 'FontSize', 9);
save_figure(fig13, 'fig13_dnn_ccdf');

%% BER 对比（Original vs 纯限幅 vs 限幅+DNN）
fprintf('计算 BER 对比...\n');
snr_range = 0:2:30;
n_snr = length(snr_range);
ber_orig = zeros(n_snr, 1);
ber_clip_only = zeros(n_snr, 1);
ber_clip_dnn = zeros(n_snr, 1);

dnn_net = net_trained.net;

% 预计算推理用辅助特征
if params.adaptive
    mod_feat = single(mod_map / 8);
    h_feat = single(abs(H_ch) / max(abs(H_ch)));
else
    mod_feat = single(ones(N, 1) * params.bps / 8);
    h_feat = single(ones(N, 1));
end

for s = 1:n_snr
    snr = snr_range(s);
    errs_orig = 0; errs_clip = 0; errs_dnn = 0; total_bits = 0;

    % 固定装载：使用 snr_work 分配的 mod_map，不随 SNR 扫描变化（非动态 AMC）
    for i = 1:min(N_test, 2000)
        if params.adaptive
            [X_i, bits_tx, ~] = ofdm_mod_adaptive(mod_map);
        else
            bits_tx = bits_test{i};
            X_i = X_test(:, i);
        end
        [x_orig, ~, x_orig_cp] = ofdm_transmitter(X_i, params);

        % 限幅
        A = CR_dnn * rms(x_orig);
        x_cl = min(abs(x_orig), A) .* exp(1j * angle(x_orig));
        x_cl_cp = cp_add(x_cl, params);

        % 统一噪声注入：以原始信号功率为参考
        ref_power = mean(abs(x_orig_cp).^2);

        % 通过信道
        if params.adaptive
            y_mp_orig = conv(x_orig_cp, h_ch);
            y_mp_orig = y_mp_orig(1:length(x_orig_cp));
            [y, ~] = channel_awgn_fixed(y_mp_orig, snr, ref_power);
            y_mp_cl = conv(x_cl_cp, h_ch);
            y_mp_cl = y_mp_cl(1:length(x_cl_cp));
            [y_cl, ~] = channel_awgn_fixed(y_mp_cl, snr, ref_power);
            H_eq = H_ch;
        else
            [y, ~] = channel_awgn_fixed(x_orig_cp, snr, ref_power);
            [y_cl, ~] = channel_awgn_fixed(x_cl_cp, snr, ref_power);
            H_eq = [];
        end

        % Original
        X_hat = ofdm_receiver(y, params, H_eq);
        if params.adaptive
            bits_rx = ofdm_demod_adaptive(X_hat, mod_map);
        else
            bits_rx = ofdm_demod(X_hat, params.mod_type);
        end
        errs_orig = errs_orig + sum(bits_tx ~= bits_rx);

        % 纯限幅
        X_hat_cl = ofdm_receiver(y_cl, params, H_eq);
        if params.adaptive
            bits_cl = ofdm_demod_adaptive(X_hat_cl, mod_map);
        else
            bits_cl = ofdm_demod(X_hat_cl, params.mod_type);
        end
        errs_clip = errs_clip + sum(bits_tx ~= bits_cl);

        % 限幅 + DNN 逐子载波补偿
        X_in_sc = single([real(X_hat_cl).'; imag(X_hat_cl).'; mod_feat.'; h_feat.']);
        X_in_dl = dlarray(X_in_sc, 'CB');
        X_pred_dl = predict(dnn_net, X_in_dl);
        X_pred = double(extractdata(X_pred_dl));
        X_rec = X_pred(1,:).' + 1j * X_pred(2,:).';
        if params.adaptive
            bits_dnn = ofdm_demod_adaptive(X_rec, mod_map);
        else
            bits_dnn = ofdm_demod(X_rec, params.mod_type);
        end
        errs_dnn = errs_dnn + sum(bits_tx ~= bits_dnn);

        total_bits = total_bits + length(bits_tx);
    end

    ber_orig(s)      = errs_orig / total_bits;
    ber_clip_only(s) = errs_clip / total_bits;
    ber_clip_dnn(s)  = errs_dnn / total_bits;
    fprintf('  SNR=%2d dB: Original=%.2e, Clip=%.2e, Clip+DNN=%.2e\n', ...
        snr, ber_orig(s), ber_clip_only(s), ber_clip_dnn(s));
end

fig14 = figure('Position', [100, 100, 800, 500]);
ax14 = axes(fig14);
semilogy(ax14, snr_range, ber_orig, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Original (无 PAPR 降低)');
hold(ax14, 'on');
semilogy(ax14, snr_range, ber_clip_only, '--^', 'Color', colors(2,:), 'LineWidth', 1.5, ...
    'MarkerFaceColor', colors(2,:), 'DisplayName', sprintf('Clipping (CR=%.1f)', CR_dnn));
semilogy(ax14, snr_range, ber_clip_dnn, '-s', 'Color', colors(8,:), 'LineWidth', 2, ...
    'MarkerFaceColor', colors(8,:), 'DisplayName', sprintf('Clip(CR=%.1f)+DNN', CR_dnn));
hold(ax14, 'off');
set(ax14, 'YScale', 'log');
xlabel(ax14, 'SNR (dB)'); ylabel(ax14, '误比特率');
if params.adaptive
    title(ax14, sprintf('BER：DNN 补偿 vs 纯限幅（自适应调制，%s）', params.channel_model));
else
    title(ax14, 'BER：DNN 补偿 vs 纯限幅 vs 无限幅');
end
ylim(ax14, [1e-5 1]); grid(ax14, 'on');
lg14 = legend(ax14, 'Location', 'eastoutside');
set(lg14, 'FontSize', 9);
save_figure(fig14, 'fig14_dnn_ber');

%% 星座图（SNR=20dB）
fprintf('绘制星座图...\n');
snr_const = 20;

if params.adaptive
    [X_1, ~, ~] = ofdm_mod_adaptive(mod_map);
else
    X_1 = X_test(:, 1);
end

[~, ~, x_orig_cp] = ofdm_transmitter(X_1, params);
A_1 = CR_dnn * rms(x_orig_cp(params.N_cp*L+1:end));
x_cl_1 = min(abs(x_orig_cp(params.N_cp*L+1:end)), A_1) .* exp(1j * angle(x_orig_cp(params.N_cp*L+1:end)));
x_cl_1_cp = cp_add(x_cl_1, params);

ref_power_const = mean(abs(x_orig_cp).^2);

if params.adaptive
    y_mp = conv(x_orig_cp, h_ch); y_mp = y_mp(1:length(x_orig_cp));
    [y_orig, ~] = channel_awgn_fixed(y_mp, snr_const, ref_power_const);
    y_mp_cl = conv(x_cl_1_cp, h_ch); y_mp_cl = y_mp_cl(1:length(x_cl_1_cp));
    [y_cl, ~] = channel_awgn_fixed(y_mp_cl, snr_const, ref_power_const);
    H_eq = H_ch;
else
    [y_orig, ~] = channel_awgn_fixed(x_orig_cp, snr_const, ref_power_const);
    [y_cl, ~] = channel_awgn_fixed(x_cl_1_cp, snr_const, ref_power_const);
    H_eq = [];
end

X_hat_orig = ofdm_receiver(y_orig, params, H_eq);
X_hat_cl_1 = ofdm_receiver(y_cl, params, H_eq);

X_in_sc = single([real(X_hat_cl_1).'; imag(X_hat_cl_1).'; mod_feat.'; h_feat.']);
X_in_dl = dlarray(X_in_sc, 'CB');
X_pred_dl = predict(dnn_net, X_in_dl);
X_pred = double(extractdata(X_pred_dl));
X_rec_1 = X_pred(1,:).' + 1j * X_pred(2,:).';

% 自动计算星座图范围
all_syms = [X_hat_orig; X_hat_cl_1; X_rec_1];
ax_lim = max(abs([real(all_syms); imag(all_syms)])) * 1.15;

fig15 = figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
scatter(real(X_hat_orig), imag(X_hat_orig), 10, 'b', 'filled');
title(sprintf('Original (SNR=%ddB)', snr_const)); xlabel('I'); ylabel('Q');
axis equal; grid on; xlim([-ax_lim ax_lim]); ylim([-ax_lim ax_lim]);

subplot(1,3,2);
scatter(real(X_hat_cl_1), imag(X_hat_cl_1), 10, 'r', 'filled');
title(sprintf('Clipping CR=%.1f (SNR=%ddB)', CR_dnn, snr_const)); xlabel('I'); ylabel('Q');
axis equal; grid on; xlim([-ax_lim ax_lim]); ylim([-ax_lim ax_lim]);

subplot(1,3,3);
scatter(real(X_rec_1), imag(X_rec_1), 10, [0 0.6 0], 'filled');
title(sprintf('Clip+DNN (SNR=%ddB)', snr_const)); xlabel('I'); ylabel('Q');
axis equal; grid on; xlim([-ax_lim ax_lim]); ylim([-ax_lim ax_lim]);
save_figure(fig15, 'fig15_dnn_constellation');

%% 保存
save('results/data/dnn_results.mat', 'train_info', 'papr_orig', 'papr_clip', ...
    'papr_golay', 'ber_orig', 'ber_clip_only', 'ber_clip_dnn', 'snr_range', ...
    'params', 'CR_dnn', 'net_trained');
fprintf('实验 6 完成。\n');
