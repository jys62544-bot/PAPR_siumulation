% SIM_06_DNN_PAPR
% 训练 DNN 自编码器用于 PAPR 降低，并与经典方法对比。用的MATLAB Deep Learning Toolbox。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);
N = params.N_fft;
N_train = 50000;
N_test = 5000;
fprintf('实验 6：基于 DNN 的 PAPR 降低\n');

%% 生成训练集和测试集
fprintf('正在生成 %d 个训练符号 + %d 个测试符号...\n', N_train, N_test);
X_train = zeros(N, N_train);
for i = 1:N_train
    bits = randi([0 1], N * params.bps, 1);
    X_train(:, i) = ofdm_mod(bits, params.mod_type);
end

X_test = zeros(N, N_test);
bits_test = cell(N_test, 1);
for i = 1:N_test
    bits_test{i} = randi([0 1], N * params.bps, 1);
    X_test(:, i) = ofdm_mod(bits_test{i}, params.mod_type);
end

%% 训练网络
fprintf('训练 DNN 自编码器...\n');
opts.epochs = 50;
opts.batch_size = 128;
opts.alpha = 0.1;
opts.lambda_mse = 1.0;
opts.lambda_papr = 0.01;
opts.lr = 1e-3;
[~, ~, net_trained, train_info] = dnn_papr_reduction(X_train, params, opts);

%% 图 1：训练损失曲线
fig12 = figure('Position', [100, 100, 700, 400]);
subplot(1,2,1);
plot(1:opts.epochs, train_info.loss_total, '-b', 'LineWidth', 1.5);
xlabel('轮次'); ylabel('总损失');
title('训练损失'); grid on;

subplot(1,2,2);
yyaxis left;
plot(1:opts.epochs, train_info.loss_mse, '-b', 'LineWidth', 1.5);
ylabel('MSE 损失');
yyaxis right;
plot(1:opts.epochs, train_info.loss_papr, '-r', 'LineWidth', 1.5);
ylabel('PAPR 损失');
xlabel('轮次'); title('各损失分量'); grid on;
legend('MSE', 'PAPR', 'Location', 'northeast');
save_figure(fig12, 'fig12_dnn_training');

%% 对测试集应用训练好的网络
fprintf('将 DNN 应用于测试集...\n');
opts_apply.mode = 'apply';
opts_apply.net = net_trained;
opts_apply.alpha = opts.alpha;
[x_dnn_test, X_dnn_test] = dnn_papr_reduction(X_test, params, opts_apply);

%% 图 2：CCDF 对比（DNN vs 经典方法）
fprintf('计算 CCDF 对比...\n');
papr_orig_test = zeros(N_test, 1);
papr_dnn_test = zeros(N_test, 1);
papr_slm_test = zeros(N_test, 1);
papr_pts_test = zeros(N_test, 1);

P_slm = exp(1j * 2 * pi * rand(N, 16));
P_slm(:,1) = ones(N, 1);
W_pts = [1, -1, 1j, -1j];

for i = 1:N_test
    [x_orig, ~] = ofdm_transmitter(X_test(:, i), params);
    papr_orig_test(i) = compute_papr(x_orig);
    papr_dnn_test(i)  = compute_papr(x_dnn_test(:, i));
    [~, ~, papr_slm_test(i)] = slm(X_test(:,i), params, 16, P_slm);
    [~, ~, papr_pts_test(i)] = pts(X_test(:,i), params, 4, 'interleaved', W_pts);
end

fig13 = figure('Position', [100, 100, 800, 600]);
hold on;
[pa, cc] = compute_ccdf(papr_orig_test);
plot(pa, cc, '-k', 'LineWidth', 1.5, 'DisplayName', '原始信号');
[pa, cc] = compute_ccdf(papr_slm_test);
plot(pa, cc, '--', 'Color', colors(5,:), 'LineWidth', 1.5, 'DisplayName', 'SLM (U=16)');
[pa, cc] = compute_ccdf(papr_pts_test);
plot(pa, cc, '-.', 'Color', colors(6,:), 'LineWidth', 1.5, 'DisplayName', 'PTS (V=4)');
[pa, cc] = compute_ccdf(papr_dnn_test);
plot(pa, cc, '-', 'Color', colors(8,:), 'LineWidth', 2, 'DisplayName', 'DNN');
hold off;

set(gca, 'YScale', 'log');
xlabel('PAPR_0 (dB)'); ylabel('Pr(PAPR > PAPR_0)');
title('CCDF：DNN 与经典方法对比');
legend('Location', 'southwest');
xlim([4 14]); ylim([1e-3 1]); grid on;
save_figure(fig13, 'fig13_dnn_ccdf');

%% 图 3：BER 对比
fprintf('计算 BER 对比...\n');
snr_range = 0:2:30;
n_snr = length(snr_range);
ber_orig = zeros(n_snr, 1);
ber_dnn = zeros(n_snr, 1);

for s = 1:n_snr
    snr = snr_range(s);
    errs_orig = 0; errs_dnn = 0; total_bits = 0;
    for i = 1:min(N_test, 2000)
        bits_tx = bits_test{i};
        [x_orig, ~] = ofdm_transmitter(X_test(:,i), params);

        % 原始信号
        [y, ~] = channel_awgn(x_orig, snr);
        X_hat = ofdm_receiver(y, params);
        bits_rx = ofdm_demod(X_hat, params.mod_type);
        errs_orig = errs_orig + sum(bits_tx ~= bits_rx);

        % DNN：通过解码器网络恢复原始符号
        [y_dnn, ~] = channel_awgn(x_dnn_test(:,i), snr);
        X_hat_dnn = ofdm_receiver(y_dnn, params);
        X_hat_real = [real(X_hat_dnn); imag(X_hat_dnn)];
        X_hat_dl = dlarray(X_hat_real, 'CB');
        X_rec_dl = predict(net_trained.decoder, X_hat_dl);
        X_rec_data = extractdata(X_rec_dl);
        X_rec = X_rec_data(1:N) + 1j * X_rec_data(N+1:end);
        bits_dnn = ofdm_demod(X_rec, params.mod_type);
        errs_dnn = errs_dnn + sum(bits_tx ~= bits_dnn);

        total_bits = total_bits + length(bits_tx);
    end
    ber_orig(s) = errs_orig / total_bits;
    ber_dnn(s)  = errs_dnn / total_bits;
end

fig14 = figure('Position', [100, 100, 700, 500]);
semilogy(snr_range, ber_orig, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
    'DisplayName', '原始信号');
hold on;
semilogy(snr_range, ber_dnn, '-s', 'Color', colors(8,:), 'LineWidth', 2, ...
    'MarkerFaceColor', colors(8,:), 'DisplayName', 'DNN');
hold off;
xlabel('SNR (dB)'); ylabel('误比特率');
title('BER：DNN 与原始信号对比');
legend('Location', 'southwest');
ylim([1e-5 1]); grid on;
save_figure(fig14, 'fig14_dnn_ber');

%% 图 4：星座图
fig15 = figure('Position', [100, 100, 900, 400]);
subplot(1,2,1);
scatter(real(X_test(:,1)), imag(X_test(:,1)), 10, 'b', 'filled');
title('原始 16QAM 星座'); xlabel('同相分量'); ylabel('正交分量');
axis equal; grid on; xlim([-1.2 1.2]); ylim([-1.2 1.2]);

subplot(1,2,2);
scatter(real(X_dnn_test(:,1)), imag(X_dnn_test(:,1)), 10, 'r', 'filled');
title('DNN 修改后的 16QAM 星座'); xlabel('同相分量'); ylabel('正交分量');
axis equal; grid on; xlim([-1.2 1.2]); ylim([-1.2 1.2]);
save_figure(fig15, 'fig15_dnn_constellation');

save('results/data/dnn_results.mat', 'train_info', 'papr_orig_test', ...
     'papr_dnn_test', 'ber_orig', 'ber_dnn', 'snr_range', 'params');
fprintf('实验 6 完成。\n');
