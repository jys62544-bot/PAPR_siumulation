% SIM_02_BER_COMPARISON
% 对比所有 PAPR 降低算法在 AWGN 信道下的 BER vs SNR 曲线， 展示各 PAPR 降低方法引入的失真代价。
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);

N = params.N_fft;
snr_range = params.snr_range;
n_snr = length(snr_range);
max_bits = 1e6;       % 每个 SNR 点的最大比特数
min_errors = 100;     % 保证统计可靠性所需的最小误比特数

algo_names = {'Original', 'Clipping (CR=1.2)', 'SLM (U=8)', ...
              'PTS (V=4)', 'Tone Res. (10%)', 'mu-law (mu=10)', ...
              'Golay Coding'};
n_algo = length(algo_names);
ber_results = zeros(n_snr, n_algo);

fprintf('实验 2：BER 对比（AWGN 信道）\n');

% 预生成固定 SLM 相位序列和 PTS 相位因子
P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N, 1);
W_pts = [1, -1, 1j, -1j];
V_pts = 4;

% PTS 分块方案
block_size_pts = N / V_pts;
partition_pts = zeros(block_size_pts, V_pts);
for v = 1:V_pts
    partition_pts(:, v) = v:V_pts:N;
end

for s = 1:n_snr
    snr = snr_range(s);
    fprintf('SNR = %d dB: ', snr);

    bit_errors = zeros(1, n_algo);
    total_bits = 0;
    total_bits_golay = 0;
    bit_errors_golay = 0;

    while total_bits < max_bits
        % 生成一个 OFDM 符号
        bits_tx = randi([0 1], N * params.bps, 1);
        X = ofdm_mod(bits_tx, params.mod_type);
        [x_orig, ~, x_orig_cp] = ofdm_transmitter(X, params);
        bits_per_sym = length(bits_tx);

        % 算法 1：原始信号
        [y, ~] = channel_awgn(x_orig_cp, snr);
        X_hat = ofdm_receiver(y, params);
        bits_rx = ofdm_demod(X_hat, params.mod_type);
        bit_errors(1) = bit_errors(1) + sum(bits_tx ~= bits_rx);

        % 算法 2：限幅滤波
        [x_cl, ~] = clipping_filtering(x_orig, params, 1.2);
        x_cl_cp = cp_add(x_cl, params);
        [y_cl, ~] = channel_awgn(x_cl_cp, snr);
        X_hat_cl = ofdm_receiver(y_cl, params);
        bits_cl = ofdm_demod(X_hat_cl, params.mod_type);
        bit_errors(2) = bit_errors(2) + sum(bits_tx ~= bits_cl);

        % 算法 3：SLM
        [x_slm, u_sel, ~, ~] = slm(X, params, 8, P_slm);
        x_slm_cp = cp_add(x_slm, params);
        [y_slm, ~] = channel_awgn(x_slm_cp, snr);
        X_hat_slm = ofdm_receiver(y_slm, params);
        X_hat_slm = X_hat_slm .* conj(P_slm(:, u_sel));
        bits_slm = ofdm_demod(X_hat_slm, params.mod_type);
        bit_errors(3) = bit_errors(3) + sum(bits_tx ~= bits_slm);

        % 算法 4：PTS
        [x_pts, b_opt, ~] = pts(X, params, V_pts, 'interleaved', W_pts);
        x_pts_cp = cp_add(x_pts, params);
        [y_pts, ~] = channel_awgn(x_pts_cp, snr);
        X_hat_pts = ofdm_receiver(y_pts, params);
        for v = 1:V_pts
            X_hat_pts(partition_pts(:, v)) = X_hat_pts(partition_pts(:, v)) / b_opt(v);
        end
        bits_pts = ofdm_demod(X_hat_pts, params.mod_type);
        bit_errors(4) = bit_errors(4) + sum(bits_tx ~= bits_pts);

        % 算法 5：预留音调
        [x_tr, ~, reserved_idx] = tone_reservation(X, params, 0.10, 10);
        x_tr_cp = cp_add(x_tr, params);
        [y_tr, ~] = channel_awgn(x_tr_cp, snr);
        X_hat_tr = ofdm_receiver(y_tr, params);
        bits_tr = ofdm_demod(X_hat_tr, params.mod_type);
        bit_errors(5) = bit_errors(5) + sum(bits_tx ~= bits_tr);

        % 算法 6：μ 律压扩
        [x_mu, expand_mu] = companding_mu(x_orig, 10);
        x_mu_cp = cp_add(x_mu, params);
        [y_mu, ~] = channel_awgn(x_mu_cp, snr);
        y_mu_exp = expand_mu(y_mu(params.N_cp * params.L + 1 : end));
        X_hat_mu = ofdm_receiver(y_mu_exp, params);
        bits_mu = ofdm_demod(X_hat_mu, params.mod_type);
        bit_errors(6) = bit_errors(6) + sum(bits_tx ~= bits_mu);

        % 算法 7：Golay 编码（编码法自身生成序列，BER 通过加噪后解调对比）
        [x_gl, X_gl, ~] = golay_coding(N, params, 16);
        x_gl_cp = cp_add(x_gl, params);
        bits_gl_tx = ofdm_demod(X_gl, params.mod_type);  % Golay 编码的等效比特
        [y_gl, ~] = channel_awgn(x_gl_cp, snr);
        X_hat_gl = ofdm_receiver(y_gl, params);
        bits_gl_rx = ofdm_demod(X_hat_gl, params.mod_type);
        n_gl_bits = length(bits_gl_tx);
        bit_errors_golay = bit_errors_golay + sum(bits_gl_tx ~= bits_gl_rx);
        total_bits_golay = total_bits_golay + n_gl_bits;

        total_bits = total_bits + bits_per_sym;

        % 所有算法均累积足够误比特数时提前终止
        if all(bit_errors >= min_errors) && bit_errors_golay >= min_errors
            break;
        end
    end

    ber_results(s, 1:6) = bit_errors(1:6) / total_bits;
    ber_results(s, 7) = bit_errors_golay / total_bits_golay;
    fprintf('BER = [');
    fprintf('%.2e ', ber_results(s,:));
    fprintf(']\n');
end

% 保存数据
save('results/data/ber_results.mat', 'ber_results', 'algo_names', 'snr_range', 'params');

% 绘制 BER 曲线
fig = figure('Position', [100, 100, 900, 600]);
ax = axes(fig);
hold(ax, 'on');
for a = 1:n_algo
    semilogy(ax, snr_range, ber_results(:, a), [lstyles{a} markers{a}], ...
        'Color', colors(a,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors(a,:), 'DisplayName', algo_names{a});
end
hold(ax, 'off');

xlabel(ax, 'SNR (dB)');
ylabel(ax, '误比特率');
title(ax, 'BER 性能对比（AWGN 信道）');
ylim(ax, [1e-5 1]);
grid(ax, 'on');

% 图例放在绘图区域外右侧，避免遮挡
lg = legend(ax, 'Location', 'eastoutside');
set(lg, 'FontSize', 9);

save_figure(fig, 'fig02_ber_awgn');
fprintf('实验 2 完成。\n');
