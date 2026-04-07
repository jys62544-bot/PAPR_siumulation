% SIM_02_BER_COMPARISON
% 对比 PAPR 降低算法在多径信道下的 BER vs SNR 曲线。
% 支持逐子载波自适应调制（固定装载方案）。
%
% 公平比较原则：
%   1) 所有算法使用同一组用户数据比特
%   2) 噪声功率统一基于原始信号功率和目标 SNR 计算（非各算法自身功率）
%   3) TR 只统计数据子载波比特（预留载波的数据率损失已体现）
%   4) Golay 用 Rate-1/2 编码：取 X 前 N/2 符号经 Golay 对编码到 N 个子载波；
%      BER 只统计前 N/2 符号对应的信息比特，数据率减半是方法固有代价
%   5) SLM/PTS 假设 side information 完美已知（学术标准做法）
clear; clc; close all;
addpath('core', 'analysis', 'papr_reduction');
[colors, markers, lstyles] = plot_config();
params = get_default_params();
rng(42);

N = params.N_fft;
snr_range = params.snr_range;
n_snr = length(snr_range);
max_bits = 1e6;
min_errors = 100;

% 所有算法使用同一数据块
algo_names = {'Original', 'Clipping (CR=1.2)', 'SLM (U=8)', ...
              'PTS (V=4)', 'Tone Res. (10%)', 'mu-law (mu=10)', 'Golay (N=16)'};
n_algo = length(algo_names);
ber_results = zeros(n_snr, n_algo);

fprintf('实验 2：BER 对比\n');

% 预生成固定 SLM 相位序列和 PTS 相位因子
P_slm = exp(1j * 2 * pi * rand(N, 8));
P_slm(:,1) = ones(N, 1);
W_pts = [1, -1, 1j, -1j];
V_pts = 4;

block_size_pts = N / V_pts;
partition_pts = zeros(block_size_pts, V_pts);
for v = 1:V_pts
    partition_pts(:, v) = v:V_pts:N;
end

%% 生成信道（固定信道实现，所有 SNR 共用）
if params.adaptive
    dummy_x = zeros(params.N_os + params.N_cp * params.L, 1);
    [~, H_ch, h_ch] = channel_multipath(dummy_x, params, params.channel_model);
    % 固定装载：用 snr_work 分配调制方式，所有 SNR 点共用（非动态 AMC）
    snr_work = 15;
    [mod_map, ~] = adaptive_modulation(H_ch, snr_work, params.snr_thresholds);
    bits_per_sym = sum(mod_map);
    fprintf('使用 %s 多径信道 + 固定装载 (工作 SNR=%d dB, avg bps=%.2f)\n', ...
        params.channel_model, snr_work, mean(mod_map));
else
    H_ch = [];
    h_ch = [];
    mod_map = [];
    bits_per_sym = N * params.bps;
    snr_work = [];
end

for s = 1:n_snr
    snr = snr_range(s);
    fprintf('SNR = %d dB: ', snr);

    bit_errors = zeros(1, n_algo);
    total_bits = 0;
    total_bits_tr = 0;
    total_bits_golay = 0;

    while total_bits < max_bits
        % 生成一个 OFDM 符号
        if params.adaptive
            [X, bits_tx, ~] = ofdm_mod_adaptive(mod_map);
        else
            bits_tx = randi([0 1], N * params.bps, 1);
            X = ofdm_mod(bits_tx, params.mod_type);
        end
        [x_orig, ~, x_orig_cp] = ofdm_transmitter(X, params);

        % 计算原始信号的参考功率（用于统一噪声注入）
        ref_power = mean(abs(x_orig_cp).^2);

        % ========== 算法 1：原始信号 ==========
        if params.adaptive
            y_mp = conv(x_orig_cp, h_ch);
            y_mp = y_mp(1:length(x_orig_cp));
            [y, ~] = channel_awgn_fixed(y_mp, snr, ref_power);
            H_eq = H_ch;
        else
            [y, ~] = channel_awgn_fixed(x_orig_cp, snr, ref_power);
            H_eq = [];
        end
        X_hat = ofdm_receiver(y, params, H_eq);
        if params.adaptive
            bits_rx = ofdm_demod_adaptive(X_hat, mod_map);
        else
            bits_rx = ofdm_demod(X_hat, params.mod_type);
        end
        bit_errors(1) = bit_errors(1) + sum(bits_tx ~= bits_rx);

        % ========== 算法 2：限幅滤波 ==========
        [x_cl, ~] = clipping_filtering(x_orig, params, 1.2);
        x_cl_cp = cp_add(x_cl, params);
        if params.adaptive
            y_mp_cl = conv(x_cl_cp, h_ch);
            y_mp_cl = y_mp_cl(1:length(x_cl_cp));
            [y_cl, ~] = channel_awgn_fixed(y_mp_cl, snr, ref_power);
        else
            [y_cl, ~] = channel_awgn_fixed(x_cl_cp, snr, ref_power);
        end
        X_hat_cl = ofdm_receiver(y_cl, params, H_eq);
        if params.adaptive
            bits_cl = ofdm_demod_adaptive(X_hat_cl, mod_map);
        else
            bits_cl = ofdm_demod(X_hat_cl, params.mod_type);
        end
        bit_errors(2) = bit_errors(2) + sum(bits_tx ~= bits_cl);

        % ========== 算法 3：SLM（side information 假设完美已知）==========
        [x_slm, u_sel, ~, ~] = slm(X, params, 8, P_slm);
        x_slm_cp = cp_add(x_slm, params);
        if params.adaptive
            y_mp_slm = conv(x_slm_cp, h_ch);
            y_mp_slm = y_mp_slm(1:length(x_slm_cp));
            [y_slm, ~] = channel_awgn_fixed(y_mp_slm, snr, ref_power);
        else
            [y_slm, ~] = channel_awgn_fixed(x_slm_cp, snr, ref_power);
        end
        X_hat_slm = ofdm_receiver(y_slm, params, H_eq);
        X_hat_slm = X_hat_slm .* conj(P_slm(:, u_sel));
        if params.adaptive
            bits_slm = ofdm_demod_adaptive(X_hat_slm, mod_map);
        else
            bits_slm = ofdm_demod(X_hat_slm, params.mod_type);
        end
        bit_errors(3) = bit_errors(3) + sum(bits_tx ~= bits_slm);

        % ========== 算法 4：PTS（side information 假设完美已知）==========
        [x_pts, b_opt, ~] = pts(X, params, V_pts, 'interleaved', W_pts);
        x_pts_cp = cp_add(x_pts, params);
        if params.adaptive
            y_mp_pts = conv(x_pts_cp, h_ch);
            y_mp_pts = y_mp_pts(1:length(x_pts_cp));
            [y_pts, ~] = channel_awgn_fixed(y_mp_pts, snr, ref_power);
        else
            [y_pts, ~] = channel_awgn_fixed(x_pts_cp, snr, ref_power);
        end
        X_hat_pts = ofdm_receiver(y_pts, params, H_eq);
        for v = 1:V_pts
            X_hat_pts(partition_pts(:, v)) = X_hat_pts(partition_pts(:, v)) / b_opt(v);
        end
        if params.adaptive
            bits_pts = ofdm_demod_adaptive(X_hat_pts, mod_map);
        else
            bits_pts = ofdm_demod(X_hat_pts, params.mod_type);
        end
        bit_errors(4) = bit_errors(4) + sum(bits_tx ~= bits_pts);

        % ========== 算法 5：预留音调（数据率损失 ~10%）==========
        [x_tr, ~, reserved_idx] = tone_reservation(X, params, 0.10, 30);
        x_tr_cp = cp_add(x_tr, params);
        if params.adaptive
            y_mp_tr = conv(x_tr_cp, h_ch);
            y_mp_tr = y_mp_tr(1:length(x_tr_cp));
            [y_tr, ~] = channel_awgn_fixed(y_mp_tr, snr, ref_power);
        else
            [y_tr, ~] = channel_awgn_fixed(x_tr_cp, snr, ref_power);
        end
        X_hat_tr = ofdm_receiver(y_tr, params, H_eq);
        if params.adaptive
            bits_tr = ofdm_demod_adaptive(X_hat_tr, mod_map);
            % 只统计数据子载波上的比特（排除预留子载波）
            data_idx_tr = setdiff(1:N, reserved_idx);
            data_mask_tr = false(length(bits_tx), 1);
            bit_ptr = 0;
            for kk = 1:N
                bps_k = mod_map(kk);
                if ismember(kk, data_idx_tr)
                    data_mask_tr(bit_ptr+1 : bit_ptr+bps_k) = true;
                end
                bit_ptr = bit_ptr + bps_k;
            end
            bit_errors(5) = bit_errors(5) + sum(bits_tx(data_mask_tr) ~= bits_tr(data_mask_tr));
        else
            bits_tr = ofdm_demod(X_hat_tr, params.mod_type);
            data_idx_tr = setdiff(1:N, reserved_idx);
            bps_fixed = params.bps;
            data_mask_tr = false(length(bits_tx), 1);
            for kk = data_idx_tr
                data_mask_tr((kk-1)*bps_fixed+1 : kk*bps_fixed) = true;
            end
            bit_errors(5) = bit_errors(5) + sum(bits_tx(data_mask_tr) ~= bits_tr(data_mask_tr));
        end
        total_bits_tr = total_bits_tr + sum(data_mask_tr);

        % ========== 算法 6：μ 律压扩 ==========
        [x_mu, expand_mu] = companding_mu(x_orig, 10);
        x_mu_cp = cp_add(x_mu, params);
        if params.adaptive
            y_mp_mu = conv(x_mu_cp, h_ch);
            y_mp_mu = y_mp_mu(1:length(x_mu_cp));
            [y_mu, ~] = channel_awgn_fixed(y_mp_mu, snr, ref_power);
            % 多径下：去CP → FFT → 信道均衡 → IFFT → 逆压扩 → FFT → 解调
            y_mu_nocp = y_mu(params.N_cp * params.L + 1 : end);
            Y_mu_full = fft(y_mu_nocp, params.N_os) / sqrt(params.N_os);
            X_mu_eq = [Y_mu_full(1:N/2); Y_mu_full(params.N_os - N/2 + 1:end)];
            X_mu_eq = X_mu_eq ./ H_eq;
            X_mu_eq_pad = [X_mu_eq(1:N/2); zeros((params.L-1)*N, 1); X_mu_eq(N/2+1:end)];
            x_mu_td = sqrt(params.N_os) * ifft(X_mu_eq_pad, params.N_os);
            x_mu_expanded = expand_mu(x_mu_td);
            X_mu_final_full = fft(x_mu_expanded, params.N_os) / sqrt(params.N_os);
            X_hat_mu = [X_mu_final_full(1:N/2); X_mu_final_full(params.N_os - N/2 + 1:end)];
        else
            [y_mu, ~] = channel_awgn_fixed(x_mu_cp, snr, ref_power);
            y_mu_exp = expand_mu(y_mu(params.N_cp * params.L + 1 : end));
            X_hat_mu = ofdm_receiver(y_mu_exp, params, H_eq);
        end
        if params.adaptive
            bits_mu = ofdm_demod_adaptive(X_hat_mu, mod_map);
        else
            bits_mu = ofdm_demod(X_hat_mu, params.mod_type);
        end
        bit_errors(6) = bit_errors(6) + sum(bits_tx ~= bits_mu);

        % ========== 算法 7：Golay 互补序列编码（Rate-1/2）==========
        % 编码：用同一 X 的前 N/2 符号作为信息数据，经 Golay 对编码到 N 个子载波
        [x_gl, ~, gl_info] = golay_coding(X, params, 16);
        x_gl_cp = cp_add(x_gl, params);
        if params.adaptive
            y_mp_gl = conv(x_gl_cp, h_ch);
            y_mp_gl = y_mp_gl(1:length(x_gl_cp));
            [y_gl, ~] = channel_awgn_fixed(y_mp_gl, snr, ref_power);
        else
            [y_gl, ~] = channel_awgn_fixed(x_gl_cp, snr, ref_power);
        end
        X_hat_gl = ofdm_receiver(y_gl, params, H_eq);
        % 解码：前 N/2 子载波的接收符号除以 a 序列（a 是 ±1，所以乘即可）
        X_hat_gl_data = X_hat_gl(1:N/2) .* gl_info.a_seq;
        % 只统计前 N/2 符号对应的比特（信息部分）
        if params.adaptive
            % 自适应模式：只取前 N/2 子载波的 mod_map
            bits_gl = ofdm_demod_adaptive(X_hat_gl_data, mod_map(1:N/2));
            % 信息比特 = bits_tx 中前 N/2 子载波对应的比特
            n_info_bits = sum(mod_map(1:N/2));
            bits_gl_tx_part = bits_tx(1:n_info_bits);
        else
            bits_gl = ofdm_demod(X_hat_gl_data, params.mod_type);
            n_info_bits = N/2 * params.bps;
            bits_gl_tx_part = bits_tx(1:n_info_bits);
        end
        bit_errors(7) = bit_errors(7) + sum(bits_gl_tx_part ~= bits_gl);
        total_bits_golay = total_bits_golay + n_info_bits;

        total_bits = total_bits + bits_per_sym;

        % 提前停止
        if all(bit_errors(1:n_algo) >= min_errors)
            break;
        end
    end

    ber_results(s, 1:4) = bit_errors(1:4) / total_bits;
    ber_results(s, 5) = bit_errors(5) / max(total_bits_tr, 1);
    ber_results(s, 6) = bit_errors(6) / total_bits;
    ber_results(s, 7) = bit_errors(7) / max(total_bits_golay, 1);
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
set(ax, 'YScale', 'log');

xlabel(ax, 'SNR (dB)');
ylabel(ax, '误比特率');
if params.adaptive
    title(ax, sprintf('BER 对比（固定装载，%s 信道，工作 SNR=%d dB）', ...
        params.channel_model, snr_work));
else
    title(ax, 'BER 性能对比（AWGN 信道）');
end
ylim(ax, [1e-5 1]);
grid(ax, 'on');

lg = legend(ax, 'Location', 'eastoutside');
set(lg, 'FontSize', 9);

save_figure(fig, 'fig02_ber_comparison');
fprintf('实验 2 完成。\n');
