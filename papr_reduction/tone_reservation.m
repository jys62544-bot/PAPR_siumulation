function [x_reduced, C, reserved_idx] = tone_reservation(X, params, reserved_ratio, max_iter)
% 预留子载波法降低 PAPR，保留一部分子载波用于承载峰值抵消信号，
% 采用经典的基于限幅的梯度投影算法（Tellado 算法）进行优化。
%
% 预留子载波上的原始数据会被置零，因此有效数据率降低为
% (1 - reserved_ratio) 倍。在 BER 统计时应排除预留子载波上的比特。
%   输入：
%     X              - 频域数据（N_fft x 1）
%     params         - get_default_params() 返回的参数结构体
%     reserved_ratio - 预留子载波比例
%     max_iter       - 最大迭代次数
%   输出：
%     x_reduced    - PAPR 降低后的时域信号（N_fft*L x 1）
%     C            - 频域峰值抵消信号（N_fft x 1）
%     reserved_idx - 预留子载波的索引

    N = params.N_fft;
    L = params.L;
    N_os = N * L;

    N_reserved = round(N * reserved_ratio);
    reserved_idx = round(linspace(1, N, N_reserved + 2));
    reserved_idx = reserved_idx(2:end-1);

    data_idx = setdiff(1:N, reserved_idx);

    % 将预留子载波上的数据置零，仅保留数据子载波
    X_data = X;
    X_data(reserved_idx) = 0;

    % 过采样 IFFT 得到初始时域信号
    X_padded = [X_data(1:N/2); zeros((L-1)*N, 1); X_data(N/2+1:end)];
    x_current = sqrt(N_os) * ifft(X_padded, N_os);

    C = zeros(N, 1);

    % 构造预留子载波的核函数（峰值集中的脉冲）
    % kernel: 仅在预留子载波位置有能量的时域脉冲
    p_freq = zeros(N, 1);
    p_freq(reserved_idx) = 1;
    P_padded = [p_freq(1:N/2); zeros((L-1)*N, 1); p_freq(N/2+1:end)];
    p_time = sqrt(N_os) * ifft(P_padded, N_os);
    % 将核函数归一化，使其峰值为 1
    [~, n_p] = max(abs(p_time));
    p_time = p_time / p_time(n_p);

    % 迭代限幅-投影降低峰值
    % 目标门限：基于初始信号的 RMS
    init_rms = rms(x_current);
    target_amp = 3.0 * init_rms;  % 目标幅度门限

    mu_step = 0.5;  % 步长因子，防止过冲

    for iter = 1:max_iter
        peak_amp = max(abs(x_current));
        if peak_amp <= target_amp
            break;
        end

        % 找到最大峰值位置
        [~, n_peak] = max(abs(x_current));
        peak_val = x_current(n_peak);

        % 需要在该位置抵消的量：将峰值降低到目标门限
        clip_excess = peak_val - target_amp * exp(1j * angle(peak_val));

        % 生成圆周移位的核函数，使其峰值对准 n_peak
        shift = n_peak - n_p;
        p_shifted = circshift(p_time, shift);

        % 从当前信号中减去缩放后的核函数
        x_current = x_current - mu_step * clip_excess * p_shifted;
    end

    % 提取预留子载波上的峰值抵消信号
    X_full = fft(x_current, N_os) / sqrt(N_os);
    X_result = zeros(N, 1);
    X_result(1:N/2) = X_full(1:N/2);
    X_result(N/2+1:end) = X_full(N_os - N/2 + 1:end);
    C = X_result - X_data;  % 预留子载波上的抵消信号

    x_reduced = x_current;
end
