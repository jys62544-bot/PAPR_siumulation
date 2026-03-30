function [x_reduced, C, reserved_idx] = tone_reservation(X, params, reserved_ratio, max_iter)
% 预留子载波法降低 PAPR，保留一部分子载波用于承载峰值抵消信号，采用基于迭代限幅的梯度投影算法进行优化。
%   输入：
%     X              - 频域数据（N_fft x 1）
%     params         - get_default_params() 返回的参数结构体
%     reserved_ratio - 预留子载波比例（如 0.10 表示 10%）
%     max_iter       - 最大迭代次数（如 10）
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

    X_data = X;
    X_data(reserved_idx) = 0;

    X_padded = [X_data(1:N/2); zeros((L-1)*N, 1); X_data(N/2+1:end)];
    x_current = sqrt(N_os) * ifft(X_padded, N_os);

    C = zeros(N, 1);

    % 预计算各预留子载波的时域基向量
    F_reserved = zeros(N_os, N_reserved);
    for i = 1:N_reserved
        k = reserved_idx(i);
        if k <= N/2
            k_os = k;
        else
            k_os = N_os - N + k;
        end
        e_k = zeros(N_os, 1);
        e_k(k_os) = 1;
        F_reserved(:, i) = sqrt(N_os) * ifft(e_k, N_os);
    end

    % 迭代投影降低峰值
    for iter = 1:max_iter
        [~, n_peak] = max(abs(x_current));

        target = 0.9 * max(abs(x_current));

        if abs(x_current(n_peak)) > target
            clip_val = x_current(n_peak) - target * exp(1j * angle(x_current(n_peak)));

            proj_coeffs = conj(F_reserved(n_peak, :)).' / N_os;
            c_update = F_reserved * proj_coeffs * clip_val;

            x_current = x_current - c_update;

            C_update = zeros(N, 1);
            C_full_update = fft(c_update, N_os) / sqrt(N_os);
            C_update(1:N/2) = C_full_update(1:N/2);
            C_update(N/2+1:end) = C_full_update(N_os - N/2 + 1:end);
            C = C + C_update;
        end
    end

    x_reduced = x_current;
end
