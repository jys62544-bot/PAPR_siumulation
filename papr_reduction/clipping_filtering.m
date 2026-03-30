function [x_clipped, X_clipped] = clipping_filtering(x, params, CR, n_iter)
% 限幅滤波法降低 PAPR，将时域信号幅度限制在 CR * rms(x)，然后进行频域滤波以消除带外频谱再生。可选择迭代 n_iter 次限幅+滤波。
%   输入：
%     x      - 过采样时域 OFDM 信号
%     params - get_default_params() 返回的参数结构体
%     CR     - 限幅比
%     n_iter - 限幅+滤波迭代次数，默认为 1
%   输出：
%     x_clipped - 限幅并滤波后的时域信号
%     X_clipped - N_fft 数据子载波上的频域表示

    if nargin < 4
        n_iter = 1;
    end

    N = params.N_fft;
    L = params.L;
    N_os = N * L;

    x_clipped = x;

    for iter = 1:n_iter
        A = CR * rms(x_clipped);
        mag = abs(x_clipped);
        clipped_mag = min(mag, A);
        x_clipped = clipped_mag .* exp(1j * angle(x_clipped));

        X_full = fft(x_clipped, N_os);
        X_filtered = zeros(N_os, 1);
        X_filtered(1:N/2) = X_full(1:N/2);
        X_filtered(N_os - N/2 + 1:end) = X_full(N_os - N/2 + 1:end);
        x_clipped = ifft(X_filtered, N_os);
    end

    X_full_final = fft(x_clipped, N_os);
    X_clipped = [X_full_final(1:N/2); X_full_final(N_os - N/2 + 1:end)];
end
