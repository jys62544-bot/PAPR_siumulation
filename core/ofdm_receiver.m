function X_hat = ofdm_receiver(y_os, params, H)
% 从时域信号中恢复频域符号。
% 自动检测输入是否含循环前缀，若含则先去除 CP 再做 FFT。
%   输入：
%     y_os   - 接收到的过采样时域信号（可含或不含 CP）
%     params - get_default_params() 返回的参数结构体
%     H      - 数据子载波上的信道频率响应（N_fft x 1）
%   输出：
%     X_hat  - 估计的频域符号（N_fft x 1）

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    N_cp_os = params.N_cp * L;

    % 自动检测并去除循环前缀
    if length(y_os) >= N_os + N_cp_os
        y_os = y_os(N_cp_os + 1 : N_cp_os + N_os);
    end

    Y_full = fft(y_os, N_os) / sqrt(N_os);
    X_hat = [Y_full(1:N/2); Y_full(N_os - N/2 + 1:end)];

    if nargin >= 3 && ~isempty(H)
        X_hat = X_hat ./ H;
    end
end
