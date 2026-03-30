function X_hat = ofdm_receiver(y_os, params, H)
% OFDM系统组件6：从时域信号中恢复频域符号。
%   输入：
%     y_os   - 接收到的过采样时域信号
%     params - get_default_params() 返回的参数结构体
%     H      - 数据子载波上的信道频率响应（N_fft x 1）
%   输出：
%     X_hat  - 估计的频域符号（N_fft x 1）

    N = params.N_fft;
    L = params.L;
    N_os = N * L;

    Y_full = fft(y_os, N_os) / sqrt(N_os);
    X_hat = [Y_full(1:N/2); Y_full(N_os - N/2 + 1:end)];

    if nargin >= 3 && ~isempty(H)
        X_hat = X_hat ./ H;
    end
end
