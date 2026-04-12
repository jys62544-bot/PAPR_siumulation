function [x_os, X_freq, x_cp] = ofdm_transmitter(X, params)
% 生成过采样的时域 OFDM 信号。
% 通过在频域中心补零实现 L 倍过采样 IFFT，精确捕捉连续时间 PAPR 峰值。
%   输入：
%     X      - 频域符号向量
%     params - get_default_params() 返回的参数结构体
%   输出：
%     x_os   - 过采样时域信号，不含循环前缀
%     X_freq - 补零后的频域向量
%     x_cp   - 含循环前缀的过采样时域信号用于信道传输
    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    X_freq = [X(1:N/2); zeros((L-1)*N, 1); X(N/2+1:end)];
    x_os = sqrt(N_os) * ifft(X_freq, N_os);

    % 添加循环前缀
    N_cp_os = params.N_cp * L;
    x_cp = [x_os(end - N_cp_os + 1 : end); x_os];
end
