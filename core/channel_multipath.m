function [y, H, h] = channel_multipath(x, params, model_type)
% OFDM系统组件2：多径衰落信道。
%   输入：
%     x          - 发送时域信号
%     params     - get_default_params() 返回的参数结构体
%     model_type - 信道模型，我查阅文献选择了常用的EPA或ETU
%
%   输出：
%     y - 经过多径后的接收信号
%     H - N_fft 个子载波上的信道频率响应
%     h - 时域信道冲激响应

    N = params.N_fft;
    L = params.L;

    switch model_type
        case 'EPA'
            delays_ns = [0, 30, 70, 90, 110, 190, 410];
            powers_dB = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];
        case 'ETU'
            delays_ns = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000];
            powers_dB = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0];
        otherwise
            error('Unknown channel model: %s. Use EPA or ETU.', model_type);
    end

    Ts = 1 / (15.36e6);
    delay_samples = round(delays_ns * 1e-9 / Ts);

    powers_lin = 10.^(powers_dB / 10);
    powers_lin = powers_lin / sum(powers_lin);
    h_coeffs = sqrt(powers_lin/2) .* (randn(size(powers_lin)) + 1j*randn(size(powers_lin)));

    h_len = max(delay_samples) + 1;
    h = zeros(h_len, 1);
    for i = 1:length(delay_samples)
        h(delay_samples(i) + 1) = h(delay_samples(i) + 1) + h_coeffs(i);
    end

    y = conv(x, h);
    y = y(1:length(x));

    H_full = fft(h, N*L);
    H = [H_full(1:N/2); H_full(N*L - N/2 + 1:end)];
end
