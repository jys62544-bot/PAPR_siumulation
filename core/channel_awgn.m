function [y, noise_var] = channel_awgn(x, snr_dB)
% 加性高斯白噪声信道
%   输入：
%     x      - 发送信号，格式为复数列向量
%     snr_dB - 信噪比（dB）
%   输出：
%     y         - 接收信号：y = x + 噪声
%     noise_var - 噪声方差

    sig_power = mean(abs(x).^2);
    noise_var = sig_power / (10^(snr_dB/10));
    noise = sqrt(noise_var/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = x + noise;
end
