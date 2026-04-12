function [y, noise_var] = channel_awgn_fixed(x, snr_dB, ref_power)
% CHANNEL_AWGN_FIXED 使用固定参考功率的 AWGN 信道，噪声功率由 ref_power 和 snr_dB 决定，而非由 x 的实际功率决定，保证不同 PAPR 降低算法在相同 SNR 下接收到相同量级的噪声。
%   输入：
%     x         - 发送信号
%     snr_dB    - 信噪比（dB），基于 ref_power 定义
%     ref_power - 参考信号功率（标量），通常为原始信号的平均功率
%   输出：
%     y         - 接收信号
%     noise_var - 噪声方差

    noise_var = ref_power / (10^(snr_dB / 10));
    noise = sqrt(noise_var / 2) * (randn(size(x)) + 1j * randn(size(x)));
    y = x + noise;
end
