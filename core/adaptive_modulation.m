function [mod_map, bps_vec] = adaptive_modulation(H, snr_avg_db, thresholds)
% ADAPTIVE_MODULATION 根据信道频率响应自适应分配每个子载波的调制方式。
%   输入：
%     H            - N x 1 信道频率响应（复数）
%     snr_avg_db   - 系统平均 SNR（dB）
%     thresholds   - 调制切换 SNR 门限向量（默认 [6, 12, 18, 24]）
%                    对应 BPSK→QPSK→16QAM→64QAM→256QAM
%   输出：
%     mod_map - N x 1 数值向量，每个子载波的每符号比特数（1/2/4/6/8）
%     bps_vec - 同 mod_map

    if nargin < 3 || isempty(thresholds)
        thresholds = [6, 12, 18, 24];
    end

    N = length(H);

    % 计算每个子载波的有效 SNR
    % SNR_k = SNR_avg + 20*log10(|H(k)|)
    H_mag = abs(H);
    H_mag(H_mag < 1e-10) = 1e-10;  % 避免 log(0)
    snr_per_sc = snr_avg_db + 20 * log10(H_mag);

    % 按门限分配调制阶数
    mod_map = ones(N, 1);  % 默认 BPSK (1 bps)
    mod_map(snr_per_sc >= thresholds(1)) = 2;   % QPSK
    mod_map(snr_per_sc >= thresholds(2)) = 4;   % 16QAM
    mod_map(snr_per_sc >= thresholds(3)) = 6;   % 64QAM
    mod_map(snr_per_sc >= thresholds(4)) = 8;   % 256QAM

    bps_vec = mod_map;
end
