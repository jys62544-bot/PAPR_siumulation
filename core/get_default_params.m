function params = get_default_params()
% OFDM 系统默认仿真参数获取函数。
% params = get_default_params() 返回包含所有仿真参数的结构体。
% 在每个脚本中首先调用此函数，然后根据需要覆盖特定字段。

    params.N_fft     = 256;
    params.N_cp      = 64;        % 循环前缀长度（N_fft/4）
    params.mod_type  = '16QAM';   % 默认调制方式（用于非自适应子载波的回退）
    params.L         = 4;         % 过采样倍数
    params.N_sym     = 10000;     % 仿真符号数
    params.snr_range = 0:2:30;    % 信噪比范围（dB）

    % 自适应调制参数
    % 使用多径信道 + 逐子载波自适应调制
    params.channel_model  = 'ETU';         % 多径信道模型（EPA 或 ETU）
    params.snr_thresholds = [6, 12, 18, 24]; % 调制切换 SNR 门限（dB）
    % 门限含义：<6→BPSK, 6~12→QPSK, 12~18→16QAM, 18~24→64QAM, >=24→256QAM

    % 派生参数
    params.N_os = params.N_fft * params.L;  % 过采样 IFFT 点数

    % 各调制方式下每符号比特数（内部函数回退用）
    switch params.mod_type
        case 'BPSK',  params.bps = 1;
        case 'QPSK',  params.bps = 2;
        case '16QAM', params.bps = 4;
        case '64QAM', params.bps = 6;
        case '256QAM', params.bps = 8;
    end
end
