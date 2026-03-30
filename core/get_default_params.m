function params = get_default_params()
% OFDM系统组件3：OFDM 系统默认仿真参数获取函数。
% params = get_default_params() 返回包含所有仿真参数的结构体。在每个脚本中首先调用此函数，然后根据需要覆盖特定字段。

    params.N_fft     = 256;
    params.N_cp      = 64;        % 循环前缀长度默认为 N_fft/4
    params.mod_type  = '16QAM';
    params.L         = 4;         % 过采样倍数
    params.N_sym     = 10000;     % 仿真符号数
    params.snr_range = 0:2:30;    % 信噪比范围（dB）
    params.channel   = 'AWGN';

    % 派生参数
    params.N_os = params.N_fft * params.L;  % 过采样 IFFT 点数

    % 各调制方式下每符号比特数
    switch params.mod_type
        case 'QPSK',  params.bps = 2;
        case '16QAM', params.bps = 4;
        case '64QAM', params.bps = 6;
    end
end
