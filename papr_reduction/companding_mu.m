function [x_comp, expand_func] = companding_mu(x, mu)
% μ律压扩降低 PAPR，对 OFDM 信号幅度施加 μ 律压缩，同时保留相位信息，并返回用于接收端解压扩的函数句柄。
%   输入：
%     x  - 时域 OFDM 信号
%     mu - 压扩参数
%   输出：
%     x_comp      - 压扩后的时域信号
%     expand_func - 解压扩函数句柄：x_recovered = expand_func(y)

    v_max = max(abs(x));

    mag = abs(x);
    phase = angle(x);
    comp_mag = v_max * log(1 + mu * mag / v_max) / log(1 + mu);
    x_comp = comp_mag .* exp(1j * phase);

    % 解压扩使用接收信号自身的最大幅度，避免信道增益失配
    expand_func = @(y) de_compand_mu(y, mu);
end

function x_rec = de_compand_mu(y, mu)
% 接收端 μ 律解压扩（自适应 v_max）
    v_max_rx = max(abs(y)) + eps;  % eps 防止全零输入时除零
    mag_y = abs(y);
    phase_y = angle(y);
    rec_mag = v_max_rx * ((1 + mu).^(mag_y / v_max_rx) - 1) / mu;
    x_rec = rec_mag .* exp(1j * phase_y);
end
