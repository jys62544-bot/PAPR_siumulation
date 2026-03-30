function [x_comp, expand_func] = companding_a(x, A)
% A律压扩法降低 PAPR，对 OFDM 信号幅度施加 A 律压缩变换。
%   输入：
%     x - 时域 OFDM 信号
%     A - 压扩参数
%   输出：
%     x_comp      - 压扩后的信号
%     expand_func - 解压扩函数句柄

    v_max = max(abs(x));
    mag = abs(x) / v_max;
    phase = angle(x);

    denom = 1 + log(A);
    comp_mag = zeros(size(mag));

    low = mag < (1/A);
    comp_mag(low)  = A * mag(low) / denom;
    comp_mag(~low) = (1 + log(A * mag(~low))) / denom;

    x_comp = (comp_mag * v_max) .* exp(1j * phase);

    expand_func = @(y) de_compand_a(y, A, v_max);
end

function x_rec = de_compand_a(y, A, v_max)
% 接收端 A 律解压扩
    mag_y = abs(y) / v_max;
    phase_y = angle(y);
    denom = 1 + log(A);

    rec_mag = zeros(size(mag_y));
    threshold = 1 / denom;

    low = mag_y < threshold;
    rec_mag(low)  = mag_y(low) * denom / A;
    rec_mag(~low) = exp(mag_y(~low) * denom - 1) / A;

    x_rec = (rec_mag * v_max) .* exp(1j * phase_y);
end
