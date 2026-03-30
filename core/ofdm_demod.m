function bits = ofdm_demod(symbols, mod_type)
% OFDM系统组件4：OFDM系统硬判决星座解映射。
%   输入：
%     symbols  - 接收到的复数符号列向量
%     mod_type - 调制类型：'QPSK'、'16QAM' 或 '64QAM'
%   输出：
%     bits     - 检测到的二进制数据列向量

    switch mod_type
        case 'QPSK'
            I = real(symbols) * sqrt(2);
            Q = imag(symbols) * sqrt(2);
            bits = zeros(length(symbols)*2, 1);
            bits(1:2:end) = I < 0;
            bits(2:2:end) = Q < 0;

        case '16QAM'
            I = real(symbols) * sqrt(10);
            Q = imag(symbols) * sqrt(10);
            I_bits = pam4_to_gray(I);
            Q_bits = pam4_to_gray(Q);
            bits_mat = [I_bits, Q_bits];
            bits = reshape(bits_mat.', [], 1);

        case '64QAM'
            I = real(symbols) * sqrt(42);
            Q = imag(symbols) * sqrt(42);
            I_bits = pam8_to_gray(I);
            Q_bits = pam8_to_gray(Q);
            bits_mat = [I_bits, Q_bits];
            bits = reshape(bits_mat.', [], 1);

        otherwise
            error('Unsupported modulation: %s', mod_type);
    end
end

function bits2 = pam4_to_gray(val)
    level = zeros(size(val));
    level(val < -2) = -3;
    level(val >= -2 & val < 0) = -1;
    level(val >= 0 & val < 2) = 1;
    level(val >= 2) = 3;
    bits2 = zeros(length(val), 2);
    bits2(level == -3, :) = repmat([0 0], sum(level == -3), 1);
    bits2(level == -1, :) = repmat([0 1], sum(level == -1), 1);
    bits2(level ==  1, :) = repmat([1 1], sum(level ==  1), 1);
    bits2(level ==  3, :) = repmat([1 0], sum(level ==  3), 1);
end

function bits3 = pam8_to_gray(val)
% 对最近的 8-PAM 电平做硬判决并进行格雷解码
    % 判决边界：-6, -4, -2, 0, 2, 4, 6
    levels = [-7, -5, -3, -1, 1, 3, 5, 7];
    [~, idx] = min(abs(val(:) - levels), [], 2);
    decided = levels(idx).';
    bits3 = zeros(length(val), 3);
    bits3(decided == -7, :) = repmat([0 0 0], sum(decided == -7), 1);
    bits3(decided == -5, :) = repmat([0 0 1], sum(decided == -5), 1);
    bits3(decided == -3, :) = repmat([0 1 1], sum(decided == -3), 1);
    bits3(decided == -1, :) = repmat([0 1 0], sum(decided == -1), 1);
    bits3(decided ==  1, :) = repmat([1 1 0], sum(decided ==  1), 1);
    bits3(decided ==  3, :) = repmat([1 1 1], sum(decided ==  3), 1);
    bits3(decided ==  5, :) = repmat([1 0 1], sum(decided ==  5), 1);
    bits3(decided ==  7, :) = repmat([1 0 0], sum(decided ==  7), 1);
end
