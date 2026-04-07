function bits = ofdm_demod(symbols, mod_type)
% OFDM系统组件4：OFDM系统硬判决星座解映射。
%   输入：
%     symbols  - 接收到的复数符号列向量
%     mod_type - 调制类型：'BPSK'、'QPSK'、'16QAM'、'64QAM' 或 '256QAM'
%   输出：
%     bits     - 检测到的二进制数据列向量

    switch mod_type
        case 'BPSK'
            bits = double(real(symbols(:)) < 0);

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

        case '256QAM'
            I = real(symbols) * sqrt(170);
            Q = imag(symbols) * sqrt(170);
            I_bits = pam16_to_gray(I);
            Q_bits = pam16_to_gray(Q);
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

function bits4 = pam16_to_gray(val)
    levels = -15:2:15;
    [~, idx] = min(abs(val(:) - levels), [], 2);
    decided = levels(idx).';
    gray_map = [
        -15, 0,0,0,0;
        -13, 0,0,0,1;
         -9, 0,0,1,0;
        -11, 0,0,1,1;
         -1, 0,1,0,0;
         -3, 0,1,0,1;
         -7, 0,1,1,0;
         -5, 0,1,1,1;
         15, 1,0,0,0;
         13, 1,0,0,1;
          9, 1,0,1,0;
         11, 1,0,1,1;
          1, 1,1,0,0;
          3, 1,1,0,1;
          7, 1,1,1,0;
          5, 1,1,1,1];
    bits4 = zeros(length(val), 4);
    for i = 1:size(gray_map, 1)
        mask = (decided == gray_map(i, 1));
        bits4(mask, :) = repmat(gray_map(i, 2:5), sum(mask), 1);
    end
end
