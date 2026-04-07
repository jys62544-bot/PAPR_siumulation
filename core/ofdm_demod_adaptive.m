function bits = ofdm_demod_adaptive(symbols, mod_map)
% OFDM_DEMOD_ADAPTIVE 逐子载波自适应硬判决解调器。
% 根据 mod_map 对每个子载波做对应星座的硬判决解映射。
%   输入：
%     symbols - N x 1 接收到的复数符号
%     mod_map - N x 1 向量，每个子载波的 bps（1/2/4/6/8）
%   输出：
%     bits    - 解调后的比特列向量

    N = length(mod_map);
    total_bits = sum(mod_map);
    bits = zeros(total_bits, 1);

    bit_ptr = 0;
    for k = 1:N
        bps = mod_map(k);
        b = demap_symbol(symbols(k), bps);
        bits(bit_ptr + 1 : bit_ptr + bps) = b;
        bit_ptr = bit_ptr + bps;
    end
end

function b = demap_symbol(s, bps)
% 将一个星座点硬判决解映射为 bps 个比特
    switch bps
        case 1  % BPSK
            b = real(s) < 0;

        case 2  % QPSK
            I = real(s) * sqrt(2);
            Q = imag(s) * sqrt(2);
            b = [I < 0; Q < 0];

        case 4  % 16QAM
            I = real(s) * sqrt(10);
            Q = imag(s) * sqrt(10);
            b = [pam4_to_gray(I); pam4_to_gray(Q)];

        case 6  % 64QAM
            I = real(s) * sqrt(42);
            Q = imag(s) * sqrt(42);
            b = [pam8_to_gray(I); pam8_to_gray(Q)];

        case 8  % 256QAM
            I = real(s) * sqrt(170);
            Q = imag(s) * sqrt(170);
            b = [pam16_to_gray(I); pam16_to_gray(Q)];
    end
end

function b = pam4_to_gray(val)
    if val < -2, level = -3;
    elseif val < 0, level = -1;
    elseif val < 2, level = 1;
    else, level = 3;
    end
    map = [-3 0 0; -1 0 1; 1 1 1; 3 1 0];
    b = map(map(:,1)==level, 2:3).';
end

function b = pam8_to_gray(val)
    levels = [-7, -5, -3, -1, 1, 3, 5, 7];
    [~, idx] = min(abs(val - levels));
    decided = levels(idx);
    gray_map = [
        -7, 0, 0, 0;
        -5, 0, 0, 1;
        -3, 0, 1, 1;
        -1, 0, 1, 0;
         1, 1, 1, 0;
         3, 1, 1, 1;
         5, 1, 0, 1;
         7, 1, 0, 0];
    b = gray_map(gray_map(:,1)==decided, 2:4).';
end

function b = pam16_to_gray(val)
    levels = -15:2:15;
    [~, idx] = min(abs(val - levels));
    decided = levels(idx);
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
    b = gray_map(gray_map(:,1)==decided, 2:5).';
end
