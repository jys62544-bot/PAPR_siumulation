function [symbols, bits, total_bits] = ofdm_mod_adaptive(mod_map)
% OFD_MOD_ADAPTIVE 逐子载波自适应调制器。
% 根据 mod_map 为每个子载波随机生成比特并映射到对应星座点。
% 支持 BPSK(1), QPSK(2), 16QAM(4), 64QAM(6), 256QAM(8)。
%   输入：
%     mod_map - N x 1 向量，每个子载波的 bps（1/2/4/6/8）
%   输出：
%     symbols    - N x 1 复数星座点
%     bits       - 总比特列向量（按子载波顺序拼接）
%     total_bits - 该 OFDM 符号的总比特数

    N = length(mod_map);
    total_bits = sum(mod_map);
    bits = randi([0 1], total_bits, 1);
    symbols = zeros(N, 1);

    bit_ptr = 0;
    for k = 1:N
        bps = mod_map(k);
        b = bits(bit_ptr + 1 : bit_ptr + bps);
        symbols(k) = map_symbol(b, bps);
        bit_ptr = bit_ptr + bps;
    end
end

function s = map_symbol(b, bps)
% 将 bps 个比特映射到一个星座点（格雷码映射）
    switch bps
        case 1  % BPSK
            s = 1 - 2*b(1);

        case 2  % QPSK
            I = 1 - 2*b(1);
            Q = 1 - 2*b(2);
            s = (I + 1j*Q) / sqrt(2);

        case 4  % 16QAM
            I = gray_to_pam4(b(1:2));
            Q = gray_to_pam4(b(3:4));
            s = (I + 1j*Q) / sqrt(10);

        case 6  % 64QAM
            I = gray_to_pam8(b(1:3));
            Q = gray_to_pam8(b(4:6));
            s = (I + 1j*Q) / sqrt(42);

        case 8  % 256QAM
            I = gray_to_pam16(b(1:4));
            Q = gray_to_pam16(b(5:8));
            s = (I + 1j*Q) / sqrt(170);
    end
end

function val = gray_to_pam4(bits2)
    gray_idx = bits2(1)*2 + bits2(2);
    map = [-3, -1, 3, 1];
    val = map(gray_idx + 1);
end

function val = gray_to_pam8(bits3)
    gray_idx = bits3(1)*4 + bits3(2)*2 + bits3(3);
    map = [-7, -5, -1, -3, 7, 5, 1, 3];
    val = map(gray_idx + 1);
end

function val = gray_to_pam16(bits4)
    gray_idx = bits4(1)*8 + bits4(2)*4 + bits4(3)*2 + bits4(4);
    map = [-15, -13, -9, -11, -1, -3, -7, -5, 15, 13, 9, 11, 1, 3, 7, 5];
    val = map(gray_idx + 1);
end
