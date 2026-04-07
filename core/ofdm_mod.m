function symbols = ofdm_mod(bits, mod_type)
% OFDM系统组件5：OFDM_MOD 的格雷码星座映射。
%   输入：
%     bits     - 二进制数据列向量（0/1）
%     mod_type - 调制类型：'BPSK'、'QPSK'、'16QAM'、'64QAM' 或 '256QAM'
%   输出：
%     symbols  - 复数星座点列向量

    switch mod_type
        case 'BPSK'
            bps = 1;
            N_sym = floor(length(bits) / bps);
            bits = bits(1:N_sym*bps);
            symbols = 1 - 2*bits(:);

        case 'QPSK'
            bps = 2;
            N_sym = floor(length(bits) / bps);
            bits = bits(1:N_sym*bps);
            bit_groups = reshape(bits, bps, N_sym).';
            I = 1 - 2*bit_groups(:,1);
            Q = 1 - 2*bit_groups(:,2);
            symbols = (I + 1j*Q) / sqrt(2);

        case '16QAM'
            bps = 4;
            N_sym = floor(length(bits) / bps);
            bits = bits(1:N_sym*bps);
            bit_groups = reshape(bits, bps, N_sym).';
            I_bits = bit_groups(:, 1:2);
            Q_bits = bit_groups(:, 3:4);
            I = gray_to_pam4(I_bits);
            Q = gray_to_pam4(Q_bits);
            symbols = (I + 1j*Q) / sqrt(10);

        case '64QAM'
            bps = 6;
            N_sym = floor(length(bits) / bps);
            bits = bits(1:N_sym*bps);
            bit_groups = reshape(bits, bps, N_sym).';
            I_bits = bit_groups(:, 1:3);
            Q_bits = bit_groups(:, 4:6);
            I = gray_to_pam8(I_bits);
            Q = gray_to_pam8(Q_bits);
            symbols = (I + 1j*Q) / sqrt(42);

        case '256QAM'
            bps = 8;
            N_sym = floor(length(bits) / bps);
            bits = bits(1:N_sym*bps);
            bit_groups = reshape(bits, bps, N_sym).';
            I_bits = bit_groups(:, 1:4);
            Q_bits = bit_groups(:, 5:8);
            I = gray_to_pam16(I_bits);
            Q = gray_to_pam16(Q_bits);
            symbols = (I + 1j*Q) / sqrt(170);

        otherwise
            error('Unsupported modulation: %s', mod_type);
    end
end

function val = gray_to_pam4(bits2)
    gray_idx = bits2(:,1)*2 + bits2(:,2);
    map = [-3, -1, 3, 1];
    val = map(gray_idx + 1).';
end

function val = gray_to_pam8(bits3)
    gray_idx = bits3(:,1)*4 + bits3(:,2)*2 + bits3(:,3);
    map = [-7, -5, -1, -3, 7, 5, 1, 3];
    val = map(gray_idx + 1).';
end

function val = gray_to_pam16(bits4)
    gray_idx = bits4(:,1)*8 + bits4(:,2)*4 + bits4(:,3)*2 + bits4(:,4);
    map = [-15, -13, -9, -11, -1, -3, -7, -5, 15, 13, 9, 11, 1, 3, 7, 5];
    val = map(gray_idx + 1).';
end
