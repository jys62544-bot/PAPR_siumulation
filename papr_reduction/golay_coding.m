function [x_golay, X_golay, info] = golay_coding(~, params, ~)
% GOLAY_CODING  Davis-Jedwab Golay 序列编码，PMEPR 严格 ≤ 3 dB
% 基于 Davis & Jedwab (1999) "Peak-to-Mean Power Control in OFDM, Golay Complementary Sequences, and Reed-Muller Codes" 的 Corollary 9。
% 每个 OFDM 码字是 Z_4 上长度 N=2^m 的 Golay 序列，映射为 QPSK 符号。信息比特编码到 RM_4(1,m) 的 coset 中，PMEPR 严格 ≤ 2 (即 ≤ 3.01 dB)。
% 输入：
%   params   - get_default_params() 返回的参数结构体（需要 N_fft, L）
%
% 输出：
%   x_golay  - 过采样时域 QPSK Golay 序列 (N_fft*L x 1)
%   X_golay  - 频域 QPSK 符号 (N_fft x 1)
%   info     - 结构体：
%                .papr       - PAPR (dB)
%                .info_bits  - 发送的信息比特 (n_info x 1)
%                .n_info     - 信息比特数
%                .perm_idx   - 排列索引（side information）
%                .perm       - 所用排列向量
%                .d_coeffs   - Z_4 线性系数 [d_0; d_1; ...; d_m]
%                .codeword   - Z_4 码字 (N x 1)
    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    m = log2(N);
    h = 4;  % QPSK

    % 排列表
    persistent perm_table perm_table_m
    if isempty(perm_table) || perm_table_m ~= m
        perm_table = golay_build_perm_table(m);
        perm_table_m = m;
    end
    n_cosets = size(perm_table, 1);
    K_coset = floor(log2(n_cosets));  % coset 选择比特数

    % 随机生成信息比特
    n_info = K_coset + 2 * (m + 1);
    info_bits = randi([0 1], n_info, 1);

    % 编码
    [codeword, perm_idx, perm, d_coeffs] = golay_encode(info_bits, m, h, ...
                                                         K_coset, perm_table);

    % Z_4 → QPSK 符号：exp(j * 2π * c / 4) = j^c
    X_golay = exp(1j * pi/2 * codeword);

    % 过采样 IFFT
    X_padded = [X_golay(1:N/2); zeros((L-1)*N, 1); X_golay(N/2+1:end)];
    x_golay = sqrt(N_os) * ifft(X_padded, N_os);

    % 计算 PAPR
    papr_val = 10 * log10(max(abs(x_golay).^2) / mean(abs(x_golay).^2));

    % 输出信息
    info.papr      = papr_val;
    info.info_bits  = info_bits;
    info.n_info     = n_info;
    info.perm_idx   = perm_idx;
    info.perm       = perm;
    info.d_coeffs   = d_coeffs;
    info.codeword   = codeword;
end


function [codeword, perm_idx, perm, d_coeffs] = golay_encode(info_bits, m, h, K_coset, perm_table)
% GOLAY_ENCODE  将信息比特编码为 Z_h 上的 Golay 序列码字
%
% 码字构造（Corollary 9, h=4）：
%   c_j = 2 * Σ_{k=1}^{m-1} x_{π(k)} * x_{π(k+1)}
%         + Σ_{k=1}^{m} d_k * x_k + d_0   (mod h)

    N = 2^m;
    ptr = 1;

    % 1) Coset 索引 → 排列 π
    coset_bits = info_bits(ptr : ptr + K_coset - 1);
    perm_idx = bit2int(coset_bits) + 1;
    perm_idx = min(perm_idx, size(perm_table, 1));
    perm = perm_table(perm_idx, :);
    ptr = ptr + K_coset;

    % 2) 线性系数 d_0, d_1, ..., d_m ∈ Z_h
    d_coeffs = zeros(m + 1, 1);
    for k = 0:m
        bits_2 = info_bits(ptr : ptr + 1);
        d_coeffs(k + 1) = bits_2(1) * 2 + bits_2(2);
        ptr = ptr + 2;
    end

    % 3) 二进制向量表
    bin_table = zeros(N, m);
    for bit = 1:m
        bin_table(:, bit) = double(bitand(uint32(0:N-1)', uint32(2^(bit-1))) > 0);
    end

    % 4) 码字计算
    quad_sum = zeros(N, 1);
    for k = 1:m-1
        quad_sum = quad_sum + bin_table(:, perm(k)) .* bin_table(:, perm(k+1));
    end
    quad_term = mod(2 * quad_sum, h);
    linear_term = mod(bin_table * d_coeffs(2:end), h);
    codeword = mod(quad_term + linear_term + d_coeffs(1), h);
end


function bits = int2bit(val, n_bits)
    bits = zeros(n_bits, 1);
    for k = 1:n_bits
        bits(k) = double(bitand(uint32(val), uint32(2^(n_bits - k))) > 0);
    end
end


function val = bit2int(bits)
    n = length(bits);
    val = 0;
    for k = 1:n
        val = val + double(bits(k)) * 2^(n - k);
    end
end
