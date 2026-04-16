function info_bits = golay_decode(rx_symbols, params, perm_idx)
% GOLAY_DECODE  从接收 QPSK 符号中解码 Davis-Jedwab Golay 编码的信息比特
%
% 使用 Davis-Jedwab Algorithm 14 (h=4 情况)。
% Coset 索引（排列 π）假设通过 side information 完美已知。
%
% 输入：
%   rx_symbols - 均衡后的接收频域符号 (N x 1, 复数)
%   params     - 系统参数（需要 N_fft）
%   perm_idx   - 排列索引（side information，1-based）
%
% 输出：
%   info_bits  - 解码后的信息比特 (n_info x 1)

    N = params.N_fft;
    m = log2(N);
    h = 4;

    % 获取排列表（persistent 缓存）
    persistent perm_table_dec perm_table_dec_m
    if isempty(perm_table_dec) || perm_table_dec_m ~= m
        perm_table_dec = golay_build_perm_table(m);
        perm_table_dec_m = m;
    end
    K_coset = floor(log2(size(perm_table_dec, 1)));
    perm = perm_table_dec(perm_idx, :);

    % 量化接收符号到 Z_4
    % QPSK: j^c → c ∈ {0,1,2,3}
    phases = angle(rx_symbols);
    phases = mod(phases, 2*pi);
    rx_z4 = mod(round(phases / (pi/2)), h);

    % 构造二进制向量表
    bin_table = zeros(N, m);
    for bit = 1:m
        bin_table(:, bit) = double(bitand(uint32(0:N-1)', uint32(2^(bit-1))) > 0);
    end

    % Algorithm 14: 逐阶 FHT 解码 (h=4, 2 passes)
    r = rx_z4;
    d_coeffs = zeros(m + 1, 1);

    % Pass t=0: 处理 mod 2
    % 二阶项 mod 2 = 0（因为系数 2 ≡ 0 mod 2）
    % 所以 r mod 2 = Σ d_k*x_k + d_0 (mod 2) → RM(1,m) 码字
    r_mod2 = mod(r, 2);
    s = 1 - 2 * r_mod2;  % {0,1} → {+1,-1}

    H_transform = golay_fht(s);

    [~, max_pos] = max(abs(H_transform));
    max_pos = max_pos - 1;  % 0-based

    % 恢复 d_1,...,d_m 的 mod-2 部分
    d_mod2 = zeros(m, 1);
    for k = 1:m
        d_mod2(k) = double(bitand(uint32(max_pos), uint32(2^(k-1))) > 0);
    end
    d0_mod2 = double(H_transform(max_pos + 1) < 0);

    % Pass t=1: 恢复完整 Z_4 系数
    % 计算并减去二阶项
    quad_sum = zeros(N, 1);
    for k = 1:m-1
        quad_sum = quad_sum + bin_table(:, perm(k)) .* bin_table(:, perm(k+1));
    end
    r_stripped = mod(r - 2 * quad_sum, h);

    % 减去 mod-2 线性部分
    linear_mod2 = mod(bin_table * d_mod2 + d0_mod2, h);
    r2 = mod(r_stripped - linear_mod2, h);

    % r2 的值应为 {0, 2}（理想），除以 2 得 Z_2 序列
    r2_div2 = mod(round(double(r2) / 2), 2);

    % FHT 解码高位
    s2 = 1 - 2 * r2_div2;
    H_transform2 = golay_fht(s2);
    [~, max_pos2] = max(abs(H_transform2));
    max_pos2 = max_pos2 - 1;

    d_high = zeros(m, 1);
    for k = 1:m
        d_high(k) = double(bitand(uint32(max_pos2), uint32(2^(k-1))) > 0);
    end
    d0_high = double(H_transform2(max_pos2 + 1) < 0);

    % 组合：d_k = d_mod2(k) + 2*d_high(k) (mod 4)
    for k = 1:m
        d_coeffs(k + 1) = mod(d_mod2(k) + 2 * d_high(k), h);
    end
    d_coeffs(1) = mod(d0_mod2 + 2 * d0_high, h);

    % 组装信息比特
    n_info = K_coset + 2 * (m + 1);
    info_bits = zeros(n_info, 1);

    % Coset 索引比特（side info，直接还原）
    info_bits(1:K_coset) = golay_int2bit(perm_idx - 1, K_coset);

    % d_0, d_1, ..., d_m 的比特（每个 2 bits, MSB first）
    ptr = K_coset + 1;
    for k = 0:m
        val = d_coeffs(k + 1);
        info_bits(ptr)     = floor(val / 2);
        info_bits(ptr + 1) = mod(val, 2);
        ptr = ptr + 2;
    end
end


function H = golay_fht(x)
% GOLAY_FHT  Fast Hadamard Transform（原地蝶形算法）
    x = x(:);
    N = length(x);
    H = x;
    step = 1;
    while step < N
        for i = 0:2*step:N-1
            for j = 0:step-1
                a = H(i + j + 1);
                b = H(i + j + step + 1);
                H(i + j + 1)        = a + b;
                H(i + j + step + 1) = a - b;
            end
        end
        step = step * 2;
    end
end


function bits = golay_int2bit(val, n_bits)
% GOLAY_INT2BIT  整数转二进制比特向量
    bits = zeros(n_bits, 1);
    for k = 1:n_bits
        bits(k) = double(bitand(uint32(val), uint32(2^(n_bits - k))) > 0);
    end
end
