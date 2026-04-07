function [x_golay, X_golay, info] = golay_coding(X, params, n_candidates)
% GOLAY_CODING 基于 Golay 互补序列对的 OFDM PAPR 降低编码法（Rate-1/2）。
%
% 编码原理（Rate-1/2 方案）：
%   1. 将输入频域符号 X（N x 1）的前 N/2 个视为信息符号 d
%   2. 生成 Golay 互补序列对 (a, b)，各长 N/2，满足 |A(z)|^2 + |B(z)|^2 = 2·(N/2)
%   3. 编码：X_coded = [d .* a; d .* b]，利用互补性质约束 PAPR
%   4. 从 n_candidates 对 Golay 序列中选 PAPR 最低的候选
%
% 编码率 = 1/2（N/2 个信息符号 → N 个编码符号），数据率减半。
% 接收端：d_hat = X_hat(1:N/2) .* a（因为 a 是 ±1，a.*a=1）
%
% 与其他算法的公平比较：
%   - 使用同一输入 X 的前 N/2 个符号作为信息数据
%   - 数据率减半是方法固有代价（类似 TR 的频谱效率损失）
%   - BER 仅统计前 N/2 符号对应的比特
%
%   输入：
%     X            - 频域数据符号（N_fft x 1），取前 N/2 个作为信息符号
%     params       - get_default_params() 返回的参数结构体
%     n_candidates - 候选 Golay 对数（默认 16）
%   输出：
%     x_golay  - PAPR 降低后的过采样时域信号（N_fft*L x 1）
%     X_golay  - 频域编码后的符号（N_fft x 1）
%     info     - 结构体：
%                  .papr       - 最优候选的 PAPR (dB)
%                  .seq_idx    - 最优候选序列编号
%                  .a_seq      - 最优 Golay 序列 a（N/2 x 1，±1）
%                  .b_seq      - 最优 Golay 序列 b（N/2 x 1，±1）
%                  .data_idx   - 信息符号对应的子载波索引（1:N/2）

    if nargin < 3 || isempty(n_candidates)
        n_candidates = 16;
    end

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    N_half = N / 2;
    m_half = log2(N_half);  % Golay 对长度 = N/2 = 2^m_half

    % 取前 N/2 个符号作为信息数据
    d = X(1:N_half);

    best_papr = inf;
    x_golay = [];
    X_golay = [];
    best_idx = 1;
    best_a = [];
    best_b = [];

    for c = 1:n_candidates
        % 生成长度 N/2 的 Golay 互补序列对
        [a, b] = gen_golay_pair(m_half, c);
        a = a(:);
        b = b(:);

        % Rate-1/2 编码：前 N/2 子载波放 d.*a，后 N/2 子载波放 d.*b
        X_c = [d .* a; d .* b];

        % 过采样 IFFT
        X_padded = [X_c(1:N/2); zeros((L-1)*N, 1); X_c(N/2+1:end)];
        x_c = sqrt(N_os) * ifft(X_padded, N_os);

        p = 10 * log10(max(abs(x_c).^2) / mean(abs(x_c).^2));

        if p < best_papr
            best_papr = p;
            x_golay = x_c;
            X_golay = X_c;
            best_idx = c;
            best_a = a;
            best_b = b;
        end
    end

    info.papr     = best_papr;
    info.seq_idx  = best_idx;
    info.a_seq    = best_a;   % 接收端还原用
    info.b_seq    = best_b;
    info.data_idx = (1:N_half)';  % 信息符号的子载波索引
end

function [a, b] = gen_golay_pair(m, seed)
% GEN_GOLAY_PAIR 拼接式递归构造长度 2^m 的 Golay 互补序列对（±1）。
% 使用 [a, W*b] / [a, -W*b] 拼接形式，保证所有元素严格为 ±1。
    rng_state = rng;
    rng(seed);
    a = 1;
    b = 1;
    for k = 1:m
        W = 2 * randi([0 1]) - 1;  % +1 或 -1
        a_new = [a(:);  W * b(:)];
        b_new = [a(:); -W * b(:)];
        a = a_new;
        b = b_new;
    end
    a = a(:);
    b = b(:);
    rng(rng_state);
end
