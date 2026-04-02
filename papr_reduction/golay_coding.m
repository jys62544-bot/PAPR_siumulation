function [x_golay, X_golay, info] = golay_coding(N, params, n_candidates)
% GOLAY_CODING 基于格雷互补序列的 OFDM PAPR 降低编码法。利用 Golay 互补序列对的自相关互补特性，使 OFDM 信号的 PAPR受到理论上界约束
%   输入：
%     N            - 序列长度
%     params       - get_default_params() 返回的参数结构体
%     n_candidates - 候选序列数
%                    默认为 16
%   输出：
%     x_golay - 过采样时域信号
%     X_golay - 频域 Golay 编码序列
    if nargin < 3 || isempty(n_candidates)
        n_candidates = 16;
    end

    L = params.L;
    N_os = N * L;
    m = log2(N);

    best_papr = inf;
    x_golay = [];
    X_golay = [];
    best_idx = 1;

    for c = 1:n_candidates
        % 生成一对 Golay 互补序列
        [a, ~] = gen_golay_pair(m, c);
        % 将 ±1 序列映射为 QPSK 星座点
        phase_bits = randi([0 3], N, 1);
        qpsk_map = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);  % 格雷映射
        X_c = a(:) .* qpsk_map(phase_bits + 1).';
        % 过采样 IFFT
        X_padded = [X_c(1:N/2); zeros((L-1)*N, 1); X_c(N/2+1:end)];
        x_c = sqrt(N_os) * ifft(X_padded, N_os);

        p = 10 * log10(max(abs(x_c).^2) / mean(abs(x_c).^2));

        if p < best_papr
            best_papr = p;
            x_golay = x_c;
            X_golay = X_c;
            best_idx = c;
        end
    end

    info.papr = best_papr;
    info.seq_idx = best_idx;
end

function [a, b] = gen_golay_pair(m, seed)
% GEN_GOLAY_PAIR 递归构造长度 2^m 的 Golay 互补序列对。
    rng_state = rng; 
    rng(seed);
    % 基础对：长度为 1
    a = 1;
    b = 1;
    for k = 1:m
        % 随机延迟
        D = 2^(k-1);
        % 随机符号 W_k ∈ {+1, -1}
        W = 2 * randi([0 1]) - 1;
        % Golay 递归：
        a_delayed = [zeros(D, 1); b(:)];
        a_padded  = [a(:); zeros(D, 1)];

        a_new = a_padded + W * a_delayed;
        b_new = a_padded - W * a_delayed;

        a = a_new;
        b = b_new;
    end
    % 归一化为 ±1（Golay 对的元素幅度相同）
    a = sign(a);
    b = sign(b);
    rng(rng_state);
end
