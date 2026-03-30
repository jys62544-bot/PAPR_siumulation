function [x_best, b_best, papr_best] = pts(X, params, V, partition_type, W)
% 部分传输序列法降低 PAPR，将子载波划分为 V 个不相交子块，穷举所有相位旋转因子组合，选择 PAPR 最低的组合发送。
%   输入：
%     X              - 频域数据
%     params         - get_default_params() 返回的参数结构体
%     V              - 子块数量（2、4、8）
%     partition_type - 分块方式：'adjacent'、'interleaved'或 'random'
%     W              - 相位因子集合
%   输出：
%     x_best    - PAPR 最低的时域信号
%     b_best    - 最优相位因子向量
%     papr_best - 达到的 PAPR（dB）

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    block_size = N / V;

    switch partition_type
        case 'adjacent'
            partition = reshape(1:N, block_size, V);

        case 'interleaved'
            partition = zeros(block_size, V);
            for v = 1:V
                partition(:, v) = v:V:N;
            end

        case 'random'
            perm = randperm(N);
            partition = reshape(perm, block_size, V);

        otherwise
            error('Unknown partition type: %s', partition_type);
    end

    x_sub = zeros(N_os, V);
    for v = 1:V
        X_v = zeros(N, 1);
        X_v(partition(:, v)) = X(partition(:, v));

        X_v_padded = [X_v(1:N/2); zeros((L-1)*N, 1); X_v(N/2+1:end)];
        x_sub(:, v) = sqrt(N_os) * ifft(X_v_padded, N_os);
    end

    W_size = length(W);
    n_comb = W_size^(V-1);  % 第一个子块相位固定为 1

    papr_best = inf;
    b_best = ones(V, 1);
    x_best = [];

    for idx = 0:(n_comb - 1)
        b = ones(V, 1);
        temp = idx;
        for v = V:-1:2
            b(v) = W(mod(temp, W_size) + 1);
            temp = floor(temp / W_size);
        end

        x_comb = x_sub * b;
        papr_val = compute_papr(x_comb);

        if papr_val < papr_best
            papr_best = papr_val;
            b_best = b;
            x_best = x_comb;
        end
    end
end
