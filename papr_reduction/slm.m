function [x_best, u_best, papr_best, P_all] = slm(X, params, U, P_all)
% SLM 选择映射法，用于降低 PAPR，通过将数据与 U 组不同随机相位序列相乘，生成 U 个候选 OFDM 符号，选择 PAPR 最小的一组发送。此方法不引入信号失真。
%   输入：
%     X     - 频域数据
%     params - get_default_params() 返回的参数结构体
%     U     - 候选相位序列数量（2、4、8、16）
%     P_all - 预生成的相位序列矩阵
%   输出：
%     x_best    - PAPR 最低的时域信号
%     u_best    - 所选相位序列索引
%     papr_best - 达到的 PAPR（dB）
%     P_all     - 使用的相位序列矩阵

    N = params.N_fft;
    L = params.L;
    N_os = N * L;

    if nargin < 4 || isempty(P_all)
        P_all = exp(1j * 2 * pi * rand(N, U));
        P_all(:, 1) = ones(N, 1);
    end

    papr_best = inf;
    u_best = 1;
    x_best = [];

    for u = 1:U
        X_u = X .* P_all(:, u);

        X_padded = [X_u(1:N/2); zeros((L-1)*N, 1); X_u(N/2+1:end)];
        x_u = sqrt(N_os) * ifft(X_padded, N_os);

        papr_u = compute_papr(x_u);

        if papr_u < papr_best
            papr_best = papr_u;
            u_best = u;
            x_best = x_u;
        end
    end
end
