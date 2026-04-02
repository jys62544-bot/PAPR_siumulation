function x_cp = cp_add(x_os, params)
% CP_ADD 为过采样时域信号添加循环前缀。
% 供 PAPR 降低算法处理后的信号在送入信道前使用。
%   输入：
%     x_os   - 过采样时域信号（不含 CP）
%     params - get_default_params() 返回的参数结构体
%   输出：
%     x_cp   - 含循环前缀的时域信号

    N_cp_os = params.N_cp * params.L;
    x_cp = [x_os(end - N_cp_os + 1 : end); x_os];
end
