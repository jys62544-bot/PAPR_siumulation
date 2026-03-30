function [papr_axis, ccdf_values] = compute_ccdf(papr_dB_vec)
% 根据上课讲的PAPR的定义，可视化 PAPR 的互补累积分布函数
%   输入：
%     papr_dB_vec  - 各 OFDM 符号的 PAPR 值向量（dB）
%   输出：
%     papr_axis    - 排序后的 PAPR 值
%     ccdf_values  - 对应的超越概率
    papr_sorted = sort(papr_dB_vec);
    N = length(papr_sorted);
    ccdf_values = (N:-1:1).' / N;
    papr_axis = papr_sorted(:);
end
