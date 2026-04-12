function perm_table = golay_build_perm_table(m)
% GOLAY_BUILD_PERM_TABLE  预计算 {1,...,m} 的排列子集作为 Golay coset 代表
%
% 每个排列 π 定义一个 coset，其中所有码字都是 Golay 序列（PMEPR ≤ 2）。
% 排列 π 和 π 的逆序（reverse）给出相同的 coset（Davis-Jedwab Corollary 4），
% 因此不重复的 coset 数为 m!/2。
%
% 输出：取 2^K 个排列（K = floor(log2(m!/2))），以便用整数比特索引。
%
% 输入：
%   m - 子载波指数（N = 2^m）
%
% 输出：
%   perm_table - (2^K x m) 矩阵，每行是一个排列

    all_perms = perms(1:m);  % m! x m
    n_all = size(all_perms, 1);

    % 去重：对每个排列，比较其与逆序版本，只保留字典序较小的
    keep = true(n_all, 1);
    for i = 1:n_all
        if ~keep(i), continue; end
        rev_perm = all_perms(i, end:-1:1);
        for j = i+1:n_all
            if ~keep(j), continue; end
            if isequal(all_perms(j, :), rev_perm)
                keep(j) = false;
                break;
            end
        end
    end
    unique_perms = all_perms(keep, :);

    % 取 2^K 个
    K = floor(log2(size(unique_perms, 1)));
    perm_table = unique_perms(1:2^K, :);
end
