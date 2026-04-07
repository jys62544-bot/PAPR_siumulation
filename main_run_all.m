function main_run_all()
% MAIN_RUN_ALL 按顺序运行全部 6 个仿真实验。
%   结果保存至 results/figures/ 和 results/data/ 目录。
%   任何实验失败会立即中止，避免后续实验使用过期数据。

    clc; close all;
    fprintf('  OFDM PAPR系统仿真\n');
    fprintf('  %s\n', datetime('now'));
    total_tic = tic;
    run_experiment(1, 6, 'sim_01_ccdf_comparison');
    run_experiment(2, 6, 'sim_02_ber_comparison');
    run_experiment(3, 6, 'sim_03_psd_analysis');
    run_experiment(4, 6, 'sim_04_complexity');
    run_experiment(5, 6, 'sim_05_param_sweep');
    run_experiment(6, 6, 'sim_06_dnn_papr');

    total_time = toc(total_tic);
    fprintf('  全部实验完成！\n');
    fprintf('  总耗时：%.1f 秒（%.1f 分钟）\n', total_time, total_time/60);
    fprintf('  图形结果：results/figures/\n');
    fprintf('  数据结果：results/data/\n');
end

function run_experiment(num, total, script_name)
% 运行单个实验脚本，失败时中止整个流程
    fprintf('[%d/%d] 正在运行 %s...\n', num, total, script_name);
    run(script_name);
    close all;
    fprintf('\n');
end
