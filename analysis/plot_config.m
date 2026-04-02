function [colors, markers, line_styles] = plot_config()
% matlab中设置全局绘图样式
% 统一配置字体、字号、网格等绘图默认属性，
% 返回各算法对应的颜色矩阵、标记符号和线型字符串。
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultTextFontSize', 12);
    set(0, 'DefaultAxesXGrid', 'on');
    set(0, 'DefaultAxesYGrid', 'on');
    set(0, 'DefaultLineLineWidth', 1.5);
    set(0, 'DefaultLineMarkerSize', 6);
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    % 各算法对应颜色（RGB）
    colors = [
        0.000, 0.000, 0.000;  % 原始信号
        0.850, 0.325, 0.098;  % 限幅滤波
        0.000, 0.447, 0.741;  % SLM
        0.466, 0.674, 0.188;  % PTS
        0.301, 0.745, 0.933;  % 预留音调
        0.929, 0.694, 0.125;  % μ律压扩
        0.494, 0.184, 0.556;  % Golay 编码
        0.635, 0.078, 0.184;  % DNN
    ];
    markers = {'o', 's', '^', 'v', '>', 'd', 'h', 'p'};
    line_styles = {'-', '--', ':', '-', '--', '-.', '-', ':'};
end
