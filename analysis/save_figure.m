function save_figure(fig, filename)
% 将图形保存为 .fig 和 .png 两种格式
    fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    % 将所有文本对象的解释器设为 none，避免 savefig 时中文触发 tex 解析错误
    all_text = findall(fig, 'Type', 'text');
    for k = 1:numel(all_text)
        set(all_text(k), 'Interpreter', 'none');
    end
    all_axes = findall(fig, 'Type', 'axes');
    for k = 1:numel(all_axes)
        set(all_axes(k), 'TickLabelInterpreter', 'none');
    end
    all_leg = findall(fig, 'Type', 'legend');
    for k = 1:numel(all_leg)
        set(all_leg(k), 'Interpreter', 'none');
    end

    savefig(fig, fullfile(fig_dir, [filename, '.fig']));
    exportgraphics(fig, fullfile(fig_dir, [filename, '.png']), 'Resolution', 300);
    fprintf('已保存：%s（.fig + .png）\n', filename);
end
