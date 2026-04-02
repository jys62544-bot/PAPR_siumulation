function [net_trained, train_info] = dnn_papr_reduction(X_train, params, options)
% DNN_PAPR_REDUCTION  基于"限幅 + DNN 逐子载波补偿"的 PAPR 降低方法。
% 发射端限幅降 PAPR，接收端 DNN 逐子载波补偿限幅失真。
% 网络输入：单个子载波的 (I, Q)，输出：补偿后的 (I, Q)。
% 所有子载波共享同一网络，训练数据量 = N_samples * N_fft。
%
%   输入：
%     X_train - N_fft x N_samples 频域符号矩阵
%     params  - 参数结构体
%     options - 训练选项
%
%   输出：
%     net_trained - 结构体，含 .net (dlnetwork) 和 .CR (限幅比)
%     train_info  - 训练损失记录

    if nargin < 3, options = struct(); end
    if ~isfield(options, 'CR'),         options.CR = 1.2; end
    if ~isfield(options, 'epochs'),     options.epochs = 30; end
    if ~isfield(options, 'lr'),         options.lr = 1e-3; end
    if ~isfield(options, 'batch_size'), options.batch_size = 1024; end

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    CR = options.CR;
    N_samples = size(X_train, 2);

    %% 生成逐子载波训练对
    fprintf('  生成逐子载波训练对 (CR=%.2f)...\n', CR);

    % 批量计算限幅后的频域符号
    X_clip_all = zeros(N, N_samples);
    for i = 1:N_samples
        X_i = X_train(:, i);
        X_padded = [X_i(1:N/2); zeros((L-1)*N, 1); X_i(N/2+1:end)];
        x_os = sqrt(N_os) * ifft(X_padded, N_os);
        A = CR * rms(x_os);
        x_clip = min(abs(x_os), A) .* exp(1j * angle(x_os));
        Y_full = fft(x_clip, N_os) / sqrt(N_os);
        X_clip_all(:, i) = [Y_full(1:N/2); Y_full(N_os - N/2 + 1:end)];
    end

    % 将所有子载波展平为逐子载波训练对
    % 输入：限幅后的 (I, Q)，标签：原始 (I, Q)
    total_sc = N * N_samples;  % 总子载波数
    X_in_sc  = zeros(2, total_sc, 'single');  % 限幅后
    X_lbl_sc = zeros(2, total_sc, 'single');  % 原始

    for i = 1:N_samples
        idx = (i-1)*N + (1:N);
        X_in_sc(1, idx) = real(X_clip_all(:, i));
        X_in_sc(2, idx) = imag(X_clip_all(:, i));
        X_lbl_sc(1, idx) = real(X_train(:, i));
        X_lbl_sc(2, idx) = imag(X_train(:, i));
    end

    distortion_mse = mean((X_in_sc(:) - X_lbl_sc(:)).^2);
    fprintf('  逐子载波失真 MSE = %.6f, 总样本 = %d\n', distortion_mse, total_sc);

    %% 构建 DNN：小型逐子载波网络
    layers = [
        featureInputLayer(2, 'Name', 'input', 'Normalization', 'none')
        fullyConnectedLayer(64, 'Name', 'fc1')
        reluLayer('Name', 'relu1')
        fullyConnectedLayer(32, 'Name', 'fc2')
        reluLayer('Name', 'relu2')
        fullyConnectedLayer(2, 'Name', 'fc_out')
    ];
    net = dlnetwork(layerGraph(layers));

    %% 训练循环
    train_info.loss_total = zeros(options.epochs, 1);
    avg_grad = []; avg_sq = [];
    n_batches = ceil(total_sc / options.batch_size);

    fprintf('  训练逐子载波 DNN（%d 样本, %d 轮）...\n', total_sc, options.epochs);

    for epoch = 1:options.epochs
        perm = randperm(total_sc);
        epoch_loss = 0;

        for b = 1:n_batches
            b_idx = perm((b-1)*options.batch_size+1 : min(b*options.batch_size, total_sc));

            X_in  = dlarray(X_in_sc(:, b_idx), 'CB');
            X_lbl = dlarray(X_lbl_sc(:, b_idx), 'CB');

            [loss, grads] = dlfeval(@mse_loss, net, X_in, X_lbl);

            [net, avg_grad, avg_sq] = adamupdate(net, grads, avg_grad, avg_sq, ...
                (epoch-1)*n_batches+b, options.lr);

            epoch_loss = epoch_loss + double(extractdata(loss));
        end

        train_info.loss_total(epoch) = epoch_loss / n_batches;

        if mod(epoch, 5) == 0
            fprintf('  Epoch %3d/%d: MSE = %.6f\n', epoch, options.epochs, ...
                train_info.loss_total(epoch));
        end
    end

    net_out.net = net;
    net_out.CR = CR;
    net_trained = net_out;
    fprintf('  训练完成。最终 MSE = %.6f（失真上界 = %.6f）\n', ...
        train_info.loss_total(end), distortion_mse);
end

function [loss, grads] = mse_loss(net, X_in, X_lbl)
    X_pred = forward(net, X_in);
    loss = mean((X_pred - X_lbl).^2, 'all');
    grads = dlgradient(loss, net.Learnables);
end
