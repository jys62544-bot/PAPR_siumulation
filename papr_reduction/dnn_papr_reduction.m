function [net_trained, train_info] = dnn_papr_reduction(X_train, params, options, mod_map, H_ch)
% DNN_PAPR_REDUCTION  基于"限幅 + DNN 逐子载波补偿"的 PAPR 降低方法。
% 发射端限幅降 PAPR，接收端 DNN 逐子载波补偿限幅失真。
% 网络输入：单个子载波的 (I, Q, mod_order_norm, H_norm)，输出：补偿后的 (I, Q)。
% 所有子载波共享同一网络，训练数据量 = N_samples * N_fft。
%
% 训练时加入信道效应和多 SNR 噪声，使 DNN 在推理时能适应实际信道条件。
%
%   输入：
%     X_train - N_fft x N_samples 频域符号矩阵
%     params  - 参数结构体
%     options - 训练选项
%     mod_map - N_fft x 1 调制阶数向量（自适应模式），或 []（非自适应模式）
%     H_ch    - N_fft x 1 信道频率响应（自适应模式），或 []（非自适应模式）
%
%   输出：
%     net_trained - 结构体，含 .net (dlnetwork) 和 .CR (限幅比)
%     train_info  - 训练损失记录

    if nargin < 3, options = struct(); end
    if nargin < 4, mod_map = []; end
    if nargin < 5, H_ch = []; end
    if ~isfield(options, 'CR'),         options.CR = 1.2; end
    if ~isfield(options, 'epochs'),     options.epochs = 30; end
    if ~isfield(options, 'lr'),         options.lr = 1e-3; end
    if ~isfield(options, 'batch_size'), options.batch_size = 1024; end

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    CR = options.CR;
    N_samples = size(X_train, 2);

    %% 生成逐子载波训练对（含信道效应）
    fprintf('  生成逐子载波训练对 (CR=%.2f, 含信道效应)...\n', CR);

    % 获取信道时域冲激响应（用于卷积）
    if ~isempty(H_ch)
        % 从 H_ch 反推 h_ch 用于训练数据生成
        H_full = zeros(N_os, 1);
        H_full(1:N/2) = H_ch(1:N/2);
        H_full(N_os - N/2 + 1:end) = H_ch(N/2+1:end);
        h_ch = ifft(H_full, N_os);
    else
        h_ch = [];
    end

    % 训练时使用多种 SNR，覆盖实际推理条件
    snr_train_range = [0, 5, 10, 15, 20, 25, 30];
    n_snr = length(snr_train_range);

    X_rx_clip_all = zeros(N, N_samples);  % 接收端均衡后的限幅信号
    X_orig_all = X_train;                 % 标签：原始无失真频域符号

    for i = 1:N_samples
        X_i = X_train(:, i);
        % 发射端：过采样 IFFT -> 限幅 -> 加 CP
        X_padded = [X_i(1:N/2); zeros((L-1)*N, 1); X_i(N/2+1:end)];
        x_os = sqrt(N_os) * ifft(X_padded, N_os);
        A = CR * rms(x_os);
        x_clip = min(abs(x_os), A) .* exp(1j * angle(x_os));

        % 加循环前缀
        N_cp_os = params.N_cp * L;
        x_cp = [x_clip(end - N_cp_os + 1 : end); x_clip];

        % 选择当前样本使用的 SNR
        snr_idx = mod(i - 1, n_snr) + 1;
        snr_i = snr_train_range(snr_idx);

        % 通过信道
        if ~isempty(h_ch)
            y_mp = conv(x_cp, h_ch);
            y_mp = y_mp(1:length(x_cp));
        else
            y_mp = x_cp;
        end

        % 添加 AWGN
        sig_pow = mean(abs(y_mp).^2);
        noise_pow = sig_pow / (10^(snr_i / 10));
        noise = sqrt(noise_pow / 2) * (randn(size(y_mp)) + 1j * randn(size(y_mp)));
        y_noisy = y_mp + noise;

        % 接收端：去 CP -> FFT -> 信道均衡
        y_os = y_noisy(N_cp_os + 1 : N_cp_os + N_os);
        Y_full = fft(y_os, N_os) / sqrt(N_os);
        X_hat = [Y_full(1:N/2); Y_full(N_os - N/2 + 1:end)];

        if ~isempty(H_ch)
            X_hat = X_hat ./ H_ch;
        end

        X_rx_clip_all(:, i) = X_hat;
    end

    % 构造辅助特征
    if ~isempty(mod_map)
        mod_feat = single(mod_map / 8);
    else
        mod_feat = single(ones(N, 1) * params.bps / 8);
    end
    if ~isempty(H_ch)
        h_feat = single(abs(H_ch) / max(abs(H_ch)));
    else
        h_feat = single(ones(N, 1));
    end

    % 将所有子载波展平为逐子载波训练对
    total_sc = N * N_samples;
    X_in_sc  = zeros(4, total_sc, 'single');
    X_lbl_sc = zeros(2, total_sc, 'single');

    for i = 1:N_samples
        idx = (i-1)*N + (1:N);
        X_in_sc(1, idx) = real(X_rx_clip_all(:, i));
        X_in_sc(2, idx) = imag(X_rx_clip_all(:, i));
        X_in_sc(3, idx) = mod_feat;
        X_in_sc(4, idx) = h_feat;
        X_lbl_sc(1, idx) = real(X_train(:, i));
        X_lbl_sc(2, idx) = imag(X_train(:, i));
    end

    distortion_mse = mean((X_in_sc(1:2,:) - X_lbl_sc).^2, 'all');
    fprintf('  逐子载波失真 MSE = %.6f, 总样本 = %d\n', distortion_mse, total_sc);

    %% 构建 DNN
    layers = [
        featureInputLayer(4, 'Name', 'input', 'Normalization', 'none')
        fullyConnectedLayer(128, 'Name', 'fc1')
        reluLayer('Name', 'relu1')
        fullyConnectedLayer(64, 'Name', 'fc2')
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
