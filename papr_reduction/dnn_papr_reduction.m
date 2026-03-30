function [x_modified, X_modified, net_trained, train_info] = ...
    dnn_papr_reduction(X_train, params, options)
% 基于深度学习自编码器的 PAPR 降低方法，训练一个自编码器，学习频域扰动以降低 PAPR，同时最小化信号失真。
%   输入：
%     X_train - 训练数据：N_fft x N_samples 的频域符号矩阵
%     params  - get_default_params() 返回的参数结构体
%     options - 训练选项结构体，包含以下字段：
%       .alpha      - 扰动幅度预算
%       .lambda_mse - MSE 损失权重
%       .lambda_papr - PAPR 损失权重
%       .epochs     - 训练轮数
%       .lr         - 学习率
%       .batch_size - 小批量大小
%       .mode       - 'train'或 'apply'，默认：'train'
%       .net        - 预训练网络

    if nargin < 3, options = struct(); end
    if ~isfield(options, 'alpha'),       options.alpha = 0.1; end
    if ~isfield(options, 'lambda_mse'),  options.lambda_mse = 1.0; end
    if ~isfield(options, 'lambda_papr'), options.lambda_papr = 0.01; end
    if ~isfield(options, 'epochs'),      options.epochs = 50; end
    if ~isfield(options, 'lr'),          options.lr = 1e-3; end
    if ~isfield(options, 'batch_size'),  options.batch_size = 128; end
    if ~isfield(options, 'mode'),        options.mode = 'train'; end

    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    n_input = 2 * N;
    % 应用模式：直接调用已训练网络
    if strcmp(options.mode, 'apply')
        net = options.net;
        [x_modified, X_modified] = apply_network(X_train, net, params, options.alpha);
        net_trained = net;
        train_info = [];
        return;
    end

    N_samples = size(X_train, 2);
    % 将复数频域数据转换为实数表示
    X_real = [real(X_train); imag(X_train)];

    % 编码器网络结构：输入 -> 512 -> 256 -> 输出，非线性函数tanh 激活
    encoder_layers = [
        featureInputLayer(n_input, 'Name', 'enc_input')
        fullyConnectedLayer(512, 'Name', 'enc_fc1')
        batchNormalizationLayer('Name', 'enc_bn1')
        reluLayer('Name', 'enc_relu1')
        fullyConnectedLayer(256, 'Name', 'enc_fc2')
        batchNormalizationLayer('Name', 'enc_bn2')
        reluLayer('Name', 'enc_relu2')
        fullyConnectedLayer(n_input, 'Name', 'enc_fc3')
        tanhLayer('Name', 'enc_tanh')
    ];

    % 解码器网络结构：输入 -> 512 -> 256 -> 输出，线性激活
    decoder_layers = [
        featureInputLayer(n_input, 'Name', 'dec_input')
        fullyConnectedLayer(512, 'Name', 'dec_fc1')
        batchNormalizationLayer('Name', 'dec_bn1')
        reluLayer('Name', 'dec_relu1')
        fullyConnectedLayer(256, 'Name', 'dec_fc2')
        batchNormalizationLayer('Name', 'dec_bn2')
        reluLayer('Name', 'dec_relu2')
        fullyConnectedLayer(n_input, 'Name', 'dec_fc3')
    ];

    enc_net = dlnetwork(layerGraph(encoder_layers));
    dec_net = dlnetwork(layerGraph(decoder_layers));

    % 初始化训练记录
    train_info.loss_total = zeros(options.epochs, 1);
    train_info.loss_mse   = zeros(options.epochs, 1);
    train_info.loss_papr  = zeros(options.epochs, 1);

    avg_grad_enc = []; avg_sq_enc = [];
    avg_grad_dec = []; avg_sq_dec = [];

    n_batches = ceil(N_samples / options.batch_size);

    % 训练循环
    for epoch = 1:options.epochs
        perm = randperm(N_samples);
        epoch_loss = 0;
        epoch_mse = 0;
        epoch_papr = 0;

        for b = 1:n_batches
            idx = perm((b-1)*options.batch_size+1 : min(b*options.batch_size, N_samples));
            X_batch = dlarray(X_real(:, idx), 'CB');

            [loss, grad_enc, grad_dec, mse_val, papr_val] = ...
                dlfeval(@model_loss, enc_net, dec_net, X_batch, ...
                        params, options);

            % Adam 优化器更新编码器和解码器参数
            [enc_net, avg_grad_enc, avg_sq_enc] = ...
                adamupdate(enc_net, grad_enc, avg_grad_enc, avg_sq_enc, ...
                           (epoch-1)*n_batches+b, options.lr);
            [dec_net, avg_grad_dec, avg_sq_dec] = ...
                adamupdate(dec_net, grad_dec, avg_grad_dec, avg_sq_dec, ...
                           (epoch-1)*n_batches+b, options.lr);

            epoch_loss = epoch_loss + extractdata(loss);
            epoch_mse  = epoch_mse + extractdata(mse_val);
            epoch_papr = epoch_papr + extractdata(papr_val);
        end

        train_info.loss_total(epoch) = epoch_loss / n_batches;
        train_info.loss_mse(epoch)   = epoch_mse / n_batches;
        train_info.loss_papr(epoch)  = epoch_papr / n_batches;

        if mod(epoch, 10) == 0
            fprintf('Epoch %d/%d: 总损失=%.4f, MSE=%.4f, PAPR=%.4f\n', ...
                epoch, options.epochs, train_info.loss_total(epoch), ...
                train_info.loss_mse(epoch), train_info.loss_papr(epoch));
        end
    end

    net_trained.encoder = enc_net;
    net_trained.decoder = dec_net;

    [x_modified, X_modified] = apply_network(X_train, net_trained, params, options.alpha);
end

function [loss, grad_enc, grad_dec, mse_val, papr_val] = ...
    model_loss(enc_net, dec_net, X_batch, params, options)
% 自定义损失函数，包含可微分的 PAPR 替代损失
    N = params.N_fft;
    alpha = options.alpha;
    % 编码器：计算扰动，tanh 输出范围 [-1,1]
    perturbation = forward(enc_net, X_batch);
    X_modified = X_batch + alpha * perturbation;
    % 解码器：重建原始信号
    X_recovered = forward(dec_net, X_modified);
    % MSE 损失
    mse_val = mean((X_batch - X_recovered).^2, 'all');
    % 用可微分 PAPR 替代损失
    % 鼓励频域信号幅度分布更均匀
    p = 10;
    X_mag_sq = X_modified(1:N, :).^2 + X_modified(N+1:end, :).^2;  % 每子载波 |X|^2
    % p 范数近似峰值
    peak_approx = (mean(X_mag_sq.^p, 1)).^(1/p);
    avg_power = mean(X_mag_sq, 1);
    papr_surrogate = mean(peak_approx ./ (avg_power + 1e-8));
    papr_val = papr_surrogate;
    % 联合损失
    loss = options.lambda_mse * mse_val + options.lambda_papr * papr_val;
    % 计算梯度
    grad_enc = dlgradient(loss, enc_net.Learnables);
    grad_dec = dlgradient(loss, dec_net.Learnables);
end

function [x_out, X_out] = apply_network(X_data, net, params, alpha)
% 对输入数据应用已训练的编码器网络，生成修改后的频域及时域信号
    N = params.N_fft;
    L = params.L;
    N_os = N * L;
    N_samples = size(X_data, 2);

    X_real = [real(X_data); imag(X_data)];
    X_dl = dlarray(X_real, 'CB');

    perturbation = predict(net.encoder, X_dl);
    X_mod = X_dl + alpha * perturbation;
    X_mod_data = extractdata(X_mod);

    X_out = X_mod_data(1:N, :) + 1j * X_mod_data(N+1:end, :);
    x_out = zeros(N_os, N_samples);

    for i = 1:N_samples
        X_padded = [X_out(1:N/2, i); zeros((L-1)*N, 1); X_out(N/2+1:end, i)];
        x_out(:, i) = sqrt(N_os) * ifft(X_padded, N_os);
    end
end
