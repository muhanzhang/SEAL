function [auc, precision] = WLNM(train_mix, test, K, ith_experiment)   
%  Usage: the main program for Weisfeiler-Lehman Neural Machine (WLNM)
%  --Input--
%  -train_mix: a struct where train_mix.pos contains indices of positive train 
%              links, train_mix.neg contains indices of negative train links, 
%              train_mix.train is a sparse adjacency matrix of observed network 
%              (1: link, 0: otherwise)
%  -test: a struct where test.pos contains indices of positive test links, and
%         test.neg contains indices of negative test links
%  -K: number of vertices in an enclosing subgraph
%  -ith_experiment: exp index, for parallel computing
%  --Output--
%  -auc: the AUC score of WLNM
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
    train = train_mix.train;
    train_pos = train_mix.pos; train_neg = train_mix.neg;
    test_pos = test.pos; test_neg = test.neg;

    if nargin < 3
        K = 20;
    end
    if nargin < 4
        ith_experiment = 1;
    end

    [train_data, train_label] = graph2vector(train_pos, train_neg, train, K);
    [test_data, test_label] = graph2vector(test_pos, test_neg, train, K);
    
    % train a model
    model = 2;
    switch model
    case 1  % logistic regression
        addpath('software/liblinear-2.1/matlab');  % need to install liblinear
        train_data = sparse(train_data);
        test_data = sparse(test_data);
        [~, optim_c] = evalc('liblinear_train(train_label, train_data, ''-s 0 -C -q'');');
        model = liblinear_train(train_label, train_data, sprintf('-s 0 -c %d -q', optim_c(1)));
        [~, acc, scores] = liblinear_predict(test_label, test_data, model, '-b 1 -q');
        acc
        l1 = find(model.Label == 1);
        scores = scores(:, l1);
    case 2 % train a feedforward neural network in Torch
        addpath('software/liblinear-2.1/matlab');  % need to install liblinear
        train_data = sparse(train_data);
        test_data = sparse(test_data);
        if exist('tempdata') ~= 7
            !mkdir tempdata
        end
        libsvmwrite(sprintf('tempdata/traindata_%d', ith_experiment), train_label, train_data);
        libsvmwrite(sprintf('tempdata/testdata_%d', ith_experiment), test_label, test_data);  % prepare data
        cmd = sprintf('th nDNN.lua -inputdim %d -ith_experiment %d', K * (K - 1) / 2, ith_experiment);
        system(cmd, '-echo'); 
        scores = load(sprintf('tempdata/test_log_scores_%d.asc', ith_experiment));
        delete(sprintf('tempdata/traindata_%d', ith_experiment));  % to delete temporal train and test data
        delete(sprintf('tempdata/testdata_%d', ith_experiment));
        delete(sprintf('tempdata/test_log_scores_%d.asc', ith_experiment));
    case 3 % train a feedforward neural network in MATLAB
        layers = [imageInputLayer([K*(K-1)/2 1 1], 'Normalization','none')
        fullyConnectedLayer(32)
        reluLayer
        fullyConnectedLayer(32)
        reluLayer
        fullyConnectedLayer(16)
        reluLayer
        fullyConnectedLayer(2)
        softmaxLayer
        classificationLayer];
        opts = trainingOptions('sgdm', 'InitialLearnRate', 0.1, 'MaxEpochs', 200, 'MiniBatchSize', 128, ...
            'LearnRateSchedule','piecewise', 'LearnRateDropFactor', 0.9, 'L2Regularization', 0, ...
            'ExecutionEnvironment', 'cpu');
        net = trainNetwork(reshape(train_data', K*(K-1)/2, 1, 1, size(train_data, 1)), categorical(train_label), layers, opts);
        [~, scores] = classify(net, reshape(test_data', K*(K-1)/2, 1, 1, size(test_data, 1)));
        scores(:, 1) = [];
    end

    % calculate the AUC, average precision
    [~, ~, ~, auc] = perfcurve(test_label', scores', 1);
    auc
    [xs, ys, ~, aucpr] = perfcurve(test_label', scores', 1, 'XCrit', 'reca', 'YCrit', 'prec');
    precision = sum(diff(xs) .* ys(2:end))  % average precision
end
