function [auc, precision] = ensemble_heuristics(train_mix, test, all_sims)   
%  Usage: Train a logistic regression over various heuristics
%  --Input--
%  -train_mix: a struct where train_mix.pos contains indices [(i1, j1); (i2, j2); ...] 
%              of positive train links, train_mix.neg contains indices of negative 
%              train links, train_mix.train is a sparse adjacency matrix of observed 
%              network (1: link, 0: otherwise)
%  -test: a struct where test.pos contains indices of positive test links, and
%         test.neg contains indices of negative test links
%  -all_sims: similarity scores of all heuristics
%  --Output--
%  -auc: the AUC score on testing links
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%


A = train_mix.train;
data_name = train_mix.data_name;
train_pos = train_mix.pos; train_neg = train_mix.neg;
test_pos = test.pos; test_neg = test.neg;

train_idxs = [train_pos; train_neg];
test_idxs = [test_pos; test_neg];

% construct training/testing data using heuristic scores
num_heuristics = length(all_sims);
train_data = zeros(size(train_idxs, 1), num_heuristics);
test_data = zeros(size(test_idxs, 1), num_heuristics);
idxs_train = sub2ind(size(A), train_idxs(:, 1), train_idxs(:, 2));
idxs_test = sub2ind(size(A), test_idxs(:, 1), test_idxs(:, 2));
for i = 1: num_heuristics
    sim = all_sims{i};
    train_data(:, i) = sim(idxs_train);
    test_data(:, i) = sim(idxs_test);

train_data = sparse(train_data);
train_label = [ones(size(train_pos, 1), 1); zeros(size(train_neg, 1), 1)];

test_data = sparse(test_data);
test_label = [ones(size(test_pos, 1), 1); zeros(size(test_neg, 1), 1)];

% train a model, default LR
model = 1;
switch model
case 1  % logistic regression
    addpath('software/liblinear-2.1/matlab');  % need to install liblinear
    [~, optim_c] = evalc('liblinear_train(train_label, train_data, ''-s 0 -C -q'');');
    optim_c
    model = liblinear_train(train_label, train_data, sprintf('-s 0 -c %d -q', optim_c(1)));
    tic
    [~, acc, scores] = liblinear_predict(test_label, test_data, model, '-b 1 -q');
    display(sprintf('Final inference time on test set: %.4f', toc))
    acc
    l1 = find(model.Label == 1);
    scores = scores(:, l1);
case 2 % train a feedforward neural network in Torch
    addpath('software/liblinear-2.1/matlab');  % need to install liblinear
    train_data = sparse(train_data);
    train_data(:, 1) = train_data(:, 1) + eps;  % avoid zeros in the first column (otherwise libsvmwrite will ignore)
    test_data = sparse(test_data);
    test_data(:, 1) = test_data(:, 1) + eps;  % avoid zeros in the first column (otherwise libsvmwrite will ignore)
    if exist('tempdata') ~= 7
        !mkdir tempdata
    end
    libsvmwrite(sprintf('tempdata/traindata_%d', ith_experiment), train_label, train_data);
    libsvmwrite(sprintf('tempdata/testdata_%d', ith_experiment), test_label, test_data);  % prepare data
    cmd = sprintf('th nDNN.lua -inputdim %d -ith_experiment %d', emd_dim, ith_experiment);
    system(cmd, '-echo'); 
    scores = load(sprintf('tempdata/test_log_scores_%d.asc', ith_experiment));
    delete(sprintf('tempdata/traindata_%d', ith_experiment));  % to delete temporal train and test data
    delete(sprintf('tempdata/testdata_%d', ith_experiment));
    delete(sprintf('tempdata/test_log_scores_%d.asc', ith_experiment));
case 3 % train a feedforward neural network in MATLAB
    layers = [imageInputLayer([emd_dim 1 1], 'Normalization','none')
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
    net = trainNetwork(reshape(train_data', emd_dim, 1, 1, size(train_data, 1)), categorical(train_label), layers, opts);
    [~, scores] = classify(net, reshape(test_data', emb_dim, 1, 1, size(test_data, 1)));
    scores(:, 1) = [];
end

% calculate the AUC
[~, ~, ~, auc] = perfcurve(test_label', scores', 1);
auc
[xs, ys, ~, aucpr] = perfcurve(test_label', scores', 1, 'XCrit', 'reca', 'YCrit', 'prec');
precision = sum(diff(xs) .* ys(2:end))  % average precision
end
