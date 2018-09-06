function [auc, precision] = SEAL(train_mix, test, h, include_embedding, include_attribute, ith_experiment)   
%  Usage: the main program of SEAL (learning from Subgraphs, Embeddings, and Attributes for Link prediction)
%  --Input--
%  -train_mix: a struct where train_mix.pos contains indices [(i1, j1); (i2, j2); ...] 
%              of positive train links, train_mix.neg contains indices of negative 
%              train links, train_mix.train is a sparse adjacency matrix of observed 
%              network (1: link, 0: otherwise), train_mix.data_name is dataset name
%  -test: a struct where test.pos contains indices of positive test links, and
%         test.neg contains indices of negative test links
%  -h: maximum hop to extract enclosing subgraphs, h='auto' selects h from {1, 2}
%  -include_embedding: 1 to include node embeddings into node information matrix, default 0
%  -include_attribute: 1 to include node attributes into node information matrix, default 0
%  -ith_experiment: exp index, for parallel computing, default 1
%  --Output--
%  -auc: the AUC score on testing links
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
    A = train_mix.train;  % the observed network
    data_name = train_mix.data_name;
    train_pos = train_mix.pos; train_neg = train_mix.neg;  % the indices of observed links used as training data
    test_pos = test.pos; test_neg = test.neg;  % the indices of unobserved links used as testing data
    
    if nargin < 3
        h = 1;
    end

    if h == 'auto'
    % randomly sample 10% train links as validation links, 
    % select h from {1, 2} based on the validation performance of AA and CN
        [val_train, val_test] = DivideNet(A, 0.9, false);
        h_val_train = triu(sparse(val_train), 1);
        h_val_test = triu(sparse(val_test), 1);
        [~, ~, val_test_pos, val_test_neg] = sample_neg(h_val_train, h_val_test, 1, 1);
        val_test = {};
        val_test.pos = val_test_pos;
        val_test.neg = val_test_neg;
        val_auc_AA = AA(val_train, val_test)
        val_auc_CN = CN(val_train, val_test)
        if val_auc_AA >= val_auc_CN
            h = 2;
            disp(['Choose h=', num2str(h)])
        else
            h = 1;
            disp(['Choose h=', num2str(h)])
        end
    end

    if nargin < 4
        include_embedding = 0;
    end
    if nargin < 5
        include_attribute = 0;
    end
    if nargin < 6
        ith_experiment = 1;
    end

    DGCNN_path = '../../DGCNN/';  % SEAL by default uses DGCNN for graph classification, change to your DGCNN path
    data_name_i = [data_name, '_', num2str(ith_experiment)];  % the ith experiment's temporal data name
    
    train_size = size(train_pos, 1) + size(train_neg, 1)
    test_size = size(test_pos, 1) + size(test_neg, 1)

    % extract enclosing subgraphs
    [data, max_size] = graph2mat([train_pos; train_neg], [test_pos; test_neg], A, h, ith_experiment, 0, data_name, include_embedding, include_attribute);
    label = [ones(size(train_pos, 1), 1); zeros(size(train_neg, 1), 1); ...
             ones(size(test_pos, 1), 1); zeros(size(test_neg, 1), 1)];  % graph labels (classes), not to confuse with node labels

    % permutate the train set, because in DGCNN we do not permutate
    perm = randperm(train_size);
    data(:, 1:train_size) = data(:, perm);
    label(1:train_size) = label(perm);

    % save to tempdata/
    data_info = whos('data');
    data_bytes = data_info.bytes
    n_splits = ceil(data_bytes / 1e8);  % split data into < 1GB splits, so that they can be saved in default .mat (torch only reads default .mat and -v7.3 is not supported)
    split_size = ceil(length(label) / n_splits);
    system('mkdir tempdata')
    system(sprintf('mkdir tempdata/%s', data_name_i));
    system(sprintf('rm tempdata/%s/*', data_name_i));
    for i = 1: n_splits
        display(sprintf('Saving split-%d of %d...', i, n_splits))
        data_split = data((i-1)*split_size+1: min(i*split_size, length(label)));
        if size(data_split, 2) == 0  % no more split
            break
        end
        save_struct.(data_name) = data_split;
        save(['tempdata/', data_name_i, '/split_', num2str(i), '.mat'], '-struct', 'save_struct');
        save_struct2.(['l' lower(data_name)]) = label((i-1)*split_size+1: min(i*split_size, length(label)));
        save(['tempdata/', data_name_i, '/split_', num2str(i), '.mat'], '-struct', 'save_struct2', '-append');
    end

    clear data

    % convert .mat to .dat format (for Torch-based DGCNN to read)
    system(sprintf('th generate_torch_graphs.lua -dataName %s -ith_experiment %d', data_name, ith_experiment));

    % run DGCNN
    %gpu = mod(ith_experiment - 1, 4) + 1;  % determine gpu ID, assuming you have 4 gpus
    gpu = 1;

    data_pos = ['tempdata/', data_name_i];
    if include_embedding == 0 && include_attribute == 0
        cmd = sprintf('th %smain.lua -dataPos %s -dataName %s -maxNodeLabel %d -printAUC -save tempdata -testNumber %d -gpu %d -maxEpoch 50 -fixed_shuffle original -testAfterAll', DGCNN_path, data_pos, data_name_i, max_size, test_size, gpu)
    else
        cmd = sprintf('th %smain.lua -dataPos %s -dataName %s -printAUC -save tempdata -testNumber %d -gpu %d -maxEpoch 50 -fixed_shuffle original -nodeLabel original -inputChannel %d -testAfterAll', DGCNN_path, data_pos, data_name_i, test_size, gpu, max_size)
    end
    system(cmd);
    
    auc = load(['tempdata/', data_name_i, '/finalAUC'])
    scores = load(['tempdata/', data_name_i, '/scores']);
    test_label = label(train_size+1:end);
    %[~, top_idx] = sort(scores, 'descend');
    %precision = nnz(test_label(top_idx(1:100))) / 100  % calculate precision of top 10% predictions
    [xs, ys, ~, aucpr] = perfcurve(test_label', scores', 1, 'XCrit', 'reca', 'YCrit', 'prec');
    precision = sum(diff(xs) .* ys(2:end))  % average precision
end
