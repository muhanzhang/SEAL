function [node_embeddings] = generate_embeddings(A, data_name, emd_method)
%  Usage: generate node embeddings to support SEAL.m and embedding_lp.m
%  --Input--
%  -A: the observed network's adjacency matrix from which to generate 
%      node embeddings
%  -data_name: the name of the dataset
%  -emd_method: 'node2vec', 'LINE', 'SPC', default 'node2vec'
%  --Output--
%  -node_embeddings: a matrix, ith row contains the ith node's embeddings
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%

if nargin < 3
    emd_method = 'node2vec'
end

[i, j] = find(triu(A));
train = [i, j];

mkdir data/embedding;

switch emd_method
case 'LINE'
    cd data/embedding;
    train = [train; [train(:, 2), train(:, 1)]];  % LINE supports only directed edges, thus double the train
    dlmwrite(strcat(data_name, '.edgelist'), train, ' ');  % convert train to edgelist which will be read by node2vec
    cmd = sprintf('./../../software/LINE/linux/line -train %s.edgelist -output %s.emd -size 128 -samples 100 -threads 20', data_name, data_name);
    system(cmd);
    cd ../..;

    node_embeddings = dlmread(['data/embedding/', data_name, '.emd']);
    node_embeddings(1, :) = [];
    tmp = node_embeddings(:, 1); % the node indices with embeddings
    node_embeddings(:, 1) = [];
    avg_embedding = mean(node_embeddings, 1);
    node_embeddings(tmp, :) = node_embeddings;
    nodes_missing_embedding = setdiff([1:size(A, 1)], tmp);
    node_embeddings(nodes_missing_embedding, :) = repmat(avg_embedding, length(nodes_missing_embedding), 1);  % replace missing rows with average embeddings

    assert(size(node_embeddings, 1) == size(A, 1));  % ensure all nodes have embeddings

case 'node2vec'
    cd data/embedding;
    dlmwrite(strcat(data_name, '.edgelist'), train, ' ');  % convert train to edgelist which will be read by node2vec

    if strfind(data_name, 'PPI_subgraph')  % use node2vec paper's provided hyperparameters for specific datasets
        p = 4; q = 1;
    elseif strfind(data_name, 'BlogCatalog')
        p = 0.25; q = 0.25;
    elseif strfind(data_name, 'Wikipedia')
        p = 4; q = 0.5;
    else  % otherwise, use the default
        p = 1; q = 1;
    end

    cmd = sprintf('python ../../software/node2vec/src/main.py --input %s.edgelist --output %s.emd --p %f --q %f --dimensions 128 --window-size 10', data_name, data_name, p, q);
    system(cmd);
    cd ../..;

    node_embeddings = dlmread(['data/embedding/', data_name, '.emd']);
    node_embeddings(1, :) = [];
    tmp = node_embeddings(:, 1); % the node indices with embeddings
    node_embeddings(:, 1) = [];
    avg_embedding = mean(node_embeddings, 1);
    node_embeddings(tmp, :) = node_embeddings;
    nodes_missing_embedding = setdiff([1:size(A, 1)], tmp);
    node_embeddings(nodes_missing_embedding, :) = repmat(avg_embedding, length(nodes_missing_embedding), 1);  % replace missing rows with average embeddings
    assert(size(node_embeddings, 1) == size(A, 1));  % ensure all nodes have embeddings

case 'SPC'  % spectral clustering
    D_inv = full(diag((sum(A, 2)).^(-0.5)));

    D_inv(isnan(D_inv)) = 0;
    D_inv(isinf(D_inv)) = 0;
    L = eye(size(A)) - D_inv * A * D_inv;
    %[V, ~] = eig(L); % randomly select top 128
    %node_embeddings = V(:, 1:128);
    [V, D] = eigs(L, size(L, 1), 0);
    assert(D(size(D, 1), size(D, 1)) >= D(1, 1));  % ensure eigs ranked smaller to larger
    node_embeddings = V(:, 1:128);
    assert(size(node_embeddings, 1) == size(A, 1));  % ensure all nodes have embeddings

case 'MF'  % needs to first run method MF.m which save embeddings to tempdata/
    ith_experiment = str2num(data_name(end));
    V = dlmread(['tempdata/FMmodelv_', num2str(ith_experiment), '.txt']);
    W = dlmread(['tempdata/FMmodelw_', num2str(ith_experiment), '.txt']);
    node_embeddings = [W, V'];
    assert(size(node_embeddings, 1) == size(A, 1));  % ensure all nodes have embeddings

case 'none'
    node_embeddings = [];

end


