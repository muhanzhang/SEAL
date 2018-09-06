function [data, max_size] = graph2mat(train, test, A, h, ith_experiment, for_graph_kernel, data_name, include_embedding, include_attribute)
%  Usage: to extract links' enclosing subgraphs saved as adjacency
%         matrices, and generate their node information matrices
%         (if include_embedding==1 or include_attribute==1, 
%         node information matrix will contain one-hot encoding
%         of node label + node embedding vector + node
%         attribute vector; otherwise will only contain "integer"
%         node labels; used by SEAL
%         
%  --Input--
%  -train: indices of training links
%  -test: indices of testing links
%  -A: the observed network's adjacency matrix from which to
%      to extract enclosing subgraphs
%  -h: the maximum hop to extract enclosing subgraphs
%  -ith_experiment: exp index, default 1
%  -for_graph_kernel: if 1, node adajacency lists will be
%                     stored (needed by some graph kernels)
%  -data_name: the name of the dataset
%  -include_embedding: if 1, node embeddings are included
%  -include_attribute: if 1, node attributes are included
%  --Output--
%  -data: a collection of graphs in the WL kernel format
%  -max_size: either 1) maximum node label (if no embedding
%             or attribute are included, but only node label)
%             or 2) length of one-hot encoding of node label
%             + node embedding + node attribute
%
%  *author: Muhan Zhang, Washington University in St. Louis

if nargin < 8
    include_embedding = 0;
end
if nargin < 9
    include_attribute = 0;
end

all = [train; test];
train_size = size(train, 1);
test_size = size(test, 1);
all_size = train_size + test_size;

% Extract subgraphs
data = {};  % used to store graph data
one_tenth = floor(all_size / 10);
max_size = zeros(1, all_size);

data_name_i = [data_name, '_', num2str(ith_experiment)];

% generate node embeddings
node_information = [];
emd_method = 'node2vec';
if include_embedding == 1
    % negative injection to avoid overfitting
    A1 = A;
    train_neg = train(size(train, 1)/2+1:end, :);
    train_neg = [train_neg; [train_neg(:, 2), train_neg(:, 1)]];
    idx_train_neg = sub2ind(size(A1), train_neg(:, 1), train_neg(:, 2));
    A1(idx_train_neg) = 1;

    % generate the node embeddings from the injected network
    node_embeddings = generate_embeddings(A1, data_name_i, emd_method);

    % whether to include some global node centrality features (increase performance on some networks)
    for include_global_node_topological_features = 1:0
        % add some global node features into the embeddings
        G = graph(A1);
        deg = centrality(G, 'degree');
        closeness = centrality(G, 'closeness');
        betweenness = centrality(G, 'betweenness');
        pr = centrality(G, 'pagerank');
        eig = centrality(G, 'eigenvector');
        node_centralities = [deg, closeness, betweenness, pr, eig];
        node_centralities = zscore(node_centralities);  % standardization
        node_embeddings = [node_embeddings, node_centralities];
    end

    clear A1;
    node_information = [node_information, node_embeddings];
end

% load node attributes
if include_attribute == 1
    load(['data/', data_name, '.mat']);  % use node classes as node attributes, assume data_name.mat contains a variable 'group' that stores the attributes
    node_attributes = full(group);
    node_information = [node_information, node_attributes];
end

% now begin enclosing subgraph extraction
display('Subgraph Extraction Begins...')
tic;

% determine whether outer loop is using parallel computing, if not, parallelize here
worker = getCurrentWorker;  % to check if the current file is running in a worker session
if exist('worker')  % worker exists means that the outer loop is running parallelly
    workers = 0;  % then set the number of workers here to 0
else
    workers = parpool(feature('numcores'));  % otherwise, parallelize here
end

parfor (i = 1: all_size, workers)
    ind = all(i, :);
    if include_embedding == 1 || include_attribute == 1
        [sample, max_nl_size] = subgraph2mat(ind, A, h, node_information);
    else
        [sample, max_nl_size] = subgraph2mat(ind, A, h);
    end

    max_size(i) = max_nl_size;
    
    data(i).am = sample.am;
    data(i).nl = sample.nl;

    % if extracting subgraphs for graph kernel, include adjacency list in data
    if for_graph_kernel == 1
        data(i).al = sample.al;
    end
   
    %i/all_size  % display finer progress
    progress = i / one_tenth;
    if ismember(progress, [1:10])
        display(sprintf('Subgraph Extraction Progress %d0%%...', progress));
    end
end
if exist('poolobj')
    delete(poolobj)
end
max_size = max(max_size);

% convert integer node labels to one-hot embeddings
if include_embedding == 1 || include_attribute == 1
    for i = 1: all_size
        tmp = data(i).nl.values;
        tmp = tmp(:, 1);  % integer node labels
        tmp_onehot = zeros(size(tmp, 1), max_size);  % the one-hot embedding matrix
        tmp_onehot(sub2ind(size(tmp_onehot), [1:size(tmp, 1)]', tmp)) = 1;
        data(i).nl.values = [tmp_onehot, data(i).nl.values(:, 2:end)];
    end
    max_size = size(data(1).nl.values, 2);  % new size
end
toc;
end


function [sample, max_nl_size] = subgraph2mat(ind, A, h, node_information)
%  Usage: to extract the enclosing subgraph for a link 
%         Aij (i = ind(1), j = ind(2)), up to h hop 
%
%  *author: Muhan Zhang, Washington University in St. Louis

if nargin < 4
    node_information_flag = 0;  % whether to include node information
else
    node_information_flag = 1;
end

dist = 0;
fringe = [ind];
links = [ind];
nodes = [ind(1); ind(2)];
nodes_dist = [0; 0];
for dist = 1: h
    fringe = neighbors(fringe, A);
    fringe = setdiff(fringe, links, 'rows');
    if isempty(fringe)  % no more new nodes
        subgraph = A(nodes, nodes);
        subgraph(1, 2) = 0;  % ensure subgraph patterns do not contain information about link existence
        subgraph(2, 1) = 0;
        break
    end
    new_nodes = setdiff(fringe(:), nodes, 'rows');
    nodes = [nodes; new_nodes];
    nodes_dist = [nodes_dist; ones(length(new_nodes), 1) * dist];
    links = [links; fringe];
    if dist == h  % nodes enough (reach h hops), extract subgraph
        subgraph = A(nodes, nodes);  % the unweighted subgraph
        subgraph(1, 2) = 0;  % ensure subgraph patterns do not contain information about link existence
        subgraph(2, 1) = 0;
        break
    end
end

sample = {};
sample.nl = {};

% calculate node labels
%labels = nodes_dist + 1;  % use node distance as node labels
labels = node_label(subgraph, 3, h);  % node labeling method
max_nl_size = max(labels);  % determin the nl size after one-hot embedding
sample.nl.values = uint16(labels);

% whether to include node information (embeddings or attributes)
if node_information_flag == 1
    node_info = node_information(nodes, :);
    sample.nl.values = [single(sample.nl.values), node_info];  % use the embeddings/attributes as (continuous) node labels; use single precision to save space
end

sample.am = uint8(full(subgraph));  % adjacency matrix
al = cellfun(@(x) uint16(find(x)),num2cell(sample.am, 2), 'un', 0);  % change to longer integer format if your graph size > 65535
sample.al = al;  % adjacency list, needed by some graph kernels

end


function labels = node_label(subgraph, method, h)
%  Usage: give integer labels to subgraph nodes based on their roles

K = size(subgraph, 1);

if method == 1 || method == 2
    % calculate initial colors based on geometric mean distance to the link
    [dist_to_1, ~, ~] = graphshortestpath(sparse(subgraph), 1, 'Directed', false);
    [dist_to_2, ~, ~] = graphshortestpath(sparse(subgraph), 2, 'Directed', false);
    dist_to_1(isinf(dist_to_1)) = 2 * K;  % replace inf nodes (unreachable from 1 or 2) by an upperbound dist
    dist_to_2(isinf(dist_to_2)) = 2 * K;
    avg_dist = sqrt(dist_to_1 .* dist_to_2);  % use geometric mean as the average distance to the link
    [~, ~, avg_dist_colors] = unique(avg_dist);  % f mapping to initial colors
end

% switch different node labeling methods
switch method
case 1  % palette_wl with initial colors, break ties by nauty
classes = palette_wl(subgraph, avg_dist_colors);
labels = canon(full(subgraph), classes)';
case 2  % use average_dist_colors
labels = avg_dist_colors;
case 3  % node labeling method of SEAL, double-radius
subgraph_wo1 = subgraph(2:end, 2:end);  % subgraph without node 1
subgraph_wo2 = subgraph([1, 3:K], [1, 3:K]);  % subgraph without node 2
[dist_to_1, ~, ~] = graphshortestpath(sparse(subgraph_wo2), 1, 'Directed', false);
[dist_to_2, ~, ~] = graphshortestpath(sparse(subgraph_wo1), 1, 'Directed', false);
dist_to_1(1) = [];
dist_to_2(1) = [];
d = dist_to_1 + dist_to_2;
d_over_2 = floor(d / 2);
d_mod_2 = mod(d, 2);
labels = 1 + min(dist_to_1, dist_to_2) + d_over_2 .* (d_over_2 + d_mod_2 - 1);
labels(isinf(labels)) = 0;
labels(isnan(labels)) = 0;
labels = [1, 1, labels];
labels = labels + 1;  % convert labels 0~MAX to 1~(MAX+1) for one-hot encoding purpose
assert(nnz(labels == 0) == 0)
labels = labels';
end
end


function N = neighbors(fringe, A);
%  Usage: from A to find the neighbor links of all nodes in fringe
N = [];
for no = 1: size(fringe, 1)
    ind = fringe(no, :);
    i = ind(1);
    j = ind(2);
    [~, ij] = find(A(i, :));
    [ji, ~] = find(A(:, j));
    N = [N; [i * ones(length(ij), 1), ij']; [ji, j * ones(length(ji), 1)]];
    N = unique(N, 'rows', 'stable');  % eliminate repeated ones and keep in order
end
end
