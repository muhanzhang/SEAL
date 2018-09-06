function [data, label] = graph2vector(pos, neg, A, K)
%  Usage: to convert links' enclosing subgraphs (both pos 
%         and neg) into real vectors, used by WLNM
%  --Input--
%  -pos: indices of positive links
%  -neg: indices of negative links
%  -A: the observed graph's adjacency matrix from which to
%      to extract subgraph features
%  -K: the number of nodes in each link's subgraph
%  --Output--
%  -data: the constructed training data, each row is a 
%         link's vector representation
%  -label: a column vector of links' labels
%
%  *author: Muhan Zhang, Washington University in St. Louis

all = [pos; neg];
pos_size = size(pos, 1);
neg_size = size(neg, 1);
all_size = pos_size + neg_size;

% Generate labels
label = [ones(pos_size, 1); zeros(neg_size, 1)];

% Generate vector data
d = K * (K - 1) / 2;  % dim of data vectors
data = zeros(all_size, d);
one_tenth = floor(all_size / 10);
display('Subgraph Pattern Encoding Begins...')
tic;
%poolobj = parpool(feature('numcores'));  % uncomment this line and change for to parfor in next line to enable parallel computing
for i = 1: all_size
    ind = all(i, :);
    sample = subgraph2vector(ind, A, K);
    data(i, :) = sample;
    %i/all_size  % display finer progress
    progress = i / one_tenth;
    if ismember(progress, [1:10])
        display(sprintf('Subgraph Pattern Encoding Progress %d0%%...', progress));
    end
end
if exist('poolobj')
    delete(poolobj)
end
toc;
end


function sample = subgraph2vector(ind, A, K)
%  Usage: 1) to extract the enclosing subgraph for a link 
%            Aij (i = ind(1), j = ind(2)) 
%         2) to impose a vertex ordering for the vertices
%            of the enclosing subgraph using graph labeling
%         3) to construct an adjacency matrix and output
%            the reshaped vector
%
%  *author: Muhan Zhang, Washington University in St. Louis

D = K * (K - 1) / 2;  % the length of output vector

% Extract a subgraph of K nodes
links = [ind];
links_dist = [0];  % the graph distance to the initial link
dist = 0;
fringe = [ind];
nodes = [ind(1); ind(2)];
nodes_dist = [0; 0];
while 1
    dist = dist + 1;
    fringe = neighbors(fringe, A);
    fringe = setdiff(fringe, links, 'rows');
    if isempty(fringe)  % no more new neighbors, add dummy nodes
        subgraph = A(nodes, nodes);
        subgraph(1, 2) = 0;  % ensure subgraph patterns do not contain information about link existence
        subgraph(2, 1) = 0;
        break
    end
    new_nodes = setdiff(fringe(:), nodes, 'rows');
    nodes = [nodes; new_nodes];
    nodes_dist = [nodes_dist; ones(length(new_nodes), 1) * dist];
    links = [links; fringe];
    links_dist = [links_dist; ones(size(fringe, 1), 1) * dist];
    if size(nodes, 1) >= K  % nodes enough, extract subgraph
        subgraph = A(nodes, nodes);  % the unweighted subgraph
        subgraph(1, 2) = 0;  % ensure subgraph patterns do not contain information about link existence
        subgraph(2, 1) = 0;
        % delete extra links, prepare link-weighted subgraph
        Lia1 = ismember(links(:, 1), nodes);
        Lia2 = ismember(links(:, 2), nodes);
        Lia = logical(Lia1 .* Lia2);
        assert(nnz(Lia==0)==0);
        links = links(Lia, :);
        links_dist = links_dist(Lia, :);
        break
    end
end

% Calculate the link-weighted subgraph, each entry in the adjacency matrix is weighted by the inverse of its distance to the target link
links_ind = sub2ind(size(A), links(:, 1), links(:, 2));
A_copy = A / (dist + 1);  % if a link between two existing nodes < dist+1, it must be in 'links'. The only links not in 'links' are the dist+1 links between some farthest nodes in 'nodes', so here we weight them by dist+1
A_copy(links_ind) = 1 ./ links_dist;
A_copy_u = max(triu(A_copy, 1), tril(A_copy, -1)');  % for links (i, j) and (j, i), keep the smallest dist
A_copy = A_copy_u + A_copy_u';
lweight_subgraph  = A_copy(nodes, nodes);

% Calculate the graph labeling of the subgraph
order = g_label(subgraph);
if length(order) > K  % if size > K, delete the last size-K vertices and reorder
    order(K + 1: end) = [];
    subgraph = subgraph(order, order);
    lweight_subgraph = lweight_subgraph(order, order);
    order = g_label(subgraph);
end

% Generate enclosing subgraph's vector representation
ng2v = 2;  % method for transforming a g_labeled subgraph to vector 
switch ng2v
    case 1  % the simplest way -- one dimensional vector by ravelling adjacency matrix
        psubgraph = subgraph(order, order);  % g_labeled subgraph
        sample = psubgraph(triu(logical(ones(size(subgraph))), 1));
        sample(1) = eps;
    case 2  % use link distance-weighted adjcency matrix, performanc is better
        plweight_subgraph = lweight_subgraph(order, order);  % g_labeled link-weighted subgraph
        sample = plweight_subgraph(triu(logical(ones(size(subgraph))), 1));
        sample(1) = eps;  % avoid inf, and more important, avoid empty vector in libsvm format (directly deleting sample(1) results in libsvm format error)
end
if length(sample) < D  % add dummy nodes if not enough nodes extracted in subgraph
    sample = [sample; zeros(D - length(sample), 1)];
end
end


function N = neighbors(fringe, A);
%  Usage: find the neighbor links of all links in fringe from A

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


function order = g_label(subgraph, p_mo)
%  Usage: impose a vertex order for a enclosing subgraph using graph labeling

if nargin < 2
    p_mo = 7;  % default palette_wl
end

K = size(subgraph, 1);  % local variable

% calculate initial colors based on geometric mean distance to the link
[dist_to_1, ~, ~] = graphshortestpath(sparse(subgraph), 1, 'Directed', false);
[dist_to_2, ~, ~] = graphshortestpath(sparse(subgraph), 2, 'Directed', false);
dist_to_1(isinf(dist_to_1)) = 2 * K;  % replace inf nodes (unreachable from 1 or 2) by an upperbound dist
dist_to_2(isinf(dist_to_2)) = 2 * K;
avg_dist = sqrt(dist_to_1 .* dist_to_2);  % use geometric mean as the average distance to the link
[~, ~, avg_dist_colors] = unique(avg_dist);  % f mapping to initial colors

% switch different graph labeling methods
switch p_mo
case 1  % use classical wl, no initial colors
classes = wl_string_lexico(subgraph);
order = canon(full(subgraph), classes)';
case 2  % use wl_hashing, no initial colors
classes = wl_hashing(subgraph);
order = canon(full(subgraph), classes)';
case 3  % use classical wl, with initial colors
classes = wl_string_lexico(subgraph, avg_dist_colors);
order = canon(full(subgraph), classes)';
case 4  % use wl_hashing, with initial colors
classes = wl_hashing(subgraph, avg_dist_colors);
order = canon(full(subgraph), classes)';
case 5  % directly use nauty for canonical labeling
order = canon(full(subgraph), ones(K, 1))';
case 6  % no graph labeling, directly use the predefined order
order = [1: 1: K];
case 7  % palette_wl with initial colors, break ties by nauty
classes = palette_wl(subgraph, avg_dist_colors);
%classes = palette_wl(subgraph);  % no initial colors
order = canon(full(subgraph), classes)';
case 8  % random labeling
order = randperm(K);
end
end
