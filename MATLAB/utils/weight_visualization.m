function weight_visualization(weights)
%  Usage: to visualize first layer's weights learned from link subgraph patterns, 
%         visualize the adjacency matrix, and sample graphs from it
%  --Input--
%  -weights: a vector of weights to visualize
%  -n: the number of edges in the sampled graph
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%

% if no weights are input, default load weight file in tempdata/nDNN_1
global current_dataset;

if nargin < 1  
weights_file = 'tempdata/nDNN_1/l1_weights.asc';
weights = load(weights_file);  % there are 32 fully connected weight vectors

[rr, cc] = size(weights);
gr = 4;  % how many rows when displaying all 32 weights
gc = rr / 4;

K = (1 + sqrt(8 * cc + 1)) / 2;  % recover the number of vertices

As = [];
for i = 1: rr
    A = zeros(K, K);  % construct adjacency matrix
    A(triu(logical(ones(K, K)), 1)) = weights(i, :);
    A = A + A';
    As = [As, A];
end
As = mat2cell(As, K, K * ones(1, rr));
As = reshape(As, gc, gr)';
As = cell2mat(As);

imagesc(As);
colormap(flipud(gray));
return;
end


% if weights are input (only one weight vector)
[~, cc] = size(weights);
K = (1 + sqrt(8 * cc + 1)) / 2;  % recover the number of vertices

A = zeros(K, K);  % construct adjacency matrix
%weights = abs(weights);
A(triu(logical(ones(K, K)), 1)) = weights;
Ac = A + A';
%Ac = Ac - min(min(Ac));
Ac = Ac - diag(diag(Ac));

figure('position', [100 100 500 500]);
imagesc(Ac);
colormap(flipud(gray));
set(gca,'visible', 'off');


% construct the probability distribution to sample graphs from
for runthis = 1: 0
P = weights;
P(1) = 0;  % do not sample the target link
%P = P.^10;
samples = datasample(1: cc, n, 'Weights', P, 'Replace', false);
weights_sampled = zeros(1, cc);
weights_sampled(samples) = 1;
G = zeros(K, K);  % construct adjacency matrix
G(triu(logical(ones(K, K)), 1)) = weights_sampled;
G = G + G';
g = graph(G);
g = addedge(g, 1, 2, 1);
figure;
gfig = plot(g, '-or', 'Layout', 'circle');
%gfig = plot(g, '-or');
set(gfig, 'LineWidth', 1, 'MarkerSize', 6);
ghighlight = graph(1, 2, 1, K); 
highlight(gfig, [1, 2], 'NodeColor', 'b', 'Marker', 'd', 'MarkerSize', 8);
highlight(gfig, 1, 2, 'EdgeColor', 'b', 'LineStyle', '--', 'LineWidth', 1);
end



