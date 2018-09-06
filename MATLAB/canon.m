function [order] = canon(subgraph, classes)
%  Usage: invoke nauty to find the canonical labelling
%  --Input--
%  -subgraph: the adjacency matrix of the graph to find cl
%  -classes: the colors of the subgraph nodes
%  --Output--
%  -order: the position vector of vertices in the new graph
%          e.g., original label = [1 2 3 4 5], canonical 
%          label = [3 1 2 4 5], then order = [2 3 1 4 5]
%
%  *author: Muhan Zhang,  Washington University in St. Louis

K = size(subgraph, 1);
if nargin < 2
    classes = ones(K, 1);
end

% Reorder subgraph to let adjacent vertices have the same colors
% The colors must be like [1, 2, 1, 3, 3], must not be like 
% [1, 2, 1, 4, 4]. Colors must be continuous from 1 to n.
[classes, order] = sort(classes);  % to sort the colors
subgraph1 = subgraph(order, order);

% Prepare the input to canonical.c
classes = [classes; classes(end) + 1];
colors_nauty = 1 - diff(classes);
num_edges = nnz(subgraph1);
degrees = sum(subgraph1, 2);

% Check if canonical.c has been compiled to mex function
flag = exist(['canonical.' mexext]);
%flag = 0;  % let it be compiled every time
if flag == 0
    !rm canonical.mex*;
    cd software/nauty26r7;
    !cp ../../canonical.c .;
    mex canonical.c nauty.c nautil.c naugraph.c schreier.c naurng.c nausparse.c;
    !cp canonical.mex* ../../;
    cd ../..;
end

% Run nauty to find canonical labeling
subgraph1;
colors_nauty;
order;
clabels = canonical(subgraph1, num_edges, degrees, colors_nauty);
clabels = clabels + 1;
order = order(clabels);

