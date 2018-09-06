function f_patterns = frequent_patterns(data, n, pos)
%  Usage: to extract the n most frequent patterns from data (one row 
%         represents one sample)
% *author: Muhan Zhang, Washington University in St. Louis
global current_dataset;

if nargin < 2
    n = 1;
end

[~, cc] = size(data);
K = (1 + sqrt(8 * cc + 1)) / 2;  % recover the number of vertices

[udata, ~, ids] = unique(data, 'rows');  % extract unique patterns

f_patterns = [];
for i = 1: n
    id = mode(ids);  % select most frequent
    %id = randi(max(ids));  % random select
    f_patterns = [f_patterns; udata(id, :)];
    ids(ids == id) = [];
end

% draw frequent patterns
green = [0 204 102]/255;
orange = [255 128 0]/255;
red = [255 80 80]/255;
blue = [102 178 255]/255;
for i = 1: n
    G = zeros(K, K);  % construct adjacency matrix
    G(triu(logical(ones(K, K)), 1)) = f_patterns(i, :);
    G = G + G';
    g = graph(G);
    g = addedge(g, 1, 2, 1);
    figure('position', [100 100 500 500]);
    gfig = plot(g, '-o', 'EdgeColor', orange, 'NodeColor', orange, 'Layout', 'force');
    %gfig = plot(g, '-or');
    set(gfig, 'LineWidth', 1, 'MarkerSize', 5);
    
    if pos == 1  % to highlight the positive link
        highlight(gfig, [1, 2], 'NodeColor', red, 'Marker', 'd', 'MarkerSize', 6);
        highlight(gfig, 1, 2, 'EdgeColor', red, 'LineStyle', '-', 'LineWidth', 2);
    else
        highlight(gfig, [1, 2], 'NodeColor', red, 'Marker', 'd', 'MarkerSize', 6);
        highlight(gfig, 1, 2, 'EdgeColor', red, 'LineStyle', '--', 'LineWidth', 2);
    end
    set(gca,'visible','off');
    cd figure/;
    %evalc(sprintf('export_fig fpattern_%s_%d -pdf -m2 -transparent -c[20,20,20,20]', current_dataset, i));
    cd ..;
end
