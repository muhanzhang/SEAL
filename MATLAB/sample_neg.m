function [train_pos, train_neg, test_pos, test_neg] = sample_neg(train, test, k, portion, evaluate_on_all_unseen)
%  Usage: to sample negative links for train and test datasets. 
%         When sampling negative train links, assume all testing
%         links are known and thus sample negative train links 
%         only from other unknown links. Set evaluate_on_all_unseen
%         to true to do link prediction on all links unseen during 
%         training.
%  --Input--
%  -train: half train positive adjacency matrix
%  -test: half test positive adjacency matrix
%  -k: how many times of negative links (w.r.t. pos links) to 
%      sample
%  -portion: if specified, only a portion of the sampled train
%            and test links be returned
%  -evaluate_on_all_unseen: if true, will not randomly sample
%                          negative testing links, but regard
%                          all links unseen during training as 
%                          neg testing links; train negative links 
%                          are sampled in the original way
%  --Output--
%  column indices for four datasets
%%

if nargin < 3
    k = 1;
end

if nargin < 4
    portion = 1;
end

if nargin < 5
    evaluate_on_all_unseen = false;
end
    
n = size(train, 1);
[i, j] = find(train);
train_pos = [i, j];
train_size = length(i);
[i, j] = find(test);
test_pos = [i, j];
test_size = length(i);

if isempty(test)
    net = train;
else
    net = train + test;
end
assert(max(max(net)) == 1);  % ensure positive train, test not overlap
neg_net = triu(-(net - 1), 1);
[i, j] = find(neg_net);
neg_links = [i, j];

% sample negative links
if evaluate_on_all_unseen
    test_neg = neg_links;  % first let all unknown links be negative test links
    % randomly select train neg from all unknown links
    perm = randperm(size(neg_links, 1));
    train_neg = neg_links(perm(1: k * train_size), :);
    test_neg(perm(1: k * train_size), :) = [];  % remove train negative links from test negative links
else
    nlinks = size(neg_links, 1);
    ind = randperm(nlinks);
    if k * (train_size + test_size) <= nlinks
        train_ind = ind(1: k * train_size);
        test_ind = ind(k * train_size + 1: k * train_size + k * test_size);
    else  % if negative links not enough, divide them proportionally
        ratio = train_size / (train_size + test_size);
        train_ind = ind(1: floor(ratio * nlinks));
        test_ind = ind(floor(ratio * nlinks) + 1: end);
    end
    train_neg = neg_links(train_ind, :);
    test_neg = neg_links(test_ind, :);
end

if portion < 1  % only sample a portion of train and test links (for fitting into memory)
    train_pos = train_pos(1:ceil(size(train_pos, 1) * portion), :);
    train_neg = train_neg(1:ceil(size(train_neg, 1) * portion), :);
    test_pos = test_pos(1:ceil(size(test_pos, 1) * portion), :);
    test_neg = test_neg(1:ceil(size(test_neg, 1) * portion), :);
elseif portion > 1  % portion is an integer, number of selections
    train_pos = train_pos(1:portion, :);
    train_neg = train_neg(1:portion, :);
    test_pos = test_pos(1:portion, :);
    test_neg = test_neg(1:portion, :);
end


