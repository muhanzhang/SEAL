function [ auc, precision ] = CalcAUC(train, test, sim)
% calculate auc on given testing links for heurisitc-based methods
% test is a struct, test.pos contains the positive testing links
% sim is the matrix of similarity scores
% train is legacy variable (useless)

    test_pos = test.pos;
    test_neg = test.neg;
    test_pos = sub2ind(size(sim), test_pos(:, 1), test_pos(:, 2));
    test_neg = sub2ind(size(sim), test_neg(:, 1), test_neg(:, 2));

    pos_pre = sim(test_pos);
    neg_pre = sim(test_neg);
    
    labels = [ones(1, size(pos_pre, 1)), zeros(1, size(neg_pre, 1))];
    scores = [pos_pre', neg_pre'];
    [~, ~, ~, auc] = perfcurve(labels, scores, 1);
    [xs, ys, ~, aucpr] = perfcurve(labels, scores, 1, 'XCrit', 'reca', 'YCrit', 'prec');
    precision = sum(diff(xs) .* ys(2:end));  % average precision
    
end
