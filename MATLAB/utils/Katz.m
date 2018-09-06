function [ auc, precision, sim ] = Katz( train, test, lambda )
%% 计算katz指标并返回AUC值
    sim = inv( sparse(eye(size(train,1))) - lambda * train);   
    % 相似性矩阵的计算
    sim = sim - sparse(eye(size(train,1)));
    [auc, precision] = CalcAUC(train,test,sim);   
    % 评测，计算该指标对应的AUC
end
