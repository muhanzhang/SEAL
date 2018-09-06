function [ auc, precision, sim ] = CN( train, test )
%% 计算CN指标并返回AUC值
    sim = train * train;        
    % 相似度矩阵的计算
    [auc, precision] = CalcAUC(train,test,sim);
    % 评测，计算该指标对应的AUC
end
