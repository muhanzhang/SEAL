function [ auc, precision, sim ] = AA( train, test )
%% 计算AA指标并返回AUC值
    train1 = train ./ repmat(log(sum(train,2)),[1,size(train,1)]); 
    % 计算每个节点的权重，1/log(k_i),网络规模过大时需要分块处理
    train1(isnan(train1)) = 0; 
    train1(isinf(train1)) = 0;  
    % 将除数为0得到的异常值置为0
    sim = train * train1;   clear train1;  
    % 实现相似度矩阵的计算
    [auc, precision] = CalcAUC(train,test,sim);  
    % 评测，计算该指标对应的AUC
end
