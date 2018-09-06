function [  auc, precision, sim ] = SimRank( train, test, lambda)
%% 计算SimRank指标并返回AUC值
    deg = sum(train,1);     
    % 求节点的入度，构成行向量，供调用
    lastsim = sparse(size(train,1), size(train,2)); 
    % 存储前一步的迭代结果，初始化为全0矩阵
    sim = sparse(eye(size(train,1))); 
    ntrain = train.*repmat(max(1./deg,0),size(train,1),1);
    % approximate SimRank
    for iter = 1:5
        sim = max(lambda*(ntrain'*sim*ntrain),eye(size(train,1)));
    end
    [auc, precision] = CalcAUC(train,test,sim);    
end
