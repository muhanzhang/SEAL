function [ thisauc ] = TSCN( train, test, lambda )
%% 计算TSCN指标并返回AUC值
    sim = train * train;                       
    % 计算共同邻居相似度矩阵
    I = sparse(eye(size(train,1)));
    sim = inv(I - lambda*sim) * sim;
    % 相似度转移
    thisauc = CalcAUC(train,test,sim, 10000);     
    % 评测，计算该指标对应的AUC
end
