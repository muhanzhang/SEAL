function [  thisauc ] = LHNII( train, test, lambda )   
%% 计算LHN2指标并返回AUC值
    M = nnz(train)/2;
    % 网络中的边数
    D = sparse(eye(size(train,1)));                                 
    D(logical(D)) = sum(train,2);   
    % 生成度矩阵 （对角线元素为同下标节点的度）
    D = inv(D);  
    %D = inv(D+0.1*eye(size(D,1)));  
    % 求度矩阵的逆矩阵
    maxeig = max(eigs(train));  
    % 求邻接矩阵的最大特征值
    tempmatrix = (sparse(eye(size(train,1))) - lambda/maxeig * train); 
    tempmatrix = inv(tempmatrix);
    sim = 2 * M * maxeig * D * tempmatrix * D;   clear D tempmatrix;
    % 完成相似度矩阵的计算
    thisauc = CalcAUC(train,test,sim, 10000);    
    % 评测，计算该指标对应的AUC
end
