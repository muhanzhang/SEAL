function [ thisauc ] = MFI( train,test )
%% 计算MFI指标并返回AUC值
    I = sparse(eye(size(train,1)));                                 
    % 生成单位矩阵
    D = I;
    D(logical(D)) = sum(train,2);         
    % 生成度矩阵（对角线元素为同下标节点的度）
    L = D - train;          
    clear D;              
    % 拉普拉斯矩阵
    sim = inv(I + L);      
    clear I L;       
    % 相似度矩阵的计算
    thisauc = CalcAUC(train,test,sim, 10000);   
    % 评测，计算该指标对应的AUC
end
