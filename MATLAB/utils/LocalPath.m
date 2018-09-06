function [  thisauc ] = LocalPath( train, test, lambda )
%% 计算LP指标并返回AUC值
    sim = train*train;    
    % 二阶路径
    sim = sim + lambda * (train*train*train);   
    % 二阶路径 + 参数×三节路径
    thisauc = CalcAUC(train,test,sim, 10000);  
    % 评测，计算该指标对应的AUC
end
