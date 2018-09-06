function [  thisauc ] = CosPlus( train, test )
%% 计算Cos+指标并返回AUC值
    D = sparse(eye(size(train,1)));                        
    % 生成稀疏的单位矩阵
    D(logical(D)) = sum(train,2);  
    % 生成度矩阵 （对角线元素为同下标节点的度）
    pinvL = sparse(pinv( full(D - train) ));      clear D;
    % 拉普拉斯矩阵的伪逆  
    Lxx = diag(pinvL);   
    % 取对角线元素
    tmp = Lxx*Lxx';
    tmp(tmp<0)=0;
    sim = pinvL ./ (tmp).^0.5; clear tmp;                         
    % 求相似度矩阵
    sim(isnan(sim)) = 0; sim(isinf(sim)) = 0;
    thisauc = CalcAUC(train,test,sim, 10000);      
    % 评测，计算该指标对应的AUC
end
