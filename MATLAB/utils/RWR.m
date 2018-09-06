function [  auc, precision, sim ] = RWR( train, test, lambda )
%% 计算RWR指标并返回AUC值
    deg = repmat(sum(train,2),[1,size(train,2)]);
    train = max(train ./ deg,0); 	clear deg;
    % 求转移矩阵
    I = sparse(eye(size(train,1)));                                
    % 生成单位矩阵
    sim = (1 - lambda) * inv(I- lambda * train') * I;
    %sim = (1 - lambda) ./(I- lambda * train') * I;
    sim = sim+sim';                           
    % 相似度矩阵计算完成
    train = spones(train);   
    % 将邻接矩阵还原，因为无孤立点，所以不会有节点的度为0
    [auc, precision] = CalcAUC(train,test,sim);      
    % 评测，计算该指标对应的AUC
end
