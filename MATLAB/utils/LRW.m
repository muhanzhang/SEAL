function [  thisauc ] = LRW( train, test, steps, lambda )
%% 计算LRW指标并返回AUC值
    deg = repmat(sum(train,2),[1,size(train,2)]);
    train = max(train ./ deg,0); clear deg;                                
    % 求转移矩阵
    I = sparse(eye(size(train,1)));                                 
    % 生成单位矩阵
    sim = I;
    stepi = 0;
    while(stepi < steps)                                     
    % 随机游走的迭代
        sim = (1-lambda)*I + lambda * train' * sim;
        stepi = stepi + 1;
    end 
    sim = sim+sim';                                 
    % 相似度矩阵计算完成
    train = spones(train); 
    % 将邻接矩阵还原，因为无孤立点，所以不会有节点的度为0
    thisauc = CalcAUC(train,test,sim, 10000);     
    % 评测，计算该指标对应的AUC
end
