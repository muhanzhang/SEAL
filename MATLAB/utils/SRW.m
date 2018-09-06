function [  thisauc ] = SRW( train, test, steps, lambda )
%% 计算SRW指标并返回AUC值
    deg = repmat(sum(train,2),[1,size(train,2)]);
    train = max(train ./ deg,0); clear deg;
    % 求转移矩阵
    I = sparse(eye(size(train,1)));                                 
    % 生成单位矩阵
    tempsim = I;                            
    % 用来暂存每步的迭代结果
    stepi = 0; sim = sparse(size(train,1),size(train,2));           
    % 随机游走的迭代 sim用来存储每步迭代的分值之和
    while(stepi < steps)
        tempsim = (1-lambda)*I + lambda * train' * tempsim;
        stepi = stepi + 1;
        sim = sim + tempsim;
    end
    sim = sim+sim';                        
    % 相似度矩阵计算完成
    train = spones(train);   
    %将邻接矩阵还原，因为无孤立点，所以不会有节点的度为0
    thisauc = CalcAUC(train,test,sim, 10000);    
    % 评测，计算该指标对应的AUC
end
