function [ thisauc ] = TSAA( train, test, lambda )
%利用Common Neighbor算法计算的相似度矩阵，做相似性转移
    % 计算AA相似度矩阵
    train = train + train';
    train1 = train ./ repmat(log(sum(train,2)),[1,size(train,1)]);  % 计算每个节点的权重，1/log(k_i)
                                                                    % 注意：网络规模过大，repmat操作会溢出，此时需要分块处理
    train1(isnan(train1)) = 0; train1(isinf(train1)) = 0;           % 将除数为0得到的异常值置为0
    sim = train * train1;      clear train1;                        % 实现相似度矩阵的计算
    % 计算相似性转移矩阵
    I = sparse(eye(size(train,1)));
    sim = inv(I - lambda*sim) * sim;
    sim = triu(sim,1);                                              % 取上三角相似性矩阵，并将训练集中存在边的对应元素置为0
    sim = sim - sim.*train;                                         % 只保留测试集和不存在集合中的边的相似度（自环除外）
    thisauc = CalcAUC(train,test,sim);                              % 评测，计算该指标对应的AUC
end
