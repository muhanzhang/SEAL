function [ thisauc ] = TSRWR( train, test, lambda )
%利用Random walk with restart算法计算的相似度矩阵，做相似性转移
    % 计算RWR相似度矩阵
    train = train + train';
    deg = repmat(sum(train,2),[1,size(train,2)]);
    train = train ./ deg; clear deg;                                % 求转移矩阵
    I = sparse(eye(size(train,1)));                                 % 生成单位矩阵
    sim = (1 - lambda) * inv(I- lambda * train') * I;
    sim = sim+sim';                                                 % 相似度矩阵计算完成
    train = spones(train);
    % 计算相似性转移矩阵
    sim = inv(I - lambda*sim) * sim;
    sim = triu(sim,1);                                              % 取上三角相似性矩阵，并将训练集中存在边的对应元素置为0
    sim = sim - sim.*train;                                         % 只保留测试集和不存在集合中的边的相似度（自环除外）
    thisauc = CalcAUC(train,test,sim);                              % 评测，计算该指标对应的AUC
end
