function [ thisauc ] = Sorenson( train, test )
%% 计算Sorenson指标并返回AUC值
    sim = train * train;                                            
    % 计算分子
    sim = triu(sim,1);
    deg_col = repmat(sum(train,2), [1 size(train,1)]);              
    % 计算分母
    deg_col = triu(deg_col' + deg_col);
    sim = 2 * sim ./ deg_col;                             
    % 相似度矩阵计算完成
    sim(isnan(sim)) = 0; sim(isinf(sim)) = 0;
    thisauc = CalcAUC(train,test,sim, 10000);       
    % 评测，计算该指标对应的AUC
end
