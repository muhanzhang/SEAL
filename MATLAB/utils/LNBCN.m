function [ thisauc ] = LNBCN( train, test ) 
%% 计算局部朴素贝叶斯模型性CN指标并返回AUC值
    s = size(train,1)*(size(train,1)-1) / nnz(train) -1;  
    % 计算每个网络中的常量s
    tri = diag(train*train*train)/2;     
    % 计算每个点所在的三角形个数
    tri_max = sum(train,2).*(sum(train,2)-1)/2;  
    % 每个点最大可能所在的三角形个数
    R_w = (tri+1)./(tri_max+1); clear tri tri_max; 
    % 接下来几步是按照公式度量每个点的角色  
    SR_w = log(s)+log(R_w); clear s R_w;
    SR_w(isnan(SR_w)) = 0; SR_w(isinf(SR_w)) = 0;
    SR_w = repmat(SR_w,[1,size(train,1)]) .* train;   
    % 节点的角色计算完毕
    sim = spones(train) * SR_w;   clear SR_w;                       
    % 将节点对（x,y）的共同邻居的角色量化值相加即可
    thisauc = CalcAUC(train,test,sim, 10000);
    % 评测，计算该指标对应的AUC
end
