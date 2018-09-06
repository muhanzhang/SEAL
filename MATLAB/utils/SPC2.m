function [  thisauc ] = SPC( train, test,k )
%% 计算ACT指标并返回AUC值
    D = sparse(eye(size(train,1)));   
    % 生成稀疏的单位矩阵
    D(logical(D)) = sum(train,2);   
    % 生成度矩阵（对角线元素为同下标节点的度）
%     pinvL = sparse(pinv( full(D - train) ));   clear D;
%     % 拉普拉斯矩阵的伪逆
%     Lxx = diag(pinvL);     
%     % 取对角线元素
%     Lxx = repmat(Lxx, [1,size(train,1)]);   
%     % 将对角线元素向量扩展为n×n阶矩阵
    if nargin < 3 
        k=10;
    end
    invD = 1./D;
    invD(isinf(invD)) = 0;
    invD(isnan(invD)) = 0;
    L = full(eye(size(train,1)) - invD*train);
    [U,S] = eig(L);
    S=diag(S);
    
    %k=1000;
    [ S,I]=sort(S,'ascend');
    U=U(:,I(1:k));
    %U = U(:,S'<0.01);
    
    sim = U*U';
    

%     xii = sum(bsxfun(@times,U,U),2);
%     sim = bsxfun(@plus,xii,xii')-2*U*U';
    
             
    sim(isnan(sim)) = 0; 
    sim(isinf(sim)) = 0;
    sim(sim<0)=0;
    
    %sim(isinf(sim)) = 0;    %convert diagnal infs to 0, no influence
    thisauc = CalcAUC(train,test,sim, 10000);    
    % 评测，计算该指标对应的AUC
end 
