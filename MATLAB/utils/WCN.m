function [ thisauc ] = WCN( train, test,alpha )
% weighted common neghbors
    train = train.^alpha;
    [r,d] = size(train);
    sim = zeros(r);
    for i=1:r
        tmp = sum(bsxfun(@plus,train(i,:),train) .* spones(bsxfun(@times,train(i,:),train)),2);
        sim(i,:) = tmp';
    end
    sim = sim - diag(diag(sim));
    thisauc = CalcAUC(train,test,sim, 10000);
end
