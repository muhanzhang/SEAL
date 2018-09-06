function [ thisauc ] = WAA( train, test,alpha )
% weighted AA
    train = train.^alpha;
    [r,d] = size(train);
    weight = 1./log(1+sum(train,1));
    weight(isnan(weight)) = 0; 
    weight(isinf(weight)) = 0;  
    sim = zeros(r);
    for i=1:r
        tmp = sum(bsxfun(@times, bsxfun(@plus,train(i,:),train) .* spones(bsxfun(@times,train(i,:),train)), weight ),2);
        sim(i,:) = tmp';
    end
    sim = sim - diag(diag(sim));
    thisauc = CalcAUC(train,test,sim, 10000);
end
