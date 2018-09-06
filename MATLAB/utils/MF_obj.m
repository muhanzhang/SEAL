function obj = MF_obj(train,IDX,gm,x)
b = x(:,end);
x(:,end) = [];
pred = x*x'+bsxfun(@plus,b,b')+gm;
obj = (train(IDX)-pred(IDX));
obj = obj'*obj;