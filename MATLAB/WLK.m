function [auc, precision] = WLK(train_mix, test, h, ith_experiment)
%  Weisfeiler-Lehman graph kernel for link prediction
%  For usage, please refer to SEAL.m
%  The code is partly adapted from the graph kernel toolbox
%  of "Shervashidze et. al., Weisfeiler-lehman graph kernels"
% 
%  *author: Muhan Zhang, Washington University in St. Louis
%%

% please download graphkernels and libsvm to software/
addpath(genpath('software/graphkernels'));
addpath(genpath('software/libsvm-3.22/matlab'));

A = train_mix.train;
data_name = train_mix.data_name;
train_pos = train_mix.pos; train_neg = train_mix.neg;
test_pos = test.pos; test_neg = test.neg;

if nargin < 3
    h = 2;
end
if nargin < 4
    ith_experiment = 1;
end

if h == 'auto'
    % automatically select h from {1, 2} by validation performance of AA and CN
    [val_train, val_test] = DivideNet(A, 0.9, false);
    h_val_train = triu(sparse(val_train), 1);
    h_val_test = triu(sparse(val_test), 1);
    [~, ~, val_test_pos, val_test_neg] = sample_neg(h_val_train, h_val_test, 1, 1);
    val_test = {};
    val_test.pos = val_test_pos;
    val_test.neg = val_test_neg;
    val_auc_AA = AA(val_train, val_test)
    val_auc_CN = CN(val_train, val_test)
    if val_auc_AA >= val_auc_CN
        h = 2;
        disp(['Choose h=', num2str(h)])
    else
        h = 1;
        disp(['Choose h=', num2str(h)])
    end
end

% extract enclosing subgraphs
train_size = size(train_pos, 1) + size(train_neg, 1)
test_size = size(test_pos, 1) + size(test_neg, 1)
[data, ~] = graph2mat([train_pos; train_neg], [test_pos; test_neg], A, h, ith_experiment, 1, data_name); 
label = [ones(size(train_pos, 1), 1); zeros(size(train_neg, 1), 1); ...
         ones(size(test_pos, 1), 1); zeros(size(test_neg, 1), 1)];

% permutate the train set
perm = randperm(train_size);
data(:, 1:train_size) = data(:, perm);
label(1:train_size) = label(perm);

% run graph kernel
X = data;  % the graph X
lk = label;  % the labels Y

% compute kernel matrices
[Ks, runtime] = WL(X, 5, 1);  % use the Weisfeiler-Lehman subtree kernel
display('kernel computation time:')
runtime

% run svm
% specify range of c-values
n = length(lk);
cvalues = (10 .^ [-7:2:7]) / n;

% calculate the stopping index of train/val
fs = test_size; % test size
ptrain = ceil((n - fs) * 0.9);  % keep 90% remaining data as training set
pval = n - fs;  % stopping position of validataion set

cv = 1;  % do not use cv, use fixed train/val/test splits
res.optkernel=zeros(cv,1);
res.optc=zeros(cv,1);
res.accuracy=zeros(cv,1);

% cross-validation loop
opth=zeros(1,cv);
for k = 1:cv
  imresult=[];
  height = length(Ks);
  for h=1:height
    K = Ks{h};
    %K_current = K([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs],[k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]);  
    K_current = K;
    %lk_current = lk([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 
    lk_current = lk;
    K_current = makepos(K_current);
    K1 = [(1:size(K_current,1))', normalizekm(K_current)];

    for i = 1:size(cvalues,2)
      % predict on validation set
      size(lk_current(1:ptrain));
      size(K1(1:ptrain,1:ptrain+1));
      model = svmtrain(lk_current(1:ptrain,1), K1(1:ptrain,1:ptrain+1), strcat(['-t 4  -c ' num2str(cvalues(i))]));
      [predict_label, accuracy, dec_values] = svmpredict(lk_current(ptrain+1:pval,1),K1(ptrain+1:pval,1:ptrain+1), model);
      imresult(h,i)= accuracy(1);
    end
  end

  % determine optimal h and c
  [junk,position]= max(imresult(:));
  [optimalh, indoptimalc]=ind2sub(size(imresult),position);

  opth(k)=optimalh;
  res.optc(k)= cvalues(indoptimalc);
  res.optkernel(k)=optimalh;
  % train and predict on testset
  K = Ks{optimalh};
  %K_current = K([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs],[k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]);  
  %lk_current = lk([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 
  K_current = K;
  lk_current = lk;
  K_current = makepos(K_current);
  K1 = [(1:size(K_current,1))', normalizekm(K_current)];

  model = svmtrain(lk_current(1:pval,1), K1(1:pval,1:pval+1),strcat(['-b 1 -t 4 -c ' num2str(cvalues(indoptimalc))]) );
  [predict_label, accuracy, dec_values] = svmpredict(lk_current(pval+1:size(K,1),1), K1(pval+1:size(K,1),1:pval+1), model, '-b 1');
  res.accuracy(k)=accuracy(1);
  if model.Label(1) == 0
    dec_values = dec_values(:,2);
  else
    dec_values = dec_values(:,1);
  end
  [~, ~, ~, auc] = perfcurve(lk_current(pval+1:size(K,1)), dec_values', 1); 
  [xs, ys, ~, aucpr] = perfcurve(lk_current(pval+1:size(K,1)), dec_values', 1, 'XCrit', 'reca', 'YCrit', 'prec');
  precision = sum(diff(xs) .* ys(2:end));  % average precision
  res.auc(k) = auc;
  res.precision(k) = precision;
end
res.mean_acc =  mean(res.accuracy);
res.std_acc = std(res.accuracy);
res.mean_auc = mean(res.auc);
res.std_auc = std(res.auc);
res.mean_precision = mean(res.precision);
res.std_precision = std(res.precision);
auc = res.mean_auc;
precision = res.mean_precision;
res
end


function result = makepos(K)
pd = 0;
addc = 10e-7;
while (pd ==  0)

addc = addc * 10
try
if (isinf(addc) == 1)
pd = 1;
else
chol(normalizekm(K + eye(size(K,1),size(K,1)) * addc));
pd = 1;
end
catch

end

end
if (isinf(addc)==0)
result = K + eye(size(K,1),size(K,1)) * addc;
else
result = eye(size(K,1));
end
end




