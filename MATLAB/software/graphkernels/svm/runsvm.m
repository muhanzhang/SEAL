function [res] = runsvm(Ks,lk)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% runsvm(Ks,lk)
% K = 1 x h cell array of kernelmatrices (n*n)
% lk = vector of class labels (n*1)
% cv = number of folds in cross-validation

% independent scheme
% best c

addpath('~/code/libsvm');
n=length(lk) % size of the dataset
% randomly permute labels: r will be also used for permuting the kernel matrices
r = randperm(n);
lk = lk(r);

% specify range of c-values
cvalues = (10 .^ [-7:2:7]) / size(lk,1);

cv = 10;
p80 = ceil(n * (1-2/cv));
p90 = ceil(n * (1-1/cv));
fs = n - p90; % fold size


% output variables
res.optkernel=zeros(cv,1);
res.optc=zeros(cv,1);
res.accuracy=zeros(cv,1);

% cross-validation loop
opth=zeros(1,cv);
for k = 1:cv
  imresult=[];
  
  height = length(Ks);
  for h=1:height
    K = Ks{h}(r,r);
    K_current = K([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs],[k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]);  
    lk_current = lk([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 
    K_current = makepos(K_current);
    K1 = [(1:size(K_current,1))', normalizekm(K_current)];
    
    
    for i = 1:size(cvalues,2)
      % train on 80%, predict on 10% (from 81% to 90%) 
      size(lk_current(1:p80));
      size(K1(1:p80,1:p80+1));
      model = svmtrain(lk_current(1:p80,1), K1(1:p80,1:p80+1), strcat(['-t 4  -c ' num2str(cvalues(i))]));
      [predict_label, accuracy, dec_values] = svmpredict(lk_current(p80+1:p90,1),K1(p80+1:p90,1:p80+1), model);
      imresult(h,i)= accuracy(1);
    end
  end
  
  % determine optimal h and c
  [junk,position]= max(imresult(:));
  [optimalh, indoptimalc]=ind2sub(size(imresult),position);
  
  opth(k)=optimalh;
  res.optc(k)= cvalues(indoptimalc);
  res.optkernel(k)=optimalh;
  % train on 90% with optimal c, predict on 10% (from 91% to 100%)
  K = Ks{optimalh}(r,r);
  K_current = K([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs],[k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]);  
  lk_current = lk([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 
  K_current = makepos(K_current);
  K1 = [(1:size(K_current,1))', normalizekm(K_current)];
  
  model = svmtrain(lk_current(1:p90,1), K1(1:p90,1:p90+1),strcat(['-t 4  -c ' num2str(cvalues(indoptimalc))]) );
  [predict_label, accuracy, dec_values] = svmpredict(lk_current(p90+1:size(K,1),1), K1(p90+1:size(K,1),1:p90+1), model);
  res.accuracy(k)=accuracy(1)
end
res.mean_acc =  mean(res.accuracy) 
res.std_acc = std(res.accuracy)
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
