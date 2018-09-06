function [result,mean_accuracy,std_accuracy] = runIndependent(K,lk)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% K = kernel matrix (n*n)
% lk = vector of labels (n*1)
% cv = number of folds in cross-validation

% standard deviation
% independent scheme
% best c

addpath('~/code/libsvm');
% randomly permute kernel matrix and labels
r = randperm(size(K,1));
K = K(r,r);
lk = lk(r);
lk=lk';
lkoriginal = lk; 
Koriginal = K;

%% stratified cross-validation
%sum(sum(Koriginal));
%neworder = stratifiedsplit(lk)
%for i = 1:size(neworder,2) 
%m = size(neworder(i).old,1);
%r = randperm(m);
%newlk(neworder(i).new) =  lk(neworder(i).old(r));
%Knew([neworder(i).new]',[neworder(i).new]') = K([neworder(i).old(r)]',[neworder(i).old(r)]');  
%end
%
%sum(sum(Knew)) - sum(sum(Koriginal))
%dbstop
%lk = newlk'
%K = Knew;
%size(lk);
%size(K);
%dbstop 
%Koriginal = K;


% bring kernel matrix into libsvm format
p80 = ceil(size(K,2) * 0.8);
p90 = ceil(size(K,2) * 0.9);

% specify range of c-values
cvalues = (10 .^ [-7:2:7]) / size(K,2);

cv = 10;
fs = size(K,2) - p90;

% cross-validation loop
for k = 1:cv
K = Koriginal;
lk = lkoriginal;

K = K([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs],[k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]);  
lk = lk([k*fs+1:size(K,2),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 
K = makepos(K);
K1 = [(1:size(K,1))', normalizekm(K)];

%if any(strcmp('optimal',options))
imresult=[];
for i = 1:size(cvalues,2)
    % train on 80%, predict on 10% (from 81% to 90%) 
  size(lk(1:p80));
  size(K1(1:p80,1:p80+1));
  model = svmtrain(lk(1:p80,1), K1(1:p80,1:p80+1), strcat(['-t 4  -c ' num2str(cvalues(i))]));
  [predict_label, accuracy, dec_values] = svmpredict(lk(p80+1:p90,1),K1(p80+1:p90,1:p80+1), model);
  accuracy80 = accuracy;
  imresult(i)= accuracy(1);
end
  
  
  % determine optimal c
  [junk,optimalc]= max(fliplr(imresult));
  optimalc = size(cvalues,2)+1 - optimalc; 
  % train on 90% with optimal c, predict on 10% (from 91% to 100%)
  model = svmtrain(lk(1:p90,1), K1(1:p90,1:p90+1),strcat(['-t 4  -c ' num2str(cvalues(optimalc))]) );
  [predict_label, accuracy, dec_values] = svmpredict(lk(p90+1:size(K,1),1), K1(p90+1:size(K,1),1:p90+1), model);
  accuracy90 = accuracy
  result(k)=accuracy(1)


end  
mean_accuracy =  mean(result) 
std_accuracy = std(result)


end

%
%% cross-validation
%if any(strcmp('cv',options))
%options = strcat(['-t 4 -v ' num2str(cv) ' -c ' num2str(cvalues(i))])
%result(i) = svmtrain(lk, K1, options); %', num2str(cv)));
%end
%
%end

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
