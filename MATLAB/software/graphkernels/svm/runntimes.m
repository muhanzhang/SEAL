function result = runntimes(K,lk,n)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% Input: K  - m x m kernel matrix
%        lk - m x 1 array of class labels
%        n - number of times we want to run svm
% Output: result - a structure with fields accuracy, mean, std and
%                  mean, std and se are the mean, the standard
%                  deviation and the standard error of accuracy

rand('seed',666 );
accuracy = zeros(n,1);

for i = 1:n
[junk1, accuracy(i), junk2] = runIndependent(K, lk)

end

result.mean = mean(accuracy)

result.std = std(accuracy)

result.se = result.std / sqrt(n)

result.accuracy = accuracy

result
