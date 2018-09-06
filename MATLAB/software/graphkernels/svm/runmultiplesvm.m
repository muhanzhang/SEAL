function res = runmultiplesvm(Ks,lk,n)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% Input: Ks - cell array of h  m x m kernel matrices
%        lk - m x 1 array of class labels
%        n - number of times we want to run svm
% Output: res is a 1 x n+1 array of structures. Each of the first n
%         elements contains fields optkernel, optc, accuracy, mean_acc and std_acc,
%         and res(n+1) has only two fields - the mean and the std of the n
%         mean_acc's
%
s=0;
for i=1:n
  res(i)=runsvm(Ks,lk);
  meanacc(i) = res(i).mean_acc;
end
res(n+1).mean_acc = mean(meanacc);
res(n+1).std_acc = std(meanacc);
end