function [K,runtime] = allkernel(Graphs,k)

% Compute all k-node graphlet kernel for a set of graphs
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
% Input: Graphs - a 1xn array of graphs
%        k - the size of considered graphlets - 3, 4
% Output: K - nxn kernel matrix K
%         runtime - scalar

n=size(Graphs,2);
switch k
 case 3
  freq=zeros(n,4);
 case 4
  freq=zeros(n,11);
 case 5
  freq=zeros(n,34);
end

t=cputime; % for measuring runtime
switch k
  case 3
   for i=1:n
     freq(i,:)=countall3graphlets(Graphs(i).al);
     sum_freq=sum(freq(i,:));
     if sum_freq~=0 freq(i,:)=freq(i,:)/sum_freq; end 	
   end
  case 4
   for i=1:n
     freq(i,:)=countall4graphlets(Graphs(i).al);
     sum_freq=sum(freq(i,:));
     if sum_freq~=0 freq(i,:)=freq(i,:)/sum_freq; end 	
   end
end
K=freq*freq';

runtime=cputime-t; % computation time of K

end
