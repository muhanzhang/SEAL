function [K,runtime] = connectedkernel(Graphs, k)
% Compute connected k-node graphlet kernel for a set of graphs
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
% Input: Graphs - a 1xn array of graphs
%        k - the size of considered graphlets - 3, 4, or 5
% Output: K - nxn kernel matrix K
%         runtime - scalar

n=size(Graphs,2);
switch k
 case 3
  freq=zeros(n,2);
 case 4
  freq=zeros(n,6);
 case 5
  freq=zeros(n,21);
end

t=cputime; % for measuring runtime
switch k
  case 3
   for i=1:n
     freq(i,:)=countconnected3graphlets(Graphs(i).am, Graphs(i).al);
     sum_freq=sum(freq(i,:));
     if sum_freq~=0 freq(i,:)=freq(i,:)/sum_freq; end % for relative number of graphlets
   end
 case 4
   for i=1:n
     freq(i,:)=countconnected4graphlets(Graphs(i).am, Graphs(i).al);
     sum_freq=sum(freq(i,:));
     if sum_freq~=0 freq(i,:)=freq(i,:)/sum_freq; end % for relative number of graphlets
   end
 case 5
   for i=1:n
     freq(i,:)=countconnected5graphlets(Graphs(i).am, Graphs(i).al);
     sum_freq=sum(freq(i,:));
     if sum_freq~=0 freq(i,:)=freq(i,:)/sum_freq; end % for relative number of graphlets
   end
end
K=freq*freq';

runtime=cputime-t; % computation time of K

end

