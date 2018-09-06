function [count] = countall3graphlets(L)

% Count all 3-node subgraphs in an undirected graph
% without node labels and with unweighted edges
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
% Input: L - 1xn cell array - adjacency list
% Output: count - 1x4 vector of integers. count(i)= number of graphlets with 
%                 3-i+1 edges (see 3graphlets.pdf)

n=length(L); % number of nodes
count=zeros(1,4);
w=[1/6, 1/4, 1/2];
for v1=1:n
  for v2=L{v1}
    cardinalities=card_inter(L{v1}, L{v2}, length(L{v1}), length(L{v2}));
    count(1)=count(1)+w(1)*cardinalities(3);
    count(2)=count(2)+w(2)*(cardinalities(1)+cardinalities(2)-2);
    count(3)=count(3)+w(3)*(n-sum(cardinalities));
  end
end
count(4)=n*(n-1)*(n-2)/6-sum(count(1:3));
end

function [n] = card_inter(o_set1, o_set2, l1, l2)
% Find the cardinality of the intersection of two ordered sets of lengths l1 and l2 respectively
% n(1)=o_set1\o_set2, n(2)=o_set2\o_set1, n(3)=(o_set2 inter o_set1)
n=zeros(1,3);
i=1; j=1;

while i<=l1 && j <=l2
  if o_set1(i)<o_set2(j) n(1)=n(1)+1; i=i+1;
  else
    if o_set1(i)>o_set2(j) n(2)=n(2)+1; j=j+1;
    else i=i+1; j=j+1; n(3)=n(3)+1;
    end
  end
end
n(1)=n(1)+l1-i+1;
n(2)=n(2)+l2-j+1;
end
