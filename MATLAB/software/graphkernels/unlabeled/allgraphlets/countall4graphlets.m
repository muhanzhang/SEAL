function [count] = countall4graphlets(L)

% Count all 4-node subgraphs in an undirected graph
% without node labels and with unweighted edges
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze 
% Input:    L - 1xn cell array - corresponding adjacency list
% Output:   count - 1x11 vector of integers. count(1) is the number of occurences of the
%		    graphlet with 4 nodes and 6 edges (type 1), count(2) is the number of 
%		    occurences of the graphlet with 4 nodes and 5 edges (type 2) and so
%		    forth (see 4graphlets.pdf)



w = [1/12, 1/10, 1/8, 1/6, 1/8, 1/6, 1/6, 1/4, 1/4, 1/2, 0];

n = length(L); % number of nodes
count = zeros(1,11);

%precompute the number of edges
m=0;
for i=1:n
  m=m+length(L{i});
end
m=m/2;

for v1=1:n
  for v2=L{v1}
    K=0;
    internalcount=zeros(1,11);
    [inter, diff1, diff2]=card_inter(L{v1},L{v2},length(L{v1}),length(L{v2})); 
    for v3=inter
      cardinalities=card_3inter(L{v1},L{v2},L{v3},length(L{v1}),...
                                length(L{v2}), length(L{v3}));
      internalcount(1)=internalcount(1)+1/2*cardinalities(7);
      internalcount(2)=internalcount(2)+1/2*(cardinalities(4)-1);
      internalcount(2)=internalcount(2)+1/2*(cardinalities(5)-1);
      internalcount(2)=internalcount(2)+1/2*(cardinalities(6)-1);
      internalcount(3)=internalcount(3)+1/2*(cardinalities(1));
      internalcount(3)=internalcount(3)+1/2*(cardinalities(2));
      internalcount(3)=internalcount(3)+cardinalities(3);
      internalcount(7)=internalcount(7)+ (n-sum(cardinalities));
      K=K+1/2*cardinalities(7)+1/2*(cardinalities(5)-1)+1/2*(cardinalities(6)-1)+cardinalities(3);
    end
    for v3=setdiff(diff1,v2)
      cardinalities=card_3inter(L{v1},L{v2},L{v3},length(L{v1}),...
                                length(L{v2}), length(L{v3}));
      internalcount(2)=internalcount(2)+1/2*(cardinalities(7));
      internalcount(3)=internalcount(3)+1/2*(cardinalities(4));
      internalcount(3)=internalcount(3)+1/2*(cardinalities(5));
      internalcount(5)=internalcount(5)+1/2*(cardinalities(6)-1);
      internalcount(4)=internalcount(4)+1/2*(cardinalities(1)-2);
      internalcount(6)=internalcount(6)+1/2*(cardinalities(2));
      internalcount(6)=internalcount(6)+cardinalities(3);
      internalcount(8)=internalcount(8)+(n-sum(cardinalities));
      K=K+1/2*(cardinalities(7))+1/2*(cardinalities(5))+1/2*(cardinalities(6)-1)+cardinalities(3);
    end
    for v3=setdiff(diff2,v1)
      cardinalities=card_3inter(L{v1},L{v2},L{v3},length(L{v1}),...
                                length(L{v2}), length(L{v3}));
      internalcount(2)=internalcount(2)+1/2*(cardinalities(7));
      internalcount(3)=internalcount(3)+1/2*(cardinalities(4));
      internalcount(5)=internalcount(5)+1/2*(cardinalities(5)-1);
      internalcount(3)=internalcount(3)+1/2*(cardinalities(6));
      internalcount(6)=internalcount(6)+1/2*(cardinalities(1));
      internalcount(4)=internalcount(4)+1/2*(cardinalities(2)-2);
      internalcount(6)=internalcount(6)+cardinalities(3);
      internalcount(8)=internalcount(8)+(n-sum(cardinalities));
      K=K+1/2*(cardinalities(7))+1/2*(cardinalities(5)-1)+1/2*(cardinalities(6))+cardinalities(3);
    end
    internalcount(9)=internalcount(9)+(m+1-length(L{v1})-length(L{v2})-K);
    internalcount(10)=internalcount(10)+((n-length(inter)-length(diff1)-length(diff2) )*...
                                         (n-length(inter)-length(diff1)-length(diff2)-1)/2-...
                                         (m+1-length(L{v1})-length(L{v2})-K));
    count=count+internalcount.*w; 
    %ones(size(internalcount))
    % DEBUG    
    %    v1
    %    v2
    %    K
  end
end

count(11)= n*(n-1)*(n-2)*(n-3)/(4*3*2) - sum(count(1:10));
end

function [inter, diff1, diff2] = card_inter(o_set1, o_set2, l1, l2)
% Find the intersection and set differences of two ordered sets of 
% lengths l1 and l2 respectively

inter=zeros(1,min(l1,l2));
diff1=zeros(1,max(l1,l2));
diff2=zeros(1,max(l1,l2));

i=1; j=1;
while i<=l1 && j <=l2
  if o_set1(i)<o_set2(j) 
    diff1(i)=o_set1(i);
    i=i+1;
  else 
    if o_set1(i)>o_set2(j)
      diff2(j)=o_set2(j);
      j=j+1;
    else
      inter(i)=o_set1(i);
      i=i+1; j=j+1;
    end
  end
end
if i<=l1
  diff1(i:l1)=o_set1(i:l1);
else 
  if j<=l2
    diff2(j:l2)=o_set2(j:l2);
  end
end
inter=inter(inter>0);
diff1=diff1(diff1>0);
diff2=diff2(diff2>0);
end

