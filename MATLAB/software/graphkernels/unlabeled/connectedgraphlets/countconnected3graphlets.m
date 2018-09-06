function [count] = countconnected3graphlets(A, L)

% Count all 3-node connected subgraphs in an undirected graph
% without node labels and with unweighted edges
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
%
% Input: A - nxn adjacency matrix
% 	 L - 1xn cell array - corresponding adjacency list
% Output: count - 1x2 vector of integers. count(1) is the number of occurences of the
%		  graphlet with 3 nodes and 2 edges (type 1), count(2) is the number of 
%		  occurences of the graphlet with 3 nodes and 3 edges (type 2)

w = [1/2, 1/6]; % 1/number of length-2 paths in graphlets of type 1 and type 2 respectively
n = size(A,1); % number of nodes
count = zeros(1,2);

for i=1:n
  for j=L{i}
    for k=L{j}	
      if k~=i
        if A(i,k)==1
          count(2)=count(2)+w(2);
        else
          count(1)=count(1)+w(1);
        end
      end
    end
  end
end

end

