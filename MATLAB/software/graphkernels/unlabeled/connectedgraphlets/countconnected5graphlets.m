function [count] = countconnected5graphlets(A, L)

% Count all 5-node connected subgraphs in an undirected graph
% without node labels and with unweighted edges
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
%
% Input: A - nxn adjacency matrix
% 	 L - 1xn cell array - corresponding adjacency list
% Output: count - 1x21 vector of integers (21 being the number of
%                 non-isomorphic connected 5-node graphs). count(1) is the number of occurences of the
%		  graphlet with 5 nodes and 10 edges (type 1), count(2) is the number of 
%		  occurences of the graphlet with 5 nodes and 9 edges (type 2) and so
%		  forth (see g_1--g_21 in ../5graphlets.pdf)

w = [1/120, 1/72, 1/48, 1/36, 1/28, 1/20, 1/14, 1/10, 1/12, 1/8, 1/8, 1/4, 1/2, 1/12, 1/12, 1/4, 1/4, 1/2, 0,0,0];
% 1/number of length-4 paths in graphlets of type 1-21 respectively
n = size(A,1); % number of nodes
count = zeros(1,21);

for i=1:n
  for j=L{i}
    for k=L{j}
      if k~=i
        for l=L{k}
          if l~=i && l~=j
            for m=L{l}
              if m~=i && m~=j && m~=k
         
                %%% Tests: Which type of connected 5 node graphlet
                %%% is induced by nodes in this length-4 simple
                %%% path (that is, nodes i,j,k,l,m)?
                aux = A(i,k)+A(i,l)+A(i,m)+A(j,l)+A(j,m)+A(k,m);
                if aux==6 % if a graphlet has 10 edges, there is
                          % only one possibility - it is the
                          % complete graph of 5 nodes 
                  count(1)=count(1)+w(1); % count increases by
                                          % w(1), because this
                                          % graphlet will be
                                          % counted as many times,
                                          % as the number of 5-node
                                          % paths it contains
                else
                  if aux==5 % if it has 9 edges, there is only one
                            % possibility as well
                    count(2)=count(2)+w(2);
                  else
                    if aux==4 % if it has 8 edges, it can be either
                              % graphlet 3 or 4, which can be
                              % distinguished by looking at the
                              % minimum degree of the graphlet
                      aux1 = min(sum(A([i,j,k,l,m],[i,j,k,l,m]),1));
                      if aux1==2 % the sorted degree distribution of the graphlet of type 4 is (2 3 3 4 4)
                        count(4)=count(4)+w(4);
                      else % the sorted degree distribution of the graphlet of type 3 is (3 3 3 3 4)
                        count(3)=count(3)+w(3);
                      end
                    else 
                      if aux==3 % if the graphlet has 7 edges,
                                % it can be of type 5, 6, 9, or 14
                        aux1 = sort(sum(A([i,j,k,l,m],[i,j,k,l,m]),1));
                        if aux1(1)==1 % the sorted degree distribution of the graphlet of type 9 is (1 3 3 3 4)
                          count(9)=count(9)+w(9);
                        else
                          if aux1(2)==3 % the sorted degree distribution of the graphlet of type 5 is (2 3 3 3 3) 
                            count(5)=count(5)+w(5);
                          else
                            if aux1(3)==2 % the sorted degree distribution of the graphlet of type 14 is (2 2 2 4 4)
                              count(14)=count(14)+w(14);
                            else % the sorted degree distribution of the graphlet of type 6 is (2 2 3 3 4)
                              count(6)=count(6)+w(6);	
                            end
                          end
                        end
                      else 
                        if aux==2 % if the graphlet has 6 edges, if
                                  % can be of type 7, 10, 11, 15,
                                  % or 16
                          aux2=sum(A([i,j,k,l,m],[i,j,k,l,m]),1);
                          aux1 = sort(aux2);
                          if aux1(1)==1
                            if aux1(3)==2 % the sorted degree distribution of the graphlet of type 16 is (1 2 2 3 4)
                              count(16)=count(16)+w(16);
                            else % the sorted degree distribution of the graphlet of type 10 is (1 2 3 3 3)
                              count(10)=count(10)+w(10);			
                            end
                          else 
                            if aux1(4)==2 % the sorted degree distribution of the graphlet of type 11 is (2 2 2 2 4)
                              count(11)=count(11)+w(11);
                            else 
                              % there are two types of connected length-6 graphlets left: 7 and 15. both have 
                              % the same sorted degree distribution - (2 2 2 3 3), but in 7 the nodes with degree 3
                              % are connected with an edge and in 15 they are not. 
                              ind=find(aux2==3);
                              path=[i,j,k,l,m];										
                              if A(path(ind(1)), path(ind(2)))==1
                                count(7)=count(7)+w(7);
                              else
                                count(15)=count(15)+w(15);
                              end
                              clear path;
                            end	
                          end
                        else 
                          if aux==1 % if the graphlet has 5 edges,
                                    % it can be of type 8, 12, 17,
                                    % or 18
                            aux2 = sum(A([i,j,k,l,m],[i,j,k,l,m]),1);		
                            aux1 = sort(aux2);
                            if aux1(1)==2 % the sorted degree distribution of the graphlet of type 8 is (2 2 2 2 2) 
                              count(8)=count(8)+w(8);
                            else
                              if aux1(2)==1 % the sorted degree distribution of the graphlet of type 18 is (1 1 2 3 3) 
                                count(18)=count(18)+w(18);
                              else
                                % there are two types of connected length-5 graphlets (containing length-4 paths) left: 
                                % 12 and 17. Both have the same sorted degree distribution - (1 2 2 2 3), but in 17 
                                % the nodes with degrees 1 and 3 are connected with an edge, while in 12 they are not.
                                ind=[find(aux2==3), find(aux2==1)];
                                path=[i,j,k,l,m];
                                if A(path(ind(1)), path(ind(2)))==1
                                  count(17)=count(17)+w(17);
                                else
                                  count(12)=count(12)+w(12);
                                end			
                                clear path;
                              end
                            end
                          else % if the graphlet has 4 edges, there
                               % is only one possibility 
                            count(13)=count(13)+w(13);
                          end
                        end
                      end
                    end
                  end
                end
                %%% end Tests
              end
            end
          end
        end
      end
    end
  end
  
  %%%% count graphlets of type 20 %%%%
  for j=L{i}
    for k=L{j}
      if k~=i && A(i,k)==0
        for l=L{k}
          if l~=i && l~=j && A(i,l)==0 && A(j,l)==0
            for m=L{k}
              if m~=i && m~=j && m~=l && A(i,m)==0 && A(j,m)==0 && A(l,m)==0
                count(20)=count(20)+0.5;
              end
            end
          end
        end
      end
    end
  end
  %%%% end count graphlets of type 20 %%%%
  
  %%%% count graphlets of type 19 and 21 %%%%
  N=length(L{i});
  for j=1:N-3
    for k=j+1:N-2
      for l=k+1:N-1
        for m=l+1:N % with this kind of enumeration we make sure
                    % that i, j, k, l, and m are pairwise different
          aux = A(L{i}(j),L{i}(k)) + A(L{i}(j),L{i}(l)) + ...
                A(L{i}(j),L{i}(m)) + A(L{i}(k),L{i}(l)) + ...
                A(L{i}(k),L{i}(m)) + A(L{i}(l),L{i}(m));
          if aux == 1
            count(19)=count(19)+1;
          else
            if aux == 0 %%% "star"
              count(21)=count(21)+1;
            end
          end
        end
      end
    end
  end
  %%%% end count graphlets of type 19 and 21 %%%%	
end

end
