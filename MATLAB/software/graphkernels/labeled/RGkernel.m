function [K,runtime] = RGkernel(Graphs,k)
% Compute height k Ramon & Gaertner kernel for a set of graphs
% Copyright 2012 Nino Shervashidze
% For the algorithm see P. Mahe and J.-P. Vert, Graph kernels based
% on tree patterns for molecules, Arxiv.org, 2008, section 4.2 
% 
% Input: Graphs - a 1xN array of graphs
% 		  Graphs(i).al is the adjacency list of the i'th
% 		  graph, represented as a cell array of line vectors
%                 containing ordered integers.
%                 Graphs(i).nl.values is a column vector of node
%                 labels for the i'th graph.
%	 k - a natural number: maximal tree height
% Output: K - NxN kernel matrix K
%         runtime - scalar

N=size(Graphs,2);
Lists = cell(1,N);
% extract adjacency lists and node labels from the Graphs structure
for i=1:N
  Lists{i}=Graphs(i).al;
  labels{i}=Graphs(i).nl.values;
end
clear Graphs;

t=cputime; % for measuring runtime

K=zeros(N,N);

for i=1:N
  for j=i:N
    li=length(Lists{i}); % li is the number of nodes in the i'th graph
    lj=length(Lists{j});
    k_nodes=double(repmat(labels{i},1,size(labels{j},1))==...
                  repmat(labels{j}',size(labels{i},1),1));
    % k_nodes is one matrix per a pair of graphs (i,j) and
    % k_nodes(u,v) is a node kernel value between the node u of the i'th
    % graph and the node v of the j'th graph    
    K(i,j)=sum(sum(k_nodes));
    for h=2:k
      K(i,j)=0;
      k_nodes_new=zeros(length(Lists{i}), length(Lists{i}));
      % k_nodes_new will replace the former k_nodes in the next iteration
      for u=1:li
        for v=1:lj
          kuv=0;  % disp(['u,v=',num2str([u,v])]); kuv us the node
                  % kernel value for this iteration for the pair (u,v)
          if labels{i}(u)==labels{j}(v)
            if isempty(Lists{i}{u}) || isempty(Lists{j}{v}) 
              M=[]; 
            else
              M=repmat(labels{i}(Lists{i}{u}'), 1, size(labels{j}(Lists{j}{v}'),1))==...
                repmat(labels{j}(Lists{j}{v}')', size(labels{i}(Lists{i}{u}'),1),1);
            end
            [ind1,ind2,junk]=find(M); % nonzero elements = matching pairs
            if size(ind1,2)>1 ind1=ind1'; ind2=ind2'; end
            pairs=[(Lists{i}{u}(ind1'))',(Lists{j}{v}(ind2'))'];
            clear ind1 ind2 junk
            % Here we compute the maximum number of nodes in each
            % graph that can be matched with each other (that is,
            % maximum cardinality that any R can have - and this
            % maximum will be reached).
            commonlabels=intersect(labels{i}(Lists{i}{u}')', labels{j}(Lists{j}{v}')');
            max_matching_size=0;
            for l=1:length(commonlabels)
              max_matching_size=max_matching_size+min(length(find(labels{i}(Lists{i}{u}')==commonlabels(l))),...
                                                      length(find(labels{j}(Lists{j}{v}')==commonlabels(l))));
            end
            clear commonlabels
            if ~isempty(pairs)
              for l=1:size(pairs,1)
                matchings(l).pairs=l;
                kuv=kuv+k_nodes(pairs(l,1), pairs(l,2)); % R's of size 1
              end
            end
            
            % m will denote the cardinality of R's considered in
            % each iteration. Note that if isempty(pairs), max_matching_size=0.
            for m=2:max_matching_size              
              % matchings holds sets R of size m-1
              % new_matchings holds sets R of size m
              new_matchings=[]; new_matchings_counter=1;
              for l=1:length(matchings)
                % n (below) is the number (order) of the last pair in the l'th
                % matching (R) of cardinality m-1
                % Note that matchings are ordered in ascending
                % order to guarantee that we do not encounter the
                % same matching twice. For example the matching 
                % (pair No 1, pair No 3, pair No 7) can exist, but
                % (pair No 1, pair No 7, pair No 3) not. 
                n=matchings(l).pairs(end);
                while n<size(pairs,1)
                  n=n+1;
                  if (any(pairs(matchings(l).pairs,1)==pairs(n,1)) ||...
                      any(pairs(matchings(l).pairs,2)==pairs(n,2)))	
                    continue; % avoid that one node occurs in two pairs
                  end									
                  new_matchings(new_matchings_counter).pairs=[matchings(l).pairs; n];
                  prod=1;
                  for p=1:length(new_matchings(new_matchings_counter).pairs)
                    prod=prod*k_nodes(pairs(new_matchings(new_matchings_counter).pairs(p),1),...
                                     pairs(new_matchings(new_matchings_counter).pairs(p),2));
                  end
                  kuv=kuv+prod;
                  new_matchings_counter=new_matchings_counter+1;
                end
              end
              matchings=new_matchings; % at each iteration the size of matchings grows
            end
          end
          k_nodes_new(u,v)=kuv;
          K(i,j)=K(i,j)+kuv;
        end
      end
      k_nodes=k_nodes_new;
    end
    %k_nodes
    K(j,i)=K(i,j);    
    disp(['kernel value for the pair of graphs ',num2str([i,j]),' computed']);
    % K will hold the kernel
  end
end
runtime=cputime-t;
end
