% Random graph construction routine.
% Note 1: Default is Erdos-Renyi graph G(n,0.5)
% Note 2: Generates undirected, simple graphs only
%
% INPUTS:  N - number of nodes
%          p - probability, 0<=p<=1
%          E - fixed number of edges; if specified, p is irrelevant
% OUTPUTS: adj - adjacency matrix of generated graph (symmetric), nxn
%
% Other routines used: numEdges.m
% GB: last updated, Oct 20, 2012

function adj = randomGraph(n,p,E)

adj=zeros(n); % initialize adjacency matrix

switch nargin   % number of function arguments
  
 case 1  % just the number of nodes, n
  
  % 0.5 - default probability of attachment
  for i=1:n
    for j=i+1:n
      if rand<=0.5; adj(i,j)=1; adj(j,i)=1; end
    end
  end
 
 case 2 % the number of nodes and the probability of attachment, n, p
  
  for i=1:n
    for j=i+1:n
      if rand<=p; adj(i,j)=1; adj(j,i)=1; end
    end
  end
  
 case 3 % fixed number of nodes and edges, n, E
  
  while numEdges(adj) < E
    i=randi(n); j=randi(n);  % pick two random nodes
    if i==j || adj(i,j)>0; continue; end  % do not allow self-loops or double edges
    adj(i,j)=adj(i,j)+1; adj(j,i)=adj(i,j);
  end
    
 
end
net = adj;
save(sprintf('data/%srandom.mat', num2str(p, 2)), 'net');




function m = numEdges(adj)

sl=selfLoops(adj); % counting the number of self-loops

if isSymmetric(adj) && sl==0    % undirected simple graph
    m=sum(sum(adj))/2;
    
elseif isSymmetric(adj) && sl>0
    m=(sum(sum(adj))-sl)/2+sl; % counting the self-loops only once

elseif not(isSymmetric(adj))   % directed graph (not necessarily simple)
    m=sum(sum(adj));
    
end