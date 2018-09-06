function [K,runtime]=lRWkernel(Graphs, lambda, small)
% Compute a random walk kernel for a set of labeled graphs
% Copyright 2012 Karsten Borgwardt, Nino Shervashidze
% Input: Graphs - a 1xN array of graphs
% 	 Graphs(i).am is the adjacency matrix, Graphs(i).nl.values is a column vector of node labels.
%        lambda - scalar: a rule of thumb for setting it is to
%                 take the largest power of 10 which is smaller than 1/d^2,
%                 d being the largest degree in the dataset
%        small - a boolean, indicating how "small" the dataset is.
%                If graphs are roughly less than 100 nodes big and
%                there are not more than 10-20 labels, then set small=1, 0
%                otherwise. In the latter case the filtered adjacency matrices 
%                will be computed for each pair of graphs, so the computation
%                will be slower, but more memory-efficient.
% Output: K - NxN kernel matrix K
%         runtime - scalar: runtime in seconds

K=[]; runtime=0;
if nargin < 3 disp('The function requires 3 arguments'); return; end

N=size(Graphs,2);

t=cputime; % for measuring runtime

%%% PREPROCESSING (mainly in order to find out the size of the node
%%% label alphabet, L, and rename labels as 1 ,..., L)
label_lookup=containers.Map();
label_counter=1;
for i=1:N
  for j=1:length(Graphs(i).nl.values)
    str_label=num2str(Graphs(i).nl.values(j));
    % str_label is the node label of the current node of the
    % current graph converted into a string
    if ~isKey(label_lookup, str_label)
      label_lookup(str_label)=label_counter;
      Graphs(i).nl.values(j)=label_counter;
      label_counter=label_counter+1;
    else
      Graphs(i).nl.values(j)=label_lookup(str_label);
    end
  end
end
L=label_counter-1; % L is the size of the node label alphabet
labelset=[1:L];
disp(['the preprocessing step took ', num2str(cputime-t), ' sec']);
t=cputime;
K=zeros(N,N);
for i=1:N
  %	if rem((i-1),100)==0 
  %		disp(['filtering graph ',num2str(i)]); 
  %	end 
  G(i).nl.values=Graphs(i).nl.values;
  G(i).am=Graphs(i).am;
  if small
    counter=1;
    for k=1:L
      for l=1:L
        aux = double(Graphs(i).nl.values==labelset(k))*double(Graphs(i).nl.values==labelset(l))';
        G(i).filteredam(counter).am=sparse(Graphs(i).am.*aux);
        counter=counter+1;
      end
    end
  end
  
end
clear Graphs;
for i=1:N
  for j=i:N
    %disp(['processing pair ',num2str(i),',',num2str(j)]);
    K(i,j)=K(i,j)+labeledrandomwalk(G(i),G(j),labelset,lambda, small);
    K(j,i)=K(i,j);
  end
end
runtime=cputime-t;
disp(['kernel computation took ', num2str(cputime-t), ' sec']);
end

function result = labeledrandomwalk(g1,g2,labelset,lambda, small)
L=size(labelset,2);
if ~small
  %g1.am is the adjacency matrix of graph g1
  %g2.am is the adjacency matrix of graph g2
  counter=1;
  for k=1:L
    for l=1:L
      aux = double(g1.nl.values==labelset(k))*double(g1.nl.values==labelset(l))';
      g1.filteredam(counter).am=sparse(g1.am.*aux);
      counter=counter+1;
    end
  end
  counter=1;
  for k=1:L
    for l=1:L
      aux = double(g2.nl.values==labelset(k))*double(g2.nl.values==labelset(l))';
      g2.filteredam(counter).am=sparse(g2.am.*aux);
      counter=counter+1;
    end
  end
end

am1 = g1.filteredam; % filtered version of g1, where am1(1).am is the first
                     % filtered am, am1(2).am is the second filtered am, etc

am2 = g2.filteredam; % filtered version of g2, where am2(1).am is the first
                     % filtered am, am2(2).am is the second filtered am, etc


%lambda = 10e-2;

[x,rubbish] = pcg(@(x)smtfilter(x,am1,am2,lambda),ones(size(g1.am,1)*size(g2.am,1),1),1e-6,20);
result =  sum(sum(x));

end

function result = smtfilter(x,am1,am2,lambda)

yy = 0;
for i = 1:size(am1,2)
  yy = yy + vec(am1(i).am * invvec(x,size(am1(i).am,1),size(am2(i).am,1)) * am2(i).am);
end

yy = lambda * yy;

vecu = minus(x,yy);

result = vecu;
end

function v = vec(M)
[m,n] = size(M);
v = reshape(M,m*n,1);
end

function v = invvec(M,m,n)
v = reshape(M,m,n);

end
