function [K,runtime]=untilpRWkernel(Graphs, p, small, intermediate, name)
% Compute a labeled p-step random walk kernel for a set of graphs
% Copyright 2012 Nino Shervashidze
% Input: Graphs - a 1xN array of graphs
% 	 Graphs(i).am is the adjacency matrix, Graphs(i).nl.values is a column vector of node labels.
%        p - integer scalar: the length up to which we consider all walks
%        small - a boolean, indicating how "small" is the dataset.
%                If graphs are roughly less than 100 nodes big and
%                there are not more than 10-20 labels, then set small=1, 0
%                otherwise. In the latter case the filtered adjacency matrices 
%                will be computed for each pair of graphs, so the computation
%                will be slower, but more memory-efficient.
%        intermediate (MAKES SENSE ONLY IF small=1) - a boolean: 1
%                     if we want to output kernel matrices
%                     corresponding to k=1..p, 0 otherwise (default: 0)
%        name (MAKES SENSE ONLY IF small=1) - a string: the name of
%                     files where the kernel matrices will be
%                     exported to (with suffixes _1,...,_p) (default: [])
% Output: K - NxN kernel matrix K
%         runtime - scalar: runtime in seconds

K=[]; runtime=0;
if nargin < 3 disp('The function requires 3 arguments'); return; end
if nargin < 4 intermediate=0; end
if nargin < 5 name=''; intermediate=0; end % if we do not know where
                                           % to put intermediate
                                           % matrices, we do not
                                           % output them

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

L2=L^2;
if small
  for k=1:p
    disp(['p = ',num2str(k)]);
    for l=1:L2
   %   disp(['l = ', num2str(l)]);
      for i=1:N
        x=sum((G(i).filteredam(l).am)^k,2);
        for j=i:N
%          disp(['processing pair ',num2str(i),',',num2str(j)]);
          y=sum((G(j).filteredam(l).am)^k,1);
          K(i,j)=K(i,j)+sum(sum(x*y));
          K(j,i)=K(i,j);
        end
      end
    end
    if intermediate
      runtime=cputime-t;
      save([name,'_', num2str(k)], 'K','runtime');
    end
  end
else
  for i=1:N
    g1=G(i);
    counter=1;
    for k=1:L
      for l=1:L
        aux = double(g1.nl.values==labelset(k))*double(g1.nl.values==labelset(l))';
        g1.filteredam(counter).am=sparse(g1.am.*aux);
        counter=counter+1;
      end
    end
    for j=i:N
%     disp(['processing pair ',num2str(i),',',num2str(j)]);
      g2=G(j);
      counter=1;
      for k=1:L
        for l=1:L
          aux = double(g2.nl.values==labelset(k))*double(g2.nl.values==labelset(l))';
          g2.filteredam(counter).am=sparse(g2.am.*aux);
          counter=counter+1;
        end
      end
      for l=1:L2
        for k=1:p
%         disp(['p = ',num2str(k)]);
          x=sum((g1.filteredam(l).am)^k,2);
          y=sum((g2.filteredam(l).am)^k,1);
          K(i,j)=K(i,j)+sum(sum(x*y));
          K(j,i)=K(i,j);
        end
      end
    end
  end
end  
    
runtime=cputime-t;
disp(['kernel computation took ', num2str(cputime-t), ' sec']);
end
