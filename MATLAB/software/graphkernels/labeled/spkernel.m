function [K, runtime, Phi] = spkernel(Graphs, features)
% Compute shortest path kernel for a set of node-labeled graphs
% Copyright 2012 Nino Shervashidze
% Input: Graphs - a 1xN array of graphs
% Output: K - NxN kernel matrix K
%         runtime - scalar
%         features - a boolean: 1 if we want to output the feature
%                    vector representation for each graph, 0 otherwise


N=size(Graphs,2);
Ds = cell(1,N); % shortest path distance matrices for each graph
if nargin<3 features=0; end

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
disp(['the preprocessing step took ', num2str(cputime-t), ' sec']);
t=cputime;

% compute Ds and the length of the maximal shortest path over the dataset
maxpath=0;
for i=1:N
  Ds{i}=floydwarshall(Graphs(i).am);
  aux=max(Ds{i}(~isinf(Ds{i})));
  if aux > maxpath
    maxpath=aux;
  end
 % if rem(i,100)==0 disp(i); end
end

sp=sparse((maxpath+1)*L*(L+1)/2,N);
for i=1:N
  labels_aux=repmat(Graphs(i).nl.values,1,length(Graphs(i).nl.values));
  a=min(labels_aux, labels_aux');
  b=max(labels_aux, labels_aux');
  I=triu(~(isinf(Ds{i})));
  Ind=Ds{i}(I)*L*(L+1)/2+(a(I)-1).*(2*L+2-a(I))/2+b(I)-a(I)+1;
  aux=accumarray(Ind,ones(nnz(I),1));
  sp(Ind,i)=aux(Ind);
end
sp=sp(sum(sp,2)~=0,:);
K=full(sp'*sp);
if features Phi=sp; end
runtime=cputime-t;
disp(['kernel computation took ', num2str(cputime-t), ' sec']);
end

