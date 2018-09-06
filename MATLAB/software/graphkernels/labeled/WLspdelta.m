% Compute h-step Weisfeiler-Lehman shortest path delta kernel for a set of graphs
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2010 Nino Shervashidze
% Input: Graphs - a 1xN array of graphs
% 		  Graphs(i).am is the adjacency matrix of the i'th graph, 
% 		  Graphs(i).al is the adjacency list of the i'th graph, 
%                 Graphs(i).nl.values is a column vector of node
%                 labels for the i'th graph.
%                 Graphs(i).sp (may be missing) is the shortest
%                 path length matrix of the i'th graph
%                 Graphs(i) may have other fields, but they will not be
%                 used here.
%	 h - a natural number: number of iterations of WL
%	 nl - a boolean: 1 if we want to use original node labels, 0 otherwise
%        spprec - a boolean: 1 if shortest paths are precomputed, 0
%        otherwise. Default: 0
% Output: K - a h+1-element cell array of NxN kernel matrices K for
%             each iter = 0,...,h
%         runtime - scalar (total runtime in seconds)

function [K,runtime] = WLspdelta(Graphs,h,nl,spprec)
% THIS SOURCE CODE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND, AND
% ITS AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH (JMLR) AND
% JMLR'S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES,
% INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON 
% INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR 
% USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR JMLR, NOR JMLR'S
% PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR DAMAGES OF ANY KIND
% RESULTING FROM ITS USE. Without limiting the generality of the foregoing, 
% neither the author, nor JMLR, nor JMLR's publishers and distributors,
% warrant that the Source Code will be error-free, will operate without
% interruption, or will meet the needs of the user.

N=size(Graphs,2);
Lists = cell(1,N);
K = cell(1,h+1);
if nargin<4 spprec=0; end
n_nodes=0;
% compute adjacency lists and n_nodes, the total number of nodes in the dataset
for i=1:N
  Lists{i}=Graphs(i).al;
  n_nodes=n_nodes+size(Graphs(i).am,1);
end
% compute/copy shortest path matrices
maxpath = 0;
if ~spprec
  disp('computing shortest paths...');
  for i=1:N
    Ds{i}=floydwarshall(Graphs(i).am);
    aux=max(Ds{i}(~isinf(Ds{i})));
    if aux > maxpath
      maxpath=aux;
    end
  end
else
  for i=1:N
    Ds{i}=Graphs(i).sp;
    aux=max(Ds{i}(~isinf(Ds{i})));
    if aux > maxpath
      maxpath=aux;
    end
  end
end  
t=cputime; % for measuring runtime

%%% INITIALISATION
% initialize the node labels for each graph with their labels or 
% with degrees (for unlabeled graphs).
label_lookup=containers.Map();
label_counter=uint32(1);
% label_lookup is an associative array, which will contain the
% mapping from multiset labels (strings) to short labels (integers)
if nl==1
  for i=1:N
    % the type of labels{i} is uint32, meaning that it can only handle
    % 2^32 labels and compressed labels over all iterations. If
    % more is needed, switching (all occurences of uint32) to
    % uint64 is a possibility
    labels{i}=zeros(size(Graphs(i).nl.values,1),1,'uint32');
    for j=1:length(Graphs(i).nl.values)
      str_label=num2str(Graphs(i).nl.values(j));
      % str_label is the node label of the current node of the
      % current graph converted into a string
      if ~isKey(label_lookup, str_label)
        %str_label
        %label_counter
        label_lookup(str_label)=label_counter;
        labels{i}(j)=label_counter;
        label_counter=label_counter+1;
      else
        labels{i}(j)=label_lookup(str_label);
      end
    end
  end
else
  for i=1:N
    labels{i}=uint32(full(sum(Graphs(i).am,2)));
    for j=1:length(labels{i})
      str_label=num2str(labels{i}(j));
      % str_label is the node label of the current node of the
      % current graph converted into a string
      if ~isKey(label_lookup, str_label)
        label_lookup(str_label)=label_counter;
        labels{i}(j)=label_counter;                
        label_counter=label_counter+1;
      else
        labels{i}(j)=label_lookup(str_label);
      end
    end
  end
end
L=double(label_counter)-1;
clear Graphs;
disp(['Number of original labels: ',num2str(L)]);
disp(['Number of potential shortest path features: ',num2str((maxpath+1)*L*(L+1)/2)]);
sp=sparse((maxpath+1)*L*(L+1)/2,N);
for i=1:N
  labels_aux=repmat(double(labels{i}),1,length(labels{i}));
  a=min(labels_aux, labels_aux');
  b=max(labels_aux, labels_aux');
  I=triu(~(isinf(Ds{i})));
  Ind=Ds{i}(I)*L*(L+1)/2+(a(I)-1).*(2*L+2-a(I))/2+b(I)-a(I)+1;
  minind=min(Ind);
  diff=max(Ind)-minind;
  aux=accumarray(Ind,ones(nnz(I),1),[],[],[],(minind > 5000 || diff > 3000));
  % sparse of full accumarray depending on the range of values in Ind
  % (and based on empirical observations on the speed of accumarray)
  sp(Ind,i)=aux(Ind);
end
sp=sp(sum(sp,2)~=0,:);
K{1}=full(sp'*sp);

%%% MAIN LOOP
iter=1;
new_labels=labels;
while iter<=h
  disp(['iter=',num2str(iter)]);
  % create an empty lookup table
  label_lookup=containers.Map();
  label_counter=uint32(1);
  for i=1:N
    for v=1:length(Lists{i})
      % form a multiset label of the node v of the i'th graph
      % and convert it to a string
      long_label=[labels{i}(v), sort(labels{i}(Lists{i}{v}))'];
      long_label_2bytes=typecast(long_label,'uint16');
      long_label_string=char(long_label_2bytes);
      % if the multiset label has not yet occurred, add it to the
      % lookup table and assign a number to it
      if ~isKey(label_lookup, long_label_string)
        label_lookup(long_label_string)=label_counter;
        new_labels{i}(v)=label_counter;
        label_counter=label_counter+1;
      else
        new_labels{i}(v)=label_lookup(long_label_string);
      end
    end
  end
  L=double(label_counter)-1;
  disp(['Number of compressed labels: ',num2str(L)]);
  disp(['Number of potential shortest path features: ',num2str((maxpath+1)*L*(L+1)/2)]);
  labels=new_labels;
  sp=sparse((maxpath+1)*L*(L+1)/2,N);
  for i=1:N
    labels_aux=repmat(double(labels{i}),1,length(labels{i}));
    a=min(labels_aux, labels_aux');
    b=max(labels_aux, labels_aux');
    I=triu(~(isinf(Ds{i})));
    Ind=Ds{i}(I)*L*(L+1)/2+(a(I)-1).*(2*L+2-a(I))/2+b(I)-a(I)+1;
    minind=min(Ind);
    diff=max(Ind)-minind;
    aux=accumarray(Ind,ones(nnz(I),1),[],[],[],(minind > 5000 || diff > 3000));
    % sparse of full accumarray depending on the range of values in Ind
    % (and based on empirical observations on the speed of accumarray)
    sp(Ind,i)=aux(Ind);
  end
  sp=sp(sum(sp,2)~=0,:);
  K{iter+1}=K{iter}+full(sp'*sp);
  iter=iter+1;
end
runtime=cputime-t; % computation time of K
end

function [D] = floydwarshall(A, sym, w)
% Input: A - nxn adjacency matrix,
%           sym - boolean, 1 if A and w symmetric
%	    w - nxn weight matrix
% Output: D - nxn distance matrix

n = size(A,1); % number of nodes
D=zeros(n,n);

if nargin<2 % if the graph is not weighted and we have no information about sym, then 
  sym=1;
  w=A;
end

if nargin<3 % if the graph is not weighted, then
  w=A;
end

D=w.*A;
D(A+diag(repmat(Inf,n,1))==0)=Inf; 
D=full(D.*(ones(n)-eye(n))); % set the diagonal to zero

if sym % then it is a bit faster
  for k=1:n
    Daux=repmat(full(D(:,k)),1,n);
    Sumdist=Daux+Daux';
    D(Sumdist<D)=Sumdist(Sumdist<D);
  end
else  
  for k=1:n
    Daux1=repmat(full(D(:,k)),1,n);
    Daux2=repmat(full(D(k,:)),n,1);
    Sumdist=Daux1+Daux2;
    D(Sumdist<D)=Sumdist(Sumdist<D);
  end
end
end