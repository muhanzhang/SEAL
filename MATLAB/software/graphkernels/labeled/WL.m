% Compute h-step Weisfeiler-Lehman kernel for a set of graphs
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2009 Nino Shervashidze
% Input: Graphs - a 1xN array of graphs
% 		  Graphs(i).am is the adjacency matrix of the i'th graph, 
% 		  Graphs(i).al is the adjacency list of the i'th graph, 
%                 Graphs(i).nl.values is a column vector of node
%                 labels for the i'th graph.
%                 Graphs(i) may have other fields, but they will not be
%                 used here.
%	 h - a natural number: number of iterations of WL
%	 nl - a boolean: 1 if we want to use original node labels, 0 otherwise
% Output: K - a h+1-element cell array of NxN kernel matrices K for
%             each iter = 0,...,h
%         runtime - scalar (total runtime in seconds)

function [K,runtime] = WL(Graphs,h,nl)
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
n_nodes=0;
% compute adjacency lists and n_nodes, the total number of nodes in the dataset
for i=1:N
  Lists{i}=Graphs(i).al;
  n_nodes=n_nodes+size(Graphs(i).am,1);
end
phi=sparse(n_nodes,N); %each column j of phi will be the explicit feature 
% representation for the graph j

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
      phi(labels{i}(j),i)=phi(labels{i}(j),i)+1;
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
      phi(labels{i}(j),i)=phi(labels{i}(j),i)+1;
    end
  end
end
L=label_counter-1;
disp(['Number of original labels: ',num2str(L)]);
clear Graphs;
K{1}=full(phi'*phi);

%%% MAIN LOOP
iter=1;
new_labels=labels;
while iter<=h
  disp(['iter=',num2str(iter)]);
  % create an empty lookup table
  label_lookup=containers.Map();
  label_counter=uint32(1);
  % create a sparse matrix for feature representations of graphs
  phi=sparse(n_nodes,N);
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
    % fill the column for i'th graph in phi
    aux=accumarray(new_labels{i}, ones(length(new_labels{i}),1));
    phi(new_labels{i},i)=phi(new_labels{i},i)+aux(new_labels{i});
  end
  L=label_counter-1;
  disp(['Number of compressed labels: ',num2str(L)]);
  K{iter+1}=K{iter}+full(phi'*phi);
  labels=new_labels;
  iter=iter+1;
end
runtime=cputime-t; % computation time of K
end

