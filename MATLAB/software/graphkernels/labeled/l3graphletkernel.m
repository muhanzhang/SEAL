function [K,runtime] = l3graphletkernel(Graphs)
% Compute a 3-node connected labeled graphlet kernel on a set of graphs
% Author: Nino Shervashidze- nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze
%
% Input: Graphs - a 1xN array of graphs
% Output: K - NxN kernel matrix K
%         runtime - scalar
 

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
B=2*L^3; % upper bound on the size of the set of possible 
       % connected 3-node graphlets labeled with this alphabet
disp(['the preprocessing step took ', num2str(cputime-t), ' sec']);
t=cputime;

phi=sparse(B,N); %each column j of phi will be the explicit feature 
		 % representation for the graph j
graphlet_lookup=containers.Map();
graphlet_counter=1;

for g=1:N
  % just some shorthands...
  A=Graphs(g).am;
  nl=Graphs(g).nl.values;
  Al=Graphs(g).al;
  for i=1:length(Al)
    for j=Al{i}
      for k=Al{j}
        if i~=k
          if A(i,k)==0
            graphlet=[1,  min(nl(i),nl(k)),  nl(j),   max(nl(i),nl(k))];
            increment=1/2;
          else 
            graphlet=[2,  sort([nl(i),nl(j),nl(k)])];
            increment=1/6;
          end
          graphlet_2_bytes=typecast(graphlet,'uint16');
          graphletstring=char(graphlet_2_bytes);
          if ~isKey(graphlet_lookup, graphletstring)
            graphlet_lookup(graphletstring)=graphlet_counter;
            phi(graphlet_counter,g)=phi(graphlet_counter,g)+increment;
            graphlet_counter=graphlet_counter+1;
          else
            phi(graphlet_lookup(graphletstring),g)= ...
                phi(graphlet_lookup(graphletstring),g)+increment;
          end
        end
      end
    end
  end
end

K=full(phi'*phi);
runtime=cputime-t; % computation time of K
disp(['kernel computation took ', num2str(cputime-t), ' sec']);
end

