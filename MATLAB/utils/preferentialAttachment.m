% Routine implementing simple preferential attachment for network growth.
% The probability that a new vertex attaches to a given old vertex
%                           is proportional to the (total) vertex degree.
% Note 1: Vertices arrive one at a time.
% Note 2: Assume undirected simple graph.
% Source: Newman, "The Structure and Function of Complex Networks"
%         B-A., "Emergence of Scaling in Random Networks"
%
% INPUTs: n - final (desired) number of vertices, 
%         m - # edges to attach at every step
% OUTPUTs: edge list, [number of edges x 3]
%
% Other routines used: weightedRandomSample.m
% GB: last updated, November 9, 2012

function net = preferentialAttachment(n,m)

vertices = 2;
if not(vertices<=n); printf('Specify more than 2 nodes.\n');  return; end
el=[1 2 1; 2 1 1];      % start with one edge


while vertices < n
  
  vertices=vertices+1;  % add new vertex

  if m>=vertices
    for node=1:vertices-1
      el = [el; node vertices 1];
      el = [el; vertices node 1];  % add symmetric edge
    end
    continue
  end
    
  deg=[];        % compute nodal degrees for this iteration
  for v=1:vertices; deg=[deg; numel(find(el(:,1)==v))]; end
  

  % add m edges
  r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  while not(length(unique(r))==length(r))
    r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  end
  
  for node=1:length(r)
    el = [el; r(node) vertices 1];
    el = [el; vertices r(node) 1];      
  end      
  
end

el = el(el(:, 1) < el(:, 2), :);
net = full(sparse(el(:, 1), el(:, 2), el(:, 3), n, n));
net = net + net';
save(sprintf('data/%dPA.mat', m), 'net');



function s = weightedRandomSample(n,P,W)

s = [];

if abs(sum(W)-1)>10^(-8); 
    printf('The probabilities do not sum up to 1.\n');
    s = 'probabilities do not sum up to 1';
    return; 
end

% divide the unit interval into |P| segments each with length W_i
unit = [0,W(1)];
for w=2:length(W)
  unit = [unit W(w)+unit(length(unit))];
end

% draw a random number in the unit interval uniformly - where does it fall?
while length(s)<n

  lb = max(find(unit<rand));    % rand is in [unit(lb), unit(lb+1)]
  s = [s P(lb)];                % pick P(lb)

end