function [K, runtime] = gestkernel3(graph,samplesize)
% 3-node graphlet kernel

% calculate results for pairs of graphs

P = Permutationmatrix3;

t=cputime; % for measuring runtime
for i = 1:size(graph,2)
  %  %i
  %  am1 = graph(i).am;
  %  sp1 = floydwarshall(am1);
  %  sp1(find(sp1 == Inf)) = 0;
  %  graph(i).sp = sp1;
  
  m = size(graph(i).am,1);
  if m < 3
    graph(i).dv = zeros(8,1);
    
  else
    
    
    am = graph(i).am;
    mmax = max(size(am,1),size(am,2));
    
    am(mmax,mmax)=0;
    
    graph(i).dv = zeros(8,1);
    
    if (samplesize ~= -1)
      
      for j = 1:samplesize
        r= randperm(m+0);
        graph(i).dv(graphlettype(am(r(1:3),r(1:3)))) = ...
            graph(i).dv(graphlettype(am(r(1:3),r(1:3)))) + 1;
        
      end
      
      graph(i).dv = graph(i).dv .* (1/samplesize);
      
    else 
      for r = 1:m
        for s = r+1:m
          for t = s+1:m
            graph(i).dv(graphlettype(am([r s t],[r s t]))) = ...
                graph(i).dv(graphlettype(am([r s t],[r s t]))) + 1;
          end
        end
      end
      
      graph(i).dv = graph(i).dv .* (6 / (m*(m-1)*(m-2)));
      
    end
  end
end


for i = 1:size(graph,2)
  %  i
  for j = 1:size(graph,2)
    
    %eye(64,64);
    
    K(i,j) = graph(i).dv' * P * graph(j).dv; 
    
  end
end
runtime=cputime-t; % computation time of K

end

function result = graphlettype(am)
% determine graphlet type
  
  factor = 2 .^ [0:2]';
  
  upper = [am(1,2:3)';am(2,3)'];  
  result = sum(factor .* upper) +1;
end
