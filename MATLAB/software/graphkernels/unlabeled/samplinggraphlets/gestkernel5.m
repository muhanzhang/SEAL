function [K, runtime] = gestkernel5(graph,samplesize)
% 5-node graphlet kernel
  
% calculate results for pairs of graphs
  
  P = Permutationmatrix5;
  
  tt=cputime; % for measuring runtime
  for i = 1:size(graph,2)
    %  %i
    %  am1 = graph(i).am;
    %  sp1 = floydwarshall(am1);
    %  sp1(find(sp1 == Inf)) = 0;
    %  graph(i).sp = sp1;
    
    m = size(graph(i).am,1);
    if m < 5
      graph(i).dv = zeros(1024,1);
      
    else
      
      
      am = graph(i).am;
      mmax = max(size(am,1),size(am,2));
      
      am(mmax,mmax)=0;
      
      graph(i).dv = zeros(1024,1);
      
      if (samplesize ~= -1)
        
        for j = 1:samplesize
          r= randperm(m+0);
          graph(i).dv(graphlettype(am(r(1:5),r(1:5)))) = ...
              graph(i).dv(graphlettype(am(r(1:5),r(1:5)))) + 1;
          
        end
        
        graph(i).dv = graph(i).dv .* (1/samplesize);
        
      else 
        for r = 1:m
          for s = r+1:m
            for t = s+1:m
              for u = t+1:m
                for v = u+1:m	 
                  graph(i).dv(graphlettype(am([r s t u v],[r s t u v]))) = ...
                      graph(i).dv(graphlettype(am([r s t u v],[r s t u v]))) + 1;
                end
              end 
            end
          end
        end
        
        graph(i).dv = graph(i).dv .* (120 / (m*(m-1)*(m-2)*(m-3)*(m-4)));
        
      end
    end
  end
  
  
  for i = 1:size(graph,2)
    for j = 1:size(graph,2)
      
      K(i,j) = graph(i).dv' * P * graph(j).dv; 
      
    end
  end
  runtime=cputime-tt; % computation time of K
  
end

function result = graphlettype(am)
% determine graphlet type
  
  factor = 2 .^ [0:9]';
  
  upper = [am(1,2:5)';am(2,3:5)'; am(3,4:5)';am(4,5)'];  
  result = sum(factor .* upper) +1;
end
  