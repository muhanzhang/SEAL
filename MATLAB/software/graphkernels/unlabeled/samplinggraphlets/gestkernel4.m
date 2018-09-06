function     [K,      runtime]     =     gestkernel4(graph,samplesize)
% graphlet kernel
% @author Karsten Borgwardt
% @date May 11th 2007
  
% calculate results for pairs of graphs
  
  P = Permutationmatrix4;
  
  t=cputime; % for measuring runtime
  for i = 1:size(graph,2)
    %  %i
    %  am1 = graph(i).am;
    %  sp1 = floydwarshall(am1);
    %  sp1(find(sp1 == Inf)) = 0;
    %  graph(i).sp = sp1;
    
    m = size(graph(i).am,1);
    if m < 4
      graph(i).dv = zeros(64,1);
      
    else
      
      
      am = graph(i).am;
      mmax = max(size(am,1),size(am,2));
      
      am(mmax,mmax)=0;
      
      graph(i).dv = zeros(64,1);
      
      if (samplesize ~= -1)
        
        
        
        for j = 1:samplesize
          r= randperm(m+0);
          graph(i).dv(graphlettype(am(r(1:4),r(1:4)))) = ...
              graph(i).dv(graphlettype(am(r(1:4),r(1:4)))) + 1;
          
        end
        
        graph(i).dv = graph(i).dv .* (1/samplesize);
        
      else 
        for r = 1:m
          for s = r+1:m
            for t = s+1:m
              for u = t+1:m
                
                graph(i).dv(graphlettype(am([r s t u],[r s t u]))) = ...
                    graph(i).dv(graphlettype(am([r s t u],[r s t u]))) + 1;
              end 
            end
          end
        end
        %(graph(i).dv'*P)'
        graph(i).dv = graph(i).dv .* (24 / (m*(m-1)*(m-2)*(m-3)));
        
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
  
  factor = 2 .^ [0:5]';
  
  upper = [am(1,2:4)';am(2,3:4)'; am(3,4)'];  
  result = sum(factor .* upper) +1;
end
