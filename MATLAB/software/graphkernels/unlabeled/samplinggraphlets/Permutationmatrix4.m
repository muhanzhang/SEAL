function P = Permutationmatrix4
% graphlet kernel
% @author Karsten Borgwardt
% @date May 11th 2007

% calculate results for pairs of graphs
%
%for i = 1:size(graph,2)
%%  %i
%%  am1 = graph(i).am;
%%  sp1 = floydwarshall(am1);
%%  sp1(find(sp1 == Inf)) = 0;
%%  graph(i).sp = sp1;
%
%m = size(graph(i).am,1)
%am = graph(i).am;
%graph(i).dv = zeros(64,1);
%
%for j = 1:samplesize
%  r= randperm(m+0);
%  graph(i).dv(graphlettype(am(r(1:4),r(1:4)))) = ...
%      graph(i).dv(graphlettype(am(r(1:4),r(1:4)))) + 1;
%  
%end
%graph(i).dv = graph(i).dv .* (1/samplesize);
%end
%
%for i = 1:size(graph,2)
%  %  i
%    for j = 1:size(graph,2)
%
%      P = eye(64,64);
%      
%      result(i,j) = graph(i).dv' * P * graph(j).dv; 
%    
%    end
%  end
%
%  

P = zeros(64,64);

for a = 0:1

for b = 0:1

for c = 0:1

for d = 0:1

for e = 0:1

for f = 0:1
  am = [0, a, b, c; 
           a, 0 , d , e; 
           b,d,0,f; 
          c, e, f, 0] ;
  
  P(graphlettype(am),graphlettype(am([ 1 2 3 4], [ 1 2 3 4]))) = 1;
  P(graphlettype(am),graphlettype(am([ 1 2 4 3], [ 1 2 4 3]))) = 1;
  P(graphlettype(am),graphlettype(am([ 1 3 2 4], [ 1 3 2 4]))) = 1;
  P(graphlettype(am),graphlettype(am([ 1 3 4 2], [ 1 3 4 2]))) = 1;
  P(graphlettype(am),graphlettype(am([ 1 4 2 3], [ 1 4 2 3]))) = 1;
  P(graphlettype(am),graphlettype(am([ 1 4 3 2], [ 1 4 3 2]))) = 1;

  P(graphlettype(am),graphlettype(am([ 2 1 3 4], [ 2 1 3 4]))) = 1;
  P(graphlettype(am),graphlettype(am([ 2 1 4 3], [ 2 1 4 3]))) = 1;
  P(graphlettype(am),graphlettype(am([ 2 3 1 4], [ 2 3 1 4]))) = 1;
  P(graphlettype(am),graphlettype(am([ 2 3 4 1], [ 2 3 4 1]))) = 1;
  P(graphlettype(am),graphlettype(am([ 2 4 1 3], [ 2 4 1 3]))) = 1;
  P(graphlettype(am),graphlettype(am([ 2 4 3 1], [ 2 4 3 1]))) = 1;

  P(graphlettype(am),graphlettype(am([3  1 2 4], [ 3  1 2 4]))) = 1;
  P(graphlettype(am),graphlettype(am([3  1 4 2], [ 3  1 4 2]))) = 1;
  P(graphlettype(am),graphlettype(am([3  2 1 4], [ 3  2 1 4]))) = 1;
  P(graphlettype(am),graphlettype(am([3  2 4 1], [ 3  2 4 1]))) = 1;
  P(graphlettype(am),graphlettype(am([3  4 1 2], [ 3  4 1 2]))) = 1;
  P(graphlettype(am),graphlettype(am([3  4 2 1], [ 3  4 2 1]))) = 1;

  P(graphlettype(am),graphlettype(am([4  1 2 3], [ 4  1 2 3]))) = 1;
  P(graphlettype(am),graphlettype(am([4  1 3 2], [ 4  1 3 2]))) = 1;
  P(graphlettype(am),graphlettype(am([4  2 1 3], [ 4  2 1 3]))) = 1;
  P(graphlettype(am),graphlettype(am([4  2 3 1], [ 4  2 3 1]))) = 1;
  P(graphlettype(am),graphlettype(am([4  3 1 2], [ 4  3 1 2]))) = 1;
  P(graphlettype(am),graphlettype(am([4  3 2 1], [ 4  3 2 1]))) = 1;


end
end
end
end
end
end



  function result = graphlettype(am)
  % determine graphlet type
  
  factor = 2 .^ [0:5]';
  
  upper = [am(1,2:4)';am(2,3:4)'; am(3,4)'];  
  result = sum(factor .* upper) +1;