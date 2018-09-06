function P = Permutationmatrix3

P = zeros(8,8);

for a = 0:1

for b = 0:1

for c = 0:1

 am = [0, a, b;
       a, 0, c;
       b, c, 0;];
  
perm=perms([1 2 3]);
for k=1:6
  P(graphlettype(am), graphlettype(am(perm(k,:),perm(k,:)) )) = 1;    
end

end
end
end

end

function result = graphlettype(am)
  % determine graphlet type
  
  factor = 2 .^ [0:2]';
  
  upper = [am(1,2:3)';am(2,3)'];  
  result = sum(factor .* upper) +1;
end
