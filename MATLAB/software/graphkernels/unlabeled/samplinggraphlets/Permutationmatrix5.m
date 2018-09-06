function P = Permutationmatrix5

P = zeros(1024,1024);

for a = 0:1

for b = 0:1

for c = 0:1

for d = 0:1

for e = 0:1

for f = 0:1

for g = 0:1

for h = 0:1

for i = 0:1

for j = 0:1

  am = [0, a, b, c, d;
        a, 0, e, f, g;
        b, e, 0, h, i;
	   c, f, h, 0, j;	 
        d, g, i, j, 0] ;
  
perm=perms([1 2 3 4 5]);
for k=1:120
  P(graphlettype(am), graphlettype(am(perm(k,:),perm(k,:)) )) = 1;    
end

end
end
end
end
end
end
end
end
end
end



  function result = graphlettype(am)
  % determine graphlet type
  
  factor = 2 .^ [0:9]';
  
  upper = [am(1,2:5)';am(2,3:5)'; am(3,4:5)'; am(4,5)'];  
  result = sum(factor .* upper) +1;
