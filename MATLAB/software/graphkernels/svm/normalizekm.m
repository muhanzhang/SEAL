function result = normalizekm(km)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% normalizes kernelmatrix km 
% such that diag(result) = 1, i.e. K(x,y) / sqrt(K(x,x) * K(y,y))
% @author Karsten Borgwardt
% @date June 3rd 2005
% all rights reserved 

nv = sqrt(diag(km));
nm =  nv * nv';
knm = nm .^ -1;
for i = 1:size(knm,1)
for j = 1:size(knm,2)
if (knm(i,j) == Inf)
  knm(i,j) = 0;
end
end
end
result = km .* knm;
