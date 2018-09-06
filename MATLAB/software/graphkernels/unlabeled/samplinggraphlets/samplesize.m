function result = samplesize(delta, epsilon, a)
% delta = confidence level (typically 0.05 or 0.1)
% epsilon = precision level (typically 0.05 or 1)
% a = number of isomorphism classes of graphlets
%
% Karsten Borgwardt
% 4/11/2008


result = 2 * ( a* log(2) + log(1/delta) ) / (epsilon^2)
