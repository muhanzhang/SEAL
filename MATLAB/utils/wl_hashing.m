function equivalence_classes = wl_hashing(A, labels)
%  Usage: hashing-based wl
%  --Input--
%  -A: original adjacency matrix of the enclosing subgraph
%  -labels: initial colors
%  --Output--
%  -equivalence_classes: final labels
%  
%  Original paper:
%    Kersting, K., Mladenov, M., Garnett, R., and Grohe, M. Power
%    Iterated Color Refinement. (2014). AAAI Conference on Artificial
%    Intelligence (AAAI 2014).
%  *author: Roman Garnett
%%

  % if no labels given, use initially equal labels
  if (nargin < 2)
    labels = ones(size(A, 1), 1);
  end

  equivalence_classes = zeros(size(labels));

  % iterate WL transformation until stability
  while (~is_equivalent(labels, equivalence_classes))
    equivalence_classes = labels;
    labels = wl_transformation(A, labels);
  end

end




% checks whether two labelings are equivalent under permutation
function x = is_equivalent(labels_1, labels_2)

  % until proven otherwise  
  x = false;

  m = max(labels_1);

  % different number of labels
  if (m ~= max(labels_2))
    return;
  end

  % check whether every equivalence class remains unchanged
  for i = 1:m
    y = labels_2(labels_1 == i);
    if (any(y(1) ~= y))
      return;
    end
  end

  x = true;
  
end




function new_labels = wl_transformation(A, labels)

  % the ith entry is equal to the 2^(i-1)'th prime
  primes_arguments_required = [2, 3, 7, 19, 53, 131, 311, 719, 1619, ...
                      3671, 8161, 17863, 38873, 84017, 180503, 386093, ...
                      821641, 1742537, 3681131, 7754077, 16290047];

  num_labels = max(labels);

  % generate enough primes to have one for each label in the graph
  primes_argument = primes_arguments_required(ceil(log2(num_labels)) + 1);
  p = primes(primes_argument);
  log_primes = log(p(1:num_labels))';

  signatures = labels + A * log_primes(labels);

  % map signatures to integers counting from 1
  [~, ~, new_labels] = uniquetol(signatures);

end
