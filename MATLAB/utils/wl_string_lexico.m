% the classical string-based Weisfeiler-Lehman algorithm (a naive implementation)

function equivalence_classes = wl_string_lexico(A, labels)

  % if no labels given, use initially equal labels
  if (nargin < 2)
    labels = ones(size(A, 1), 1);
  end

  equivalence_classes = zeros(size(labels));

  % iterate WL transformation until stability
  i=0;
  while ~isequal(labels,  equivalence_classes)
    i = i+1
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

  num_labels = max(labels);
  upbound = num_labels + 1;
  neighbors = bsxfun(@times, A, labels');
  neighbors(neighbors == 0) = upbound;
  neighbors = sort(neighbors, 2);
  neighbors(neighbors == upbound) = 0;
  signatures = [labels, neighbors];
  % map signatures to integers counting from 1
  [~, ~, new_labels] = unique(signatures, 'rows');

end

