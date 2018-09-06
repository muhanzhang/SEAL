function [new_data,shuffle_order] = shuffle(old_data)
%SHUFFLE is a function that randomly permutes the vertex labels of the 
% input network (either an nxn adjacency matrix or mx3 edge list)
%
% [new_data] = shuffle(old_data)
% [shuffled_Model] = shuffle(True_Model);
%
%INPUTS:
%    old_data - either an nxn adjacency matrix, an mx3 edge list, or Model
%               struct
%
%OUTPUTS:
%    new_data - a randomly shuffled adjacency matrix, edge list, or Model
%               struct (same as input)
% See Also GENERATEDATA

% SHUFFLE comes with ABSOLUTELY NO WARRANTY
% Version 1.0 | December 2013 | Christopher Aicher

if isstruct(old_data),
    % Model Struct
    if isfield(old_data,'Data'),
        if isfield(old_data.Data,'Raw_Data'),
            new_data = old_data;
            [new_data.Data.Raw_Data,shuffle_order] = shuffle(old_data.Data.Raw_Data);
            if isfield(old_data,'Para'),
                if isfield(old_data.Para,'mu'),
                    new_data.Para.mu(:,shuffle_order) = old_data.Para.mu;
                end
            end
        else
            error('Missing field: .Data.Raw_Data in input');
        end
    else
        error('Missing field: .Data.Raw_Data in input');
    end
else
    [m,n] = size(old_data);
    if m == n,
        % Adj Matrix
        shuffle_order = randperm(n);
        temp = old_data(shuffle_order,:);
        new_data = temp(:,shuffle_order);
    elseif n == 3,
        % Edge List
        n = max(max(old_data(:,1:2)));
        shuffle_order = randperm(n);
        new_data = old_data;
        new_data(:,1) = shuffle_order(old_data(:,1));
        new_data(:,2) = shuffle_order(old_data(:,2));
    else
        error('old_data is not an nxn adj_mat nor an mx3 edg_list');
    end
end

end