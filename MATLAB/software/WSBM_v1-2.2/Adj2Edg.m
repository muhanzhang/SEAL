function [Edge_list] = Adj2Edg(Adj_matrix)
%ADJ2EDG converts an Adjacency Matrix into an Edge List 
%   
%   ADJ2EDG(Adj_matrix,symmetric) converts Adj_matrix into an edge list 
%   Edge_list treating NaNs as non-edges. 
%
%   Examples:
%       [Edge_list] = Adj2Edg(Adj_matrix)
% 
%   Inputs:
%       Adj_matrix - nxn matrix, with NaN for non-edges
%       symmetric - 0 directed (default), 1 symmetric
%
%   Outputs:
%       Edge_list - mx3 matrix of (parent, child, weight)
%      
%   See also EDG2ADJ

% ADJ2EDG
% Version 1.0 | December 2013 | Christopher Aicher
%   Copyright 2013-2014 Christopher Aicher
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>
[ntemp,mtemp] = size(Adj_matrix);
if ntemp ~= mtemp,
    if mtemp == 3,
        fprintf('Warning: Input already an Edge List\n');
        Edge_list = Adj_matrix;
        return
    end
    error('Adj Matrix is not square.\n');
end

[Is,Js] = find(~isnan(Adj_matrix));
Edge_list = [Is,Js,Adj_matrix(~isnan(Adj_matrix(:)))];


end