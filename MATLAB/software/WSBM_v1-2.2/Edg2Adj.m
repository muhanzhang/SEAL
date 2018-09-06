function [Adj_matrix,NumEdges] = Edg2Adj(Edge_List)
%EDG2ADJ converts an Edge List into an Adjacency Matrix 
%   
%   EDG2ADJ(Edge_List) converts Edge_List into an average adjacency matrix
%   Adj_matrix and counts the number of edges in the adjacency matrix
%   NumEdges. 
%
%   Examples:
%       [Adj_matrix,NumEdges] = Edg2Adj(Edge_List)
%
%   Inputs:
%       Edg_list - mx3 matrix of (parent, child, weight)
%   Outputs:
%       Adj_matrix - nxn matrix, weighted average with NaN for non-edges
%       NumEdges - nxn matrix, counts of the edges
%      
%   See also ADJ2EDG

% EDG2ADJ
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


[m,m2] = size(Edge_List);
if(m == m2),
    fprintf('Edge_List input appears to already be an Adj_matrix\n');
    Adj_matrix = Edge_List;
    NumEdges = ones(m,m2);
    return;
elseif m2 ~=3,
    error('Edge_List is not an mx3 matrix');
end

n = max(max(Edge_List(:,1:2)));
NumEdges = zeros(n,n);
Adj_matrix = zeros(n,n);
for ee = 1:m,
    NumEdges(Edge_List(ee,1),Edge_List(ee,2)) = ...
        NumEdges(Edge_List(ee,1),Edge_List(ee,2))+1;
    Adj_matrix(Edge_List(ee,1),Edge_List(ee,2)) = ...
        Adj_matrix(Edge_List(ee,1),Edge_List(ee,2))+Edge_List(ee,3); 
end

Adj_matrix = Adj_matrix./NumEdges;
Adj_matrix(NumEdges(:) == 0) = NaN;

end