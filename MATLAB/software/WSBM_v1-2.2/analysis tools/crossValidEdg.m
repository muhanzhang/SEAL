function [Train_List, Test_List] = crossValidEdg(Edge_List,Frac_Test)
%CROSSVALIDEDG splits an edge list into a train list and test list
%   
%   Syntax:
%    [Train_List, Test_List] = crossValidEdg(Edge_List,Frac_Test)
%
%   Input:
%       Edge_List - m x 3 matrix 
%       Frac_Test - real value between 0-1 determines the test set
%                default is 0.2 (20%)
%
%   Outputs:
%       Train_List - m_train x 3 edge list of train data
%       Test_List - m_test x 3 edge list of test data
%
% Note: Since Edge_List needs to know n = number of nodes, the last
%       edge is always in the Training Set
%       m_train = floor(Frac_Test*n^2)
%
%   Examples:
%    [E,True_Model] =  generateEdges();
%    [Train_List, Test_List] = crossValidEdg(E,0.2);
%    
%   See Also WSBM, WSBMLOOPER, PREDICTERROR 

% CROSSVALIDEDG
% Version 1.0 | Mar 2014 | Christopher Aicher
%
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

%-------------------------------------------------------------------------

if(nargin < 2), Frac_Test = 0.2; end;

[m,~] = size(Edge_List);
n = max(max(Edge_List(:,1:2)));
m_test = floor(Frac_Test*n^2);

% Setup Lists
List = randperm(n^2);
List = List(1:m_test);

Row_List = ceil(List/n);
Col_List = List-(Row_List-1)*n;

Test_List = zeros(m_test,3);
Test_List(:,1) = Row_List;
Test_List(:,2) = Col_List;

Train_List = zeros(m+m_test,3);
Train_List(1:m,:) = Edge_List;

Temp_List = Edge_List(:,1)+n*Edge_List(:,2);

to_add = 0;
to_remove = 0;
for ee = 1:m_test,
    ind = find(Temp_List == Test_List(ee,1)+n*Test_List(ee,2));
    if(isempty(ind)),
        to_add = to_add+1;
        Train_List(m+to_add,1) = Row_List(ee);
        Train_List(m+to_add,2) = Col_List(ee);
        Train_List(m+to_add,3) = NaN;
        Test_List(ee,3) = NaN;
    else
        to_remove = to_remove+1;
        Train_List(ind,3) = NaN;
        Test_List(ee,3) = mean(Edge_List(ind,3));        
    end
end
Train_List = Train_List(1:(end-to_remove),:);


end
