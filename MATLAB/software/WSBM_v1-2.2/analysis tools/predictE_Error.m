function [edge_Error] = predictE_Error(Model,Test_List,Edge_Error)
%PREDICTE_ERROR - calculates the error of Model prediction on the Test_List
%   
%   Syntax:
%    [edge_Error] = predictE_Error(Model,Test_List,Edge_Error)
%
%   Input:
%       Model - WSBM struct 
%       Test_List - m_test x 3 edge list of test data
%       Edge_Error - error function for edges
%
%   Outputs:
%       weight_Error - average error on weights in the Test_List
%
%   Examples:
%    [E,True_Model] =  generateEdges();
%    [Train_List, Test_List] = crossValidEdg(E,0.2);
%    
%
%    
%   See Also WSBM, WSBMLOOPER, CROSSVALIDEDG, PREDICTW_ERROR 

% PREDICTE_ERROR
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
if isnan(Model.Para.theta_e), edge_Error = NaN; return; end;

% Predict Edges
I = Test_List(:,1);
J = Test_List(:,2);
Mu_I = Model.Para.mu(:,I);
Mu_J = Model.Para.mu(:,J);
edge_Predict = reshape(Model.E_Distr.Predict(...
    Model.Para.theta_e(Model.R_Struct.R(:),:)),...
    size(Model.R_Struct.R,1),size(Model.R_Struct.R,2)); 
edge_Predict(isnan(edge_Predict)) = 0;
edge_list = sum(Mu_I'*edge_Predict.*Mu_J',2);
edge_Error = mean(Edge_Error(edge_list,~isnan(Test_List(:,3))));



