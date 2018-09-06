function [weight_Error] = predictW_Error(Model,Test_List,Weight_Error)
%PREDICTW_ERROR - calculates the error of Model prediction on the Test_List
%   
%   Syntax:
%    [Weight_Error] = predictW_Error(Model,Test_List,Weight_Error)
%
%   Input:
%       Model - WSBM struct 
%       Test_List - m_test x 3 edge list of test data
%       Weight_Error - error function for weights
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
%   See Also WSBM, WSBMLOOPER, CROSSVALIDEDG, PREDICTE_ERROR 

% predictW_Error
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
if isnan(Model.Para.theta_w), weight_Error = NaN; return; end;

% Predict Weights
weighted_edges = ~isnan(Test_List(:,3));
I = Test_List(weighted_edges,1);
J = Test_List(weighted_edges,2);
Mu_I = Model.Para.mu(:,I);
Mu_J = Model.Para.mu(:,J);
weight_Predict = reshape(Model.W_Distr.Predict(...
    Model.Para.theta_w(Model.R_Struct.R(:),:)),...
    size(Model.R_Struct.R,1),size(Model.R_Struct.R,2)); 
weight_Predict(isnan(weight_Predict)) = 0;
weight_list = sum(Mu_I'*weight_Predict.*Mu_J',2);
weight_Error = mean(Weight_Error(weight_list,Test_List(weighted_edges,3)));


