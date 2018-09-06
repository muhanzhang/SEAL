function [Best_Models,Scores,Models] = wsbmLooper(E,ModelInputs,scorefuncs)
%WSBMLOOPER - a wrapper for wsbm.m fitting multiple models
%
%   WSBMLOOPER runs the wsbm.m for various models specified by 
%   the cell array ModelInputs 
%   
%   Syntax:
%       [Best_Model,LogEvidences,...] = wsbmLooper(E,ModelInputs)
%       [Best_Models,Scores,...] = wsbmLooper(E,ModelInputs,scorefunc)
%       [Best_Models,Scores,...] = wsbmLooper(E,ModelInputs,scorefuncs)
%       [...,Models] = wsbmLooper(E,ModelInputs)
%
%   Inputs: 
%       E - an m by 3 network edge list (parent,child,weight) 
%           (If E is an n by n network adjacency matrix, then it will be
%           converted into an edge list by Adj2Edg.m)
%       ModelInputs - an Mx1 cell array of cells containing the inputs 
%           arguments to wsbm.m. 
%       scorefunc - a score function taking the Model struct as input
%       scorefuncs - an Ax1 cell array of score functions
%
%   Outputs:
%       Best_Model - the best model determined by maximizing the score
%       Best_Models - an 1xA cell array of best models: one per score func
%       Scores - an MxA matrix of scores (Row - Model, Column - Score func)
%       Models - an Mx1 cell array of fitted models
%
%   Examples:
%    [E,True_Model] =  generateEdges();
%    Infer1 = {4,'W_Distr','None','E_Distr','Bernoulli'};
%    Infer2 = {4,'W_Distr','Normal','E_Distr','None'};
%    Infer3 = {4,'W_Distr','Normal','E_Distr','Bernoulli'};
%    ModelInputs = {Infer1;Infer2;Infer3};
%    [Best_Model,LogEvidence] = wsbmLooper(E,ModelInputs);
%    
%    score1 = @(model) model.Para.LogEvidence; % Log Evidence
%    score2 = @(model) nmi(model,True_Model); % NMI vs True Labels
%    score3 = @(model) varInfo(model,True_Model); % VI vs True Labels
%    scorefuncs = {score1;score2;score3};
%    [Best_Models,Scores,Models] = wsbmLooper(E,ModelInputs,scorefuncs);
%    

%   See Also WSBM, SETUP_DISTR, CALC_LOGEVIDENCE 

% WSBMLOOPER
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
% Check Inputs
if nargin < 3, scorefuncs = @(model) model.Para.LogEvidence; end
if ~iscell(scorefuncs),
    %Convert scorefunc to scorefuncs
    scorefuncs = {scorefuncs};
end

% Setup Variables
num_models = numel(ModelInputs);
num_scores = numel(scorefuncs);
Scores = zeros(num_models,num_scores);
BestScores = -Inf*ones(1,num_scores);
Best_Models = cell(1,num_scores);
if nargout > 2,
    Models = cell(num_models,1);
end

% Fit Each Model
for mm = 1:num_models,
    fprintf('\n----- Fitting Model %u of %u ------\n',mm,num_models);
    try
        [~,Model] = wsbm(E,ModelInputs{mm}{:});
    catch exception
        fprintf('Error in ModelInputs{%d}:\n',mm);
        error(exception.message)
    end
    % Score Each Model 
    for aa = 1:num_scores,
        Scores(mm,aa) = scorefuncs{aa}(Model);
        if Scores(mm,aa) > BestScores(aa),
            Best_Models{aa} = Model;
            BestScores(aa) = Scores(mm);
        end
    end
    if nargout > 2,
        Models{mm} = Model;
    end
end
if num_scores == 1, 
    % If only one score function, just return the Model Struct
    Best_Models = Best_Models{1};
end
fprintf('wsbmLooper Done... \n');

