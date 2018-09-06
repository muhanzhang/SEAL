function [Labels,Model] = wsbm(E,R_Struct,varargin)
%WSBM find latent community structure in weighted networks.
%
%   WSBM is the main driver program for finding community structure, 
%   inferring the vertex-labels and edge-bundle parameters of a 
%   Weighted Stochastic Block Model (WSBM). 
%   This algorithm infers the parameters by approximating a posterior 
%   distribution using an iterative variational Bayes algorithm.
%   See Aicher, Jacobs, Clauset (2013) for the theoretical derivation of 
%   the algorithm 
% 
%   Syntax:
%       [Labels] = wsbm(E)
%       [Labels] = wsbm(E,k)
%       [Labels, Model] = wsbm(...,'ParaName',ParaValue) 
%
%   Examples:
%       Raw_Data = generateEdges();
%     % Default
%       [Labels] = wsbm(Raw_Data);
%     % Infer 2 Groups
%       [Labels] = wsbm(Raw_Data,2);
%     % Change W_Distr to Exp
%       [Labels] = wsbm(Raw_Data,2,'W_Distr','Exp');
%     % Run Code in Parallel
%       [Labels] = wsbm(Raw_Data,2,'parallel',1);
%     % Ignore E_Distr
%       [Labels] = wsbm(Raw_Data,2,'E_Distr','None');
%       [Labels] = wsbm(Raw_Data,2,'alpha',0);
%     % Increase the number of trials
%       [Labels] = wsbm(Raw_Data,2,'numTrials','500');
%     % Multiple changes at once
%       [Labels] = wsbm(Raw_Data,2,'W_Distr','Exp','alpha',0,'parallel',1);
%   See WSBMDemo.m for more examples
%
%   Inputs:
%       E - an m by 3 network edge list (parent,child,weight) 
%           (If E is an n by n network adjacency matrix, then it will be
%           converted into an edge list by ADJ2EDG)
%       k - number of blocks            (k = 4 default)
%
%   Outputs:
%       Labels - a n by 1 vector of edge labels (using the MAP estimates)
%                   - Ties are broken randomly (with a warning message)
%       Model - a MATLAB structure for advanced output
%
%   For more information, try 'type wsbm.m'
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
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   See also SETUP_DISTR, INSTALLMEXFILES, WSBMDEMO, CALC_LOGEVIDENCE

% WSBM
% Version 1.0 | December 2013 | Christopher Aicher
%
% 
% ADVANCED INPUT/OUTPUT:
% -Advanced Inputs: 'ParaName' - ParaValue      * indicates default
% --Model Inputs:
%   'W_Distr' - edge-weight distr    ('Normal' *)          
%   'E_Distr' - edge-existence distr ('Bernoulli' *)        
%       To select a distribution type it's name:
%       Weighted distributions:
%           Bernoulli, Binomial, Poisson, Normal, LogNormal, 
%           Exponential, Pareto, None, 
%       Edge distributions:
%           Bernoulli, Binomial, Poisson, None, DC (Degree Corrected)
%       See SETUP_DISTR for more information
%   'R_Struct' - kxk matrix of the block structure.                
%   'alpha' - para in [0,1]. 0 only weight, 0.5* both, 1 only existence
% --Inference Options:
%   'numTrials' - number of trials with different random initial conditions
%   'algType' - 'vb'* naive bayes, 'bp' belief propagation
%   'networkType' - 'dir'* directed,'sym' = symmetric,'asym' = asymmetric
%   'nanType' - 'missing'* nans are missing, 'nonedge' = nans are nonedge 
%   'mainMaxIter' - Maximum number of iterations in main_loop
%   'mainTol' - Minimum (Max Norm) convergence tolerance in main_loop
%   'muMaxIter' - Maximum number of iterations in mu_loop
%   'muTol' - Minimum (Max Norm) convergence tolerance in mu_loop
% --Extra Options:
%   'verbosity' - 0 silent, 1* default, 2 verbose, 3 very verbose 
%   'parallel' - boolean, run in parallel? 0* No (Need Parallel ToolBox)
%   'save' - boolean, save temp results? 0* No 
%   'outputpath' - string to where to save temp results (Only if save = 1)
%   'seed' - seed for algtype (mu_0 or mes_0) (Sets numTrials = 1) 
%   'mexfile' - boolean, run using MEX files? 1* Yes (Need MEX Files)
% --Prior Options:
%   'mu_0' - kxn matrix prior vector for vertex-label parameters(sums to 1)
%
% -Advanced Outputs:
%  'Model' - struct with the following fields
%       'name' - name of model <W_Distr-E_Distr-R_Struct>
%       'Data' - struct with Data related variables
%       'W_Distr' - struct from SETUP_DISTR
%       'E_Distr' - struct from SETUP_DISTR
%       'R_Struct' - struct with edge-bundle (R) variables
%       'Para' - struct with inferred hyperparamters (tau,mu) and parameter
%                estimates (theta)
%       'Options' - struct with inference option information
%       'Flags' - struct with convergence flags
%
%
% For an overview see the README.txt file
%

%-------------------------------------------------------------------------%
% WSBM CODE
%-------------------------------------------------------------------------%

% Call wsbm_driver.m 
if nargin > 1,
    Model = wsbm_driver(E,R_Struct,varargin{:});
else
    Model = wsbm_driver(E);
end

Labels = Model.Para.mu';
[n,k] = size(Labels);
if sum(sum(Labels >= 1/k-10^-3)) > n
    fprintf('Breaking %u Ties Randomly\n',sum(sum(Labels >= 1/k-10^-3,2) > 1));
    e = [zeros(size(Labels,1),1) cumsum(Labels,2)];
    e(:,end) = 1; e(e>1) = 1;
    Labels = diff(e,1,2);
    Labels = mnrnd(1,Labels);
end
[~,Labels] = max(Labels,[],2);

%-------------------------------------------------------------------------%
% END OF WSBM CODE
%-------------------------------------------------------------------------%
%EOF