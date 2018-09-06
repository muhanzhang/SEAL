function [Labels,Model] = biwsbm(E,K_0,K_1,types,varargin)
%BIWSBM finds latent community structure in bipartite networks.
%
%   BIWSBM is a wrapper around the WSBM algorithm used for the special case 
%   of bipartite networks. It takes advantage of our prior knowledge of
%   the bipartite sstructure.
%   This algorithm infers the parameters by approximating a posterior 
%   distribution using an iterative variational Bayes algorithm.
% 
%   Syntax:
%       [Labels] = biwsbm(E,K_A,K_B,types)
%       [Labels,Model] = biwsbm(...,'ParaName',ParaValue) 
%
%   Examples:
%       A = Edg2Adj(generateEdges());
%       A = [zeros(size(A)),A;A,zeros(size(A))];
%       K_0 = 2; 
%       K_1 = 2;
%       types = zeros(length(A),1);
%       types(length(A)/2+1:end) = 1;
%       [Labels,Model] = biwsbm(A,K_0,K_1,types);
%   Other Examples:
%       % Ignore Edge Information:
%       [Labels,Model] = biwsbm(A,K_0,K_1,types,E_Distr,'None');
%       % Ignore Weight Information:
%       [Labels,Model] = biwsbm(A,K_0,K_1,types,W_Distr,'None');
%       % Use Belief Propagation Algorithm:
%       [Labels,Model] = biwsbm(A,K_0,K_1,types,algType,'bp');
%       % Change number of trials:
%       [Labels,Model] = biwsbm(A,K_0,K_1,types,numTrials,123);
%       % Run algorithm in parallel:
%       [Labels,Model] = biwsbm(A,K_0,K_1,types,parallel,1);
%
%   Inputs:
%       E - an m by 3 network edge list (parent,child,weight) 
%           (If E is an n by n network adjacency matrix, then it will be
%           converted into an edge list by ADJ2EDG)
%       K_0 - number of type 0 communities to infer
%       K_1 - number of type 1 communities to infer
%       types - an n by 1 binary vector of vertex types.
%           types(i) = 0 or 1 depending on the vertex type
%       ... - varargin (additional optional parameters)
%             Try `type biwsbm.m' for more information
%
%   Outputs:
%       Labels - a n by 1 vector of edge labels (using the MAP estimates)
%                   - Ties are broken randomly (with a warning message)
%       Model - a MATLAB structure for advanced output
%
%   For more information, try 'type biwsbm.m'
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
%   See also SETUP_DISTR, INSTALLMEXFILES, WSBM, CALC_LOGEVIDENCE

% BIWSBM
% Version 1.0 | April 2014 | Christopher Aicher
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
% BIWSBM CODE
%-------------------------------------------------------------------------%

if nargin < 4,
    error('Not enough input arguments');
end

% Modify the Prior Distribution to include Type Information
K = K_0 + K_1;
n = length(types);
mu_0 = zeros(K,n);
mu_0(1:K_0,types == 0) = 1/K_0;
mu_0(K_0+1:end,types == 1) = 1/K_1;
if sum((types ~= 0) & (types ~= 1)),
    fprintf('Warning: Variable types is not binary...\n'); 
    mu_0(:,(types ~= 0) & (types ~= 1)) = 1/K;
end
% Ignore edges between the same type.
R = zeros(K,K);
for ii = 1:K_0,
    for jj = 1:K_1,
        R(ii,jj+K_0) = ii+(jj-1)*K_0;
        R(jj+K_0,ii) = ii+(jj-1)*K_0+K_0*K_1;
    end
end

% Call wsbm_driver.m
Model = wsbm_driver(E,R,'mu_0',mu_0,varargin{:});

% Create Labels from Vertex Posterior
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
% END OF BIWSBM CODE
%-------------------------------------------------------------------------%
%EOF