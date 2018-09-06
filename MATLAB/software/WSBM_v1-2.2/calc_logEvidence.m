function [LogEvidence] = calc_logEvidence(Data,W_Distr,E_Distr,R_Struct,Para,Options)
%CALC_LOGEVIDENCE approximates the log-evidence of model, logPr(E|Model) 
%
%   CALC_LOGEVIDENCE(Model) approximates the log-evidence of the model 
%   using a variational Bayes lower bound described in 
%   Aicher, Jacobs, Clauset (2013). This approximate log-evidence is used 
%   for fitting a weighted stochastic block model in the WSBM.m algorithm 
%   and can be used for model selection (e.g. the number of blocks k).   
%
%   Examples:
%       [Labels_1, Model_1] = wsbm(Raw_Data,k1);
%       LogEvidence_1 = calc_logEvidence(Model_1)
%       [Labels_2, Model_2] = wsbm(Raw_Data,k2);
%       LogEvidence_2 = calc_logEvidence(Model_2)
%
%   Input:
%       Model - The annonated output of WSBM
%
%   Outputs:
%       LogEvidence - logPr(E|Model) 
%           the log-evidence or marginal log-likelihood of the model           
%           Note that LogEvidence can be positive for continuous r.v.s
%           (log-probability density, not a log-probability)
%
%   For more information, try 'type calc_logEvidence'   
%
%   See also WSBM

% CALC_LOGEVIDENCE
% Version 1.0 | December 2013 | Christopher Aicher
%   -Currently using VB approximation for BP.
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
%
% Input Requires the following fields:
%     Data, W_Distr, E_Distr, R_Struct, Para, Options
%
% Output Calculation:
%     log Pr(E|Model) = 1-alpha_w (sum_r log(Z(tau_w)/Z(tau_w0)) + logHw)
%                       + alpha_w (sum_r log(Z(tau_e)/Z(tau_e0)) + logHe)
%                       + sum_i sum_k mu_ik log(mu_0k/mu_ik)
%
%

%-------------------------------------------------------------------------%
% CALC_LOGEVIDENCE CODE
%-------------------------------------------------------------------------%
if nargin == 1, % Parse a Model Struct
    Model = Data;
    Data = Model.Data; 
    W_Distr = Model.W_Distr;
    E_Distr = Model.E_Distr;
    R_Struct = Model.R_Struct;
    Para = Model.Para;
    Options = Model.Options;
end

switch lower(Options.algType)
    case {'vb'}
        % Weighted Distribution Component (tau_w)
        if ~isempty(Para.tau_w),
            logW = Data.logHw;
            logW_0 = W_Distr.logZ(W_Distr.tau_0);
            % Which Edge Bundles were updated:
            updated = Para.tau_w(:,end) > 10^-6+W_Distr.tau_0(end);
            logW = logW+sum(W_Distr.logZ(Para.tau_w(updated,:)));
            if ~isnan(logW_0) && ~isinf(logW_0),
                logW = logW-sum(updated)*logW_0;
            else
                if Options.verbosity > 1,
                    fprintf('logW_0 is not nice %u\n',logW_0);
                end
            end
        else
            logW = 0;
        end
        % Edge Distribution Component (tau_e)
        if ~isempty(Para.tau_e),
            logE = Data.logHe;
            logE_0 = E_Distr.logZ(E_Distr.tau_0);
            % Which Edge Bundles were updated:
            updated = Para.tau_e(:,end) > 10^-6+E_Distr.tau_0(end);
            logE = logE+sum(E_Distr.logZ(Para.tau_e(updated,:)));
            if ~isnan(logE_0) && ~isinf(logE_0),
                logE = logE-sum(updated)*logE_0;
            else
                if Options.verbosity > 1,
                    fprintf('logE_0 is not nice %u\n',logE_0);
                end
            end
        else
            logE = 0;
        end
        % Latent Group Assignment Component (mu)
        mu_0 = Options.mu_0;
        Index = Para.mu > 10^-6;
        logMu = -sum(sum(Para.mu(Index).*log(Para.mu(Index)./mu_0(Index))));
        % Calculate Log Evidence
        LogEvidence = (1-Options.alpha)*logW+Options.alpha*logE+logMu;
    case {'bp'} % Should Fix This Later
        % Weighted Distribution Component (tau_w)
        if ~isempty(Para.tau_w),
            logW = Data.logHw;
            logW_0 = W_Distr.logZ(W_Distr.tau_0);
            % Which Edge Bundles were updated:
            updated = Para.tau_w(:,end) > 10^-6+W_Distr.tau_0(end);
            logW = logW+sum(W_Distr.logZ(Para.tau_w(updated,:)));
            if ~isnan(logW_0) && ~isinf(logW_0),
                logW = logW-sum(updated)*logW_0;
            else
                if Options.verbosity > 1,
                    fprintf('logW_0 is not nice %u\n',logW_0);
                end
            end
        else
            logW = 0;
        end
        % Edge Distribution Component (tau_e)
        if ~isempty(Para.tau_e),
            logE = Data.logHe;
            logE_0 = E_Distr.logZ(E_Distr.tau_0);
            % Which Edge Bundles were updated:
            updated = Para.tau_e(:,end) > 10^-6+E_Distr.tau_0(end);
            logE = logE+sum(E_Distr.logZ(Para.tau_e(updated,:)));
            if ~isnan(logE_0) && ~isinf(logE_0),
                logE = logE-sum(updated)*logE_0;
            else
                if Options.verbosity > 1,
                    fprintf('logE_0 is not nice %u\n',logE_0);
                end
            end
        else
            logE = 0;
        end
        % Latent Group Assignment Component (mu)
        mu_0 = Options.mu_0;
        Index = Para.mu > 10^-6;
        logMu = -sum(sum(Para.mu(Index).*log(Para.mu(Index)./mu_0(Index))));
        % Calculate Log Evidence
        LogEvidence = (1-Options.alpha)*logW+Options.alpha*logE+logMu;
        
    otherwise
      error('Unrecognized Options.algType');
end
  

%-------------------------------------------------------------------------%
end % END OF CALC_LOGEVIDENCE.m
%-------------------------------------------------------------------------%
%EOF
