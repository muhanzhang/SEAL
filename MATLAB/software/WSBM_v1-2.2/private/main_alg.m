function [Para,Flags] = main_alg(Data,W_Distr,E_Distr,R_Struct,Seed,Options)
%MAIN_ALG is a single run of a variational algorithm to infer the 
% parameters of the WSBM using a specified weight (W_Distr) and 
% edge (E_Distr) distribution and partitioning (R_Struct). 
%
%   The algorithm consists of two nested-loops. 
%       1. The main loop alternates two main steps:
%           -Updating the distribution posterior (tau) 
%           -Updating the vertex-label posterior (mu) 
%         Until convergence 
%       2. The mu update contains the second loop that is run iteratively 
%         until convergence. The exact update equation (VB, BP, etc.)
%         is dependent on the algorithm type specified in Options.algType. 
%
%   Syntax:
%       [Para,Flags] = main_alg(Data,W_Distr,E_Distr,R_Struct,Seed,Options)
%
%   Input:
%       Data     - struct of data,   
%       W_Distr  - struct of weight distr,  
%       E_Distr  - struct of edge distr,
%       R_Struct - struct of R
%       Seed     - kxn mat ~ seed value for mu
%       Options  - struct of options
%
%   Output:
%    Para  - Struct containing the fields
%       mu - k by n matrix, mu(z,i) = prob vertex i is in group z
%       tau_w - r by t_w matrix, tau of W_Distr for each edge bundle
%       tau_e - r by t_e matrix, tau of E_Distr for each edge bundle
%       theta_w - r by dim(theta) matrix, estimate for each bundle
%       theta_e - r by dim(theta) matrix, estimate for each bundle
%       predict_w - r by 1 vector of predicted weighted edge values
%       predict_e - r by 1 vector of predicted edge-existence prob.
%       LogEvidence - Marginal Log-likelihood (aka Log-Evidence) 
%           For more information on `theta' and `predict',
%               see the associated help file in setup_distr.m
%    Flags - Struct containing the fields
%       mainConv - bool, convergence flag for main loop
%       mainDiff - max abs difference in mu values for main loop
%       mainIter - number of iterations
%       muConv - bool, convergence flag for mu loop
%       muDiff - max abs difference in mu values for mu loop	
%
%
%   See Also WSBM, Private/VB_WSBM, Private/BP_WSBM, 
%       Private/CALC_T_E_BRA, Private/CALC_T_W_BRA

% v0.1 - Christopher Aicher - 2/20/2013
% v0.2 - Christopher Aicher - 7/25/2013 Replaced inner loop with Mex Files
% v1.0 - Christopher Aicher - 11/30/2013 Revision
% v1.1 - Christopher Aicher - 12/14/2014 Documentation + Allow non-MEXFiles
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
        
%-------------------------------------------------------------------------%
% MAIN_ALG CODE
%-------------------------------------------------------------------------%

% Setup Flags
Flags = struct('mainConv',0,'mainDiff',Inf,'mainIter',0,'muConv',0,'muDiff',Inf);

% Setup Initial Parameters
mu = Seed*1; % Make A Deep Copy
mu_old = mu*1; % Make A Deep Copy
if ~isempty(W_Distr.tau_0)
    tau_w0 = ones(R_Struct.r,1)*W_Distr.tau_0;
    tau_w = tau_w0;
else
    tau_w = zeros(R_Struct.r,0); % Empty Case
end
if ~isempty(E_Distr.tau_0),
    tau_e0 = ones(R_Struct.r,1)*E_Distr.tau_0;
    tau_e = tau_e0;
else
    tau_e = zeros(R_Struct.r,0); % Empty Case
end
Para = struct('mu',[],'tau_w',[],'tau_e',[],'theta_w',[],...
    'theta_e',[],'predict_w',[],'predict_e',[],'LogEvidence',[]);
if strcmpi(Options.algType,'bp'),
    mes = cell(Data.n,1);
    update_mes();
end

% Main Loop
for main_loop = 1:Options.mainMaxIter,
    if Options.verbosity > 1,
        fprintf('Main Loop: %u of %u\n',main_loop,Options.mainMaxIter); 
    end

    % Update Tau
    if(Options.mexFile == 1),
        update_tau();
    else
        update_tau_NoMex();
    end
    
    % Update Eta_bra
    [Eta_w_bra,Eta_e_bra] = update_eta();
    
    % Update Mu
    switch lower(Options.algType),
        case {'vb'},
            if(Options.mexFile == 1),
                vb_update_mu();
            else
                vb_update_mu_NoMex();
            end
        case {'bp'},
            if(Options.mexFile == 1),
                bp_update_mu();
            else
                error('BP Algorithm needs MEX files to work. See INSTALLMEXFILES.m to install MEX files.');
            end
        otherwise
            error('Unrecognized algType: ''%s''',Options.algType);
    end
    
    
    % Print Info
    if Options.verbosity > 1,
        Para.mu = mu;
        Para.tau_w = tau_w;
        Para.tau_e = tau_e;
%         if strncmpi(Options.algType,'bp',2),
%             Para.mes = mes;    
%             Para.mus = mus;
%         end 
% 		calc_logEvidence(Data,F_Distr,E_Distr,R_Struct,Para,Options);
    end
    
    % Check Main Loop Convergence
    Flags.mainDiff = max(max(abs(mu_old-mu)));
    if Options.verbosity > 2,
        fprintf('Main Max Diff: %4.4f\n',Flags.mainDiff); 
    end
    if Flags.mainDiff < Options.mainTol,
        Flags.mainConv = 1;
        break;
    end
    mu_old = mu*1; % Make A Deep Copy
end % End Main Loop
% Display Warning if mu does not converge in main loop
if (Options.verbosity > 0) && (~Flags.mainConv),
    fprintf('Algorithm did not converge in %u iterations\n',Options.mainMaxIter);
    fprintf('Max Difference was %4.4e\n', Flags.mainDiff); 
%     keyboard;
end

% Update Tau One Last Time
if(Options.mexFile == 1),
    update_tau();
else
    update_tau_NoMex();
end

% Format Output
Para.mu = mu;
% if strncmpi(Options.algType,'bp',2),
%     Para.mes = mes;
%     Para.mus = mus;
% end 
Para.tau_w = tau_w;
Para.tau_e = tau_e;
if ~isempty(tau_w)
    Para.theta_w = W_Distr.Theta(tau_w);
    Para.predict_w = W_Distr.Predict(Para.theta_w);
    Para.theta_w(tau_w(:,end) < 10^-2+tau_w0(:,end)) = NaN;
    Para.predict_w(tau_w(:,end) < 10^-2+tau_w0(:,end)) = NaN;
else
    Para.theta_w = NaN;
    Para.predict_w = NaN(R_Struct.r,1);
end
if ~isempty(tau_e)
    Para.theta_e = E_Distr.Theta(tau_e);
    Para.predict_e = E_Distr.Predict(Para.theta_e);
    Para.theta_e(tau_e(:,end) < 10^-2+tau_e0(:,end)) = NaN;
    Para.predict_e(tau_e(:,end) < 10^-2 +tau_e0(:,end)) = 0.5;
else
    Para.theta_e = NaN;
    Para.predict_e = NaN(R_Struct.r,1);
end


%-------------------------------------------------------------------------%
% NESTED FUNCTIONS
%-------------------------------------------------------------------------%

% UPDATE_TAU --------------------------------------------------------------
function update_tau()
%UPDATE_TAU calculate the prior of theta's parameters 

    if ~isempty(W_Distr.tau_0),
        % Calculate <T_w>
        T_w_bra = calc_T_w_bra(mu,Data.T_w_out,R_Struct.R,R_Struct.r);
        % Update tau_w
        tau_w = tau_w0 + T_w_bra';
    end
    if ~isempty(E_Distr.tau_0),
        % Calculate <T_e>
        if strncmpi(E_Distr.name,'dc',2),
            %Degree Corrected Case
T_e_bra = calc_T_e_bra(mu,Data.T_e_out,R_Struct.R,R_Struct.r,Data.degrees_w,1);
T_e_last = sum(mu.*(ones(R_Struct.k,1)*Data.degrees_w(2,:)),2)*...
    sum(mu.*(ones(R_Struct.k,1)*Data.degrees_w(1,:)),2)';
        else
            %Regular Case
T_e_bra = calc_T_e_bra(mu,Data.T_e_out,R_Struct.R,R_Struct.r,Data.degrees_w,0);
T_e_last = sum(mu,2)*sum(mu,2)'; 
        end
        % Correct the last term in <T_e>
        for kk = 1:R_Struct.k,
            for k2 = 1:R_Struct.k,
                if R_Struct.R(kk,k2) > 0,
T_e_bra(end,R_Struct.R(kk,k2)) = T_e_bra(end,R_Struct.R(kk,k2))+T_e_last(kk,k2);
                end
            end
        end
        % Update tau_e
        tau_e = tau_e0 + T_e_bra';
    end
end

% UPDATE_TAU_NOMEX --------------------------------------------------------------
function update_tau_NoMex()
%UPDATE_TAU_NOMEX calculate the prior of theta's parameters w/out MEX files
    % Update tau_w
    for tw = 1:numel(W_Distr.T),
        T_w_bra = zeros(R_Struct.r,1);
        for k1 = 1:R_Struct.k,
            for k2 = 1:R_Struct.k,
T_w_bra(R_Struct.R(k1,k2)) = T_w_bra(R_Struct.R(k1,k2)) +...
    sum(sum(mu(k1,:)*Data.T_w{tw}*mu(k2,:)'));
            end
        end
        tau_w(:,tw) = tau_w0(:,tw)+T_w_bra;
    end
    % Update tau_e
    for te = 1:numel(E_Distr.T),
        T_e_bra = zeros(R_Struct.r,1);
        for k1 = 1:R_Struct.k,
            for k2 = 1:R_Struct.k,
T_e_bra(R_Struct.R(k1,k2)) = T_e_bra(R_Struct.R(k1,k2)) +...
    sum(sum(mu(k1,:)*Data.T_e{te}*mu(k2,:)'));
            end
        end
        tau_e(:,te) = tau_e0(:,te)+T_e_bra;
    end
end

% UPDATE_ETA --------------------------------------------------------------
function [Eta_w_bra, Eta_e_bra] = update_eta()
%UPDATE ETA - calculate Eta_bras
    % Calculate Eta_w_bra
    Eta_w_bra = zeros(size(tau_w));
    for te = 1:size(tau_w,2),
        Eta_w_bra(:,te) = W_Distr.Eta{te}(tau_w);
    end
    if ~isempty(tau_w),
        Eta_w_bra(tau_w(:,end) < 10^-6+W_Distr.tau_0(end),:) = 0; % If no update from prior
    end
    Eta_w_bra = Eta_w_bra*(1-Options.alpha);
       
    % Calculate Eta_e_bra
    Eta_e_bra = zeros(size(tau_e));
    for te = 1:size(tau_e,2),
        Eta_e_bra(:,te) = E_Distr.Eta{te}(tau_e);
    end
    if ~isempty(tau_e),
        Eta_e_bra(tau_e(:,end) < 10^-6+E_Distr.tau_0(end),:) = 0; % If no update from prior
    end
    Eta_e_bra = Eta_e_bra*(Options.alpha);
end

% VB_UPDATE_MU ------------------------------------------------------------ 
function vb_update_mu()
%VB_UPDATE_MU - updates mu using a mean-field approximation

    % Run mu_loop
    if strncmpi(E_Distr.name,'dc',2),
        % Degree Corrected Case
        [mu_new,Flags.muConv,Flags.muDiff] = vb_wsbm(mu,Data.T_w_in,Data.T_w_out,...
            Data.T_e_in,Data.T_e_out,R_Struct.R,Eta_w_bra,Eta_e_bra,...
            Options.muMaxIter,Options.muTol,Options.verbosity,...
            Data.degrees_w,1,Options.mu_0);
    else
        % Normal Case
        [mu_new,Flags.muConv,Flags.muDiff] = vb_wsbm(mu,Data.T_w_in,Data.T_w_out,...
            Data.T_e_in,Data.T_e_out,R_Struct.R,Eta_w_bra,Eta_e_bra,...
            Options.muMaxIter,Options.muTol,Options.verbosity,...
            Data.degrees_w,0,Options.mu_0);     
    end
    mu = mu_new;    
    % Print Convergence Info
    if Flags.muDiff > Options.muTol && Options.verbosity > 1,
        fprintf('Mu Did not converge. Mu Diff = %2.2f\n',Flags.muDiff);
    end
    
end

% VB_UPDATE_MU_NOMEX ------------------------------------------------------------ 
function vb_update_mu_NoMex()
%VB_UPDATE_MU_NOMEX - updates mu using a mean-field approximation w/out MEX
    
    Flags.muConv = 0;
    mu_new = mu;
    % Mu Loop
    for mu_loop = 1:Options.muMaxIter,
        for ii = 1:Data.n,
            % Calculate mu_i
            logMu = zeros(R_Struct.k,1);
            % Calculate <dT>*<Eta>
            for tw = 1:numel(W_Distr.T),
                % Weights in
                TW_ = Data.T_w{tw}(:,ii);
                for k1 = 1:R_Struct.k,
logMu(k1) = logMu(k1) + Eta_w_bra(R_Struct.R(:,k1),tw)'*(mu_new*TW_);
                end
                % Weights out
                TW_ = Data.T_w{tw}(ii,:)';
                TW_(ii) = 0;
                for k1 = 1:R_Struct.k,
logMu(k1) = logMu(k1) + Eta_w_bra(R_Struct.R(k1,:),tw)'*(mu_new*TW_);
                end
            end
            for te = 1:numel(E_Distr.T),
                % Edges in
                TE_ = Data.T_e{te}(:,ii);
                for k1 = 1:R_Struct.k,
logMu(k1) = logMu(k1) + Eta_e_bra(R_Struct.R(:,k1),te)'*(mu_new*TE_);
                end
                % Edges out
                TE_ = Data.T_e{te}(ii,:)';
                TE_(ii) = 0;
                for k1 = 1:R_Struct.k,
logMu(k1) = logMu(k1) + Eta_e_bra(R_Struct.R(k1,:),te)'*(mu_new*TE_);
                end
            end    
            logMu = logMu-max(logMu); % For Numerics
            mu_ii = Options.mu_0(:,ii).*exp(logMu); % Prior Information
            mu_ii = mu_ii./sum(mu_ii); % To Normalize
            mu_new(:,ii) = mu_ii;
        end
        % Check Convergence
        Flags.muDiff = max(max(abs(mu-mu_new)));
        if Flags.muDiff < Options.muTol,
            mu = mu_new;
            Flags.muConv = 1;
            break;
        end
        mu = mu_new;
    end
    % Print Convergence Info
    if Flags.muDiff > Options.muTol && Options.verbosity > 1,
        fprintf('Mu Did not converge. Mu Diff = %2.2f\n',Flags.muDiff);
    end
    
end


% UPDATE_MES --------------------------------------------------------------
function update_mes()
%UPDATE_MES - updates mes to uniform messages
    for ii = 1:Data.n,
        temp = Data.T_w_in{ii};
        mes_ii = ones(1+R_Struct.k,size(temp,2));
        mes_ii(1,:) = temp(1,:);
        mes{ii} = mes_ii;
    end
end

% BP_UPDATE_MU ------------------------------------------------------------ 
function bp_update_mu()
%BP_UPDATE_MU - updates mu using a loopy belief propagation approximation 
    update_mes();
    [mu_new,Flags.muConv,Flags.muDiff] = ...
        bp_wsbm(mu,mes,Data.T_w_in,Data.T_w_out,Data.T_e_in,Data.T_e_out,...
        R_Struct.R,Eta_w_bra,Eta_e_bra,Options.muMaxIter,Options.muTol,...
        Options.verbosity,Options.mu_0);
    mu = mu_new;    
    % Print Convergence Info
    if Flags.muDiff > Options.muTol && Options.verbosity > 1,
        fprintf('Mu Did not converge. Mu Diff = %2.2f\n',Flags.muDiff);
    end
end

%-------------------------------------------------------------------------%
end % End of MAIN_ALG.m
%-------------------------------------------------------------------------%

%EOF