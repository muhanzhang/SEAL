function [Distr] = setup_distr(distr_name,distr_prior,distr_para)
%SETUP_DISTR creates a Distr Struct.
%   This is for advanced users. 
%
%   Syntax:
%       [Distr] = setup_distr(distr_name,distr_prior,distr_para)
%
%   SETUP_DISTR is used to create Distr Structs for wsbm and generateEdges.
%   It can be used to create both 
%       Weighted distributions:
%           Bernoulli, Binomial, Poisson, Geometric, NegBinomial, 
%           Normal, LogNormal, Exponential, Pareto, None, 
%       Edge distributions:
%           Bernoulli, Binomial, Poisson, Geometric, NegBinomial,
%           None, DC (Degree Corrected)
%
%   To access distribution help, type 'help setup_distr>Normal_Family'
%    (replace "Name" with the distr name or scroll down to 'See also')
%   
%   Examples:
%       Biased Normal Prior:
%           [W_Distr] = setup_distr('Normal',[0,1,100],1);
%       Binomial with n Trials:
%           [E_Distr] = setup_distr('Binomial',[-1/2,-1/n],n);
%       Default Poisson Prior:
%           [Distr] = setup_distr('Poisson');
%
%   Inputs:
%       distr_name - string of the distribution name
%       distr_prior - tx1 vector of prior parameters (Optional)
%  !! Inappropriate prior parameters can lead to unexpected behavior !!
%       distr_para - additional distribution parameter (Optional)
%               (e.g. number of trials for Binomial distribution)
%
%   It's recommended to vary the distr_prior to test robustness
% 
%
%   Output:
%    Distr struct for ExpFamily: f(x) = h(x)*exp(T(x)*Eta(Theta))
%     tau_0  - 1xt matrix of prior parameters (See Distr.pdf for defaults)
%     logh   - log of h (func of A)
%     T      - tx1 cell array of T (func of A)
%     Eta    - tx1 cell array of Eta (func of tau)
%     logZ   - log of Z (func of tau)
%     Theta  - estimate Theta (func of tau)
%     Predict- estimate A (func of Para) (Not Ideal Estimators)
%     name   - string of distr name (i.e. 'Normal')
%     distr  - draw r.v. (func of nrow,ncol,theta) 
%     cdf    - cdf (func of x,theta)
%    
%   See Also 
%   SETUP_DISTR>BINOMIAL_FAMILY, SETUP_DISTR>POISSON_FAMILY,
%   SETUP_DISTR>NEGBINOMIAL_FAMILY,
%   SETUP_DISTR>NORMAL_FAMILY, SETUP_DISTR>LOGNORMAL_FAMILY,
%   SETUP_DISTR>EXP_FAMILY, SETUP_DISTR>PARETO_FAMILY,
%   SETUP_DISTR>DC_FAMILY, SETUP_DISTR>NO_DISTR,
%


% Version 1.0 | December 2013 | Christopher Aicher
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

%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% SETUP_DISTR CODE
%-------------------------------------------------------------------------%

%Get Distr From Name
switch lower(distr_name),
    case {'bernoulli'} 
        % Bernoulli
        if nargin < 2,
            Distr = Binomial_Family([-1/2,-1/2],1); % Default
        else
            Distr = Binomial_Family(distr_prior,1);
        end
    case {'binomial','bino'},
        % Binomial
        if nargin < 3,
            error('Number of trials needs to be specified for the Binomial distribution');
        else
            Distr = Binomial_Family(distr_prior,distr_para);
        end
    case 'poisson',
        if nargin < 2,
            Distr = Poisson_Family(); % Default
        else
            Distr = Poisson_Family(distr_prior);
        end
    case {'geometric'}
        % Geometric
        if nargin < 2,
            Distr = NegBinomial_Family([-1/2,-1/2],1); % Default
        else
            Distr = NegBinomial_Family(distr_prior,1);
        end
    case {'negbinomial','negbino','negativebinomial'},
        %NegBinomial
        if nargin < 3,
            error('Number of successes needs to be specified for NegBinomial distribution');
        else
            Distr = NegBinomial_Family(distr_prior,distr_para);
        end
    case {'normal','gaussian'},
        if nargin < 2,
            Distr = Normal_Family();
        else
            if nargin < 3,
                Distr = Normal_Family(distr_prior);
            else
                Distr = Normal_Family(distr_prior,distr_para);
            end
        end
    case 'lognormal',
        if nargin < 2,
            Distr = LogNormal_Family();
        else
            if nargin < 3,
                Distr = LogNormal_Family(distr_prior);
            else
                Distr = LogNormal_Family(distr_prior,distr_para);
            end
        end
    case {'exponential','exp'},
        if nargin < 2,
            Distr = Exp_Family();
        else
            if nargin < 3,
                Distr = Exp_Family(distr_prior);
            else
                Distr = Exp_Family(distr_prior,distr_para);
            end
        end
    case {'pareto','logexponential','logexp'},
        if nargin < 2,
            Distr = Pareto_Family();
        else
            if nargin < 3,
                Distr = Pareto_Family(distr_prior);
            else
                Distr = Pareto_Family(distr_prior,distr_para);
            end
        end
    case {'thresh','theshold'},
        if nargin < 3,
            error('Requires Threshold Value');
        else
            Distr = Binomial_Family(distr_prior,1);
            Distr.name = sprintf('Thresh_%2.2f',distr_para);
            Distr.logh = @(A) ones(size(A))*1; 
            Distr.T{1} = @(x) (x > distr_para)*1;
            Distr.cdf = @(x,theta) binocdf((x > distr_para)*1,1,theta(:,1));
        end
    case {'threshp','theshold_percentile'},
        if nargin < 3,
            error('Requires Threshold Value');
        else
            Distr = Binomial_Family(distr_prior,1);
            Distr.name = sprintf('Thresh_Percentile_%2.2f',distr_para);
            Distr.logh = @(A) ones(size(A))*1; 
            Distr.T{1} = @(x) (x > prctile(x(:),distr_para))*1;
            Distr.cdf = @(x,theta) NaN; 
        end
    case {'none'}
        Distr = No_Distr();
    case {'dc','degree-corrected','dcpoisson'},
        if nargin < 2,
            Distr = DC_Family();
            Distr.Distr_Para = 1; % Avg Degree
        else
            if naragin < 3,
                Distr = DC_Family(distr_prior);
                Distr.Distr_Para = 1;
            else
                Distr = DC_Family(distr_prior);
                Distr.Distr_Para = distr_para;
            end
        end
    otherwise,
        error('distr_name = %s is not a valid distr name',distr_name);
end
% Save F+Distr
end % End of SETUP_DISTR

%-------------------------------------------------------------------------%
% HELPER FUNCTIONS
%-------------------------------------------------------------------------%

function [Distr] = Binomial_Family(tau_0,n) 
%BINOMIAL_FAMILY creates the Distr for the binomial distribution
%Example:
%  Uniform Prior:   [Distr] = setup_distr('binomial',[0,0],n) 
%  Reference Prior: [Distr] = setup_distr('binomial',[-1/2,-n^-1],n) 
%
%Syntax:
%  [Distr] = setup_distr('binomial',distr_prior,distr_para) 
%
%Input: 
%   distr_prior - tau_0, hyperparameters (Default: 0,0)
%   distr_para  - n, number of trials
%   !! Note that A(i,j) cannot be greater than n !!
%      For that reason it is required to input n
%   The prior distribution = Beta(tau(1)+1,n tau(2)-tau(1)+1)
%  
%Output: Distr struct.
%   Theta is for p
%   Predict is the estimate of np 
%See Also setup_distr, wsbm 
        if nargin < 2, error('n not specified for binomial distribution'); end
        if nargin < 1, tau_0 = [0,0.001]; end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Binomial tau must be 1x2.');
        end
        
        name = 'Binomial';
        distr = @(nrow,ncol,theta) binoinv(rand(nrow,ncol),n,theta);
        logh = @(A)gammaln(n+1)-gammaln(n-A+1)-gammaln(A+1);
        T_1 = @(x) x;
        N = @(x) ~isnan(x);
        T = {T_1; N};
        E_1 = @(tau) psi(tau(:,1)+1)-psi(n*tau(:,2)-tau(:,1)+1);
        negPhi = @(tau) n*(psi(n*tau(:,2)-tau(:,1)+1)-psi(n*tau(:,2)+2));
        Eta = {E_1; negPhi};
        logZ = @(tau) betaln(tau(:,1)+1,n*tau(:,2)-tau(:,1)+1);
        Theta = @(tau) (tau(:,1)+1)./(n.*tau(:,2)+2);
        Predict = @(theta) n*theta(:,1); 
        cdf = @(x,theta) binocdf(x,n,theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = Poisson_Family(tau_0) 
%POISSON_FAMILY creates the Distr for the Poisson distribution
%Example:
%  Default          [Distr] = setup_distr('poisson',[0,0.001]) 
%                                       % For numerical reasons 
%  Uniform Prior:   [Distr] = setup_distr('poisson',[0,0]) 
%  Reference Prior: [Distr] = setup_distr('poisson',[-1/2,0])
%
%Syntax:
%  [Distr] = setup_distr('poisson',distr_para) 
%
%Input: 
%   distr_para - tau_0 prior parameters (Default: 1,0.001)
%   The prior distribution = Gamma(tau(1)+1,tau(2))
%  
%Output: Distr struct.
%   Theta gives mean parameter
%   Predict is the estimate of mean  
%See Also setup_distr, wsbm 
        if nargin < 1, 
            tau_0 = [0,0.001]; 
            fprintf('It''s recommended to vary the prior for the Poisson distr.\n'); 
        end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Poisson tau must be 1x2.');
        end

        name = 'Poisson';
        distr = @(nrow,ncol,theta) poissinv(rand(nrow,ncol),theta(1));
        logh = @(A)-gammaln(A+1);
        T_1 = @(x) x;
        N = @(x) ~isnan(x);
        T = {T_1; N};
        E_1 = @(tau) psi(tau(:,1)+1)-log(tau(:,2));
        negPhi = @(tau) -(tau(:,1)+1)./tau(:,2);
        Eta = {E_1; negPhi};
        logZ = @(tau) gammaln(tau(:,1)+1)-log(tau(:,2)).*(tau(:,1)+1);
        Theta = @(tau) (tau(:,1)+1)./(tau(:,2));
        Predict = @(theta) theta(:,1); 
        cdf = @(x,theta) poisscdf(x,theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = NegBinomial_Family(tau_0,r) 
%NEGBINOMIAL_FAMILY creates the Distr for the negative binomial
% distribution
%
% Here X, a negative binomial r.v., is the number of failures before r
% successes. Theta = p = probability of success
%
%Example:
%  Uniform Prior:   [Distr] = setup_distr('negbinomial',[0,0],n) 
%  Reference Prior: [Distr] = setup_distr('negbinomial',[-1/2,-r^-1/2],r) 
%
%Syntax:
%  [Distr] = setup_distr('negbinomial',distr_prior,distr_para) 
%
%Input: 
%   distr_prior - tau_0, hyperparameters (Default: -1/2,-1/(2r))
%   distr_para  - r, number of successes
%
%   The prior distribution = Beta(r*tau(2)+1, tau(1)+1)
%  
%Output: Distr struct.
%   Theta is for p
%   Predict is the estimate of np 
%See Also setup_distr, wsbm 
        if nargin < 2, error('r not specified for negbinomial distribution'); end
        if nargin < 1, tau_0 = [0,0.001]; end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For NegBinomial tau must be 1x2.');
        end
        
        name = 'NegBinomial';
        distr = @(nrow,ncol,theta) nbininv(rand(nrow,ncol),r,theta);
        logh = @(A)gammaln(r+A-1+1)-gammaln(A+1)-gammaln(r-1+1);
        T_1 = @(x) x;
        N = @(x) ~isnan(x);
        T = {T_1; N};
        E_1 = @(tau) psi(tau(:,1)+1)-psi(tau(:,1)+r*tau(:,2)+2);
        negPhi = @(tau) r*(psi(r*tau(:,2)+2)-psi(tau(:,1)+r*tau(:,2)+2));
        Eta = {E_1; negPhi};
        logZ = @(tau) betaln(tau(:,1)+1,r*tau(:,2)+1);
        Theta = @(tau) (tau(:,1)+1)./(r.*tau(:,2)+tau(:,1)+2);
        Predict = @(theta) r*(1-theta(:,1))./theta(:,1); 
        cdf = @(x,theta) nbincdf(x,r,theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = Normal_Family(tau_0,alpha_0) 
%NORMAL_FAMILY creates the Distr for the normal distribution
%Example:
%  Default          [Distr] = setup_distr('Normal',[0,0.001,0.001],1) 
%  Uniform Prior:   [Distr] = setup_distr('Normal',[0,0,0],1) 
%  Reference Prior: [Distr] = setup_distr('Normal',[0,0,0],3)
%               or: [Distr] = setup_distr('Normal',[0,0,0],5)
%  Empirical Bayes: [Distr] = setup_distr('Normal',...
%                               [mean(Data),var(Data)+mean(data)^2,w],1)
%                    Where you pick w to be small.
%Syntax:
%  [Distr] = setup_distr('Normal',distr_prior,distr_para) 
%
%Input: 
%   distr_prior   - tau_0, prior parameters (Default: 0,0,0.001)
%   distr_para - alpha, additional prior parameter (Default: 1)
%   Note that prior distribution is 
%    mean ~ Normal(tau(1)/tau(3),sigma^2/tau(3))
%    var ~ Gamma((tau(3)+alpha)/2,(tau(3)tau(2)-tau(1)^2)/2tau(3)))
%  
%Output: Distr struct.
%   Theta gives mean and variance (not standard deviation)
%   Predict is the estimate of mean 
%   distr's theta(1) is mean, theta(2) is variance 
%See Also setup_distr, wsbm 
        if nargin < 2,alpha_0 = 1; end
        if nargin < 1, tau_0 = [0,0.001,0.001]; end
        if 3 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Normal tau must be 1x3.');
        end
        name = 'Normal';
        distr = @(nrow,ncol,theta) norminv(rand(nrow,ncol),theta(1),sqrt(theta(2)));
        logh = @(A)-1/2*log(2*pi)*ones(size(A));
        T_1 = @(x) x;
        T_2 = @(x) x.^2;
        N = @(x) ~isnan(x);
        T = {T_1; T_2; N};
        E_1 = @(tau) (tau(:,3)+alpha_0).*tau(:,1)./(tau(:,3).*tau(:,2)-tau(:,1).^2);
        E_2 = @(tau) -((tau(:,3)+alpha_0).*tau(:,3))./(2*(tau(:,3).*tau(:,2)-tau(:,1).^2));
        negPhi = @(tau) -(1/2*((tau(:,2)...
     +(1+(alpha_0-1)./tau(:,3)).*tau(:,1).^2)./(tau(:,3).*tau(:,2)-tau(:,1).^2)...
     - psi((tau(:,3)+alpha_0)/2) - log(2) + log(tau(:,2)-tau(:,1).^2./tau(:,3))));
        Eta = {E_1;E_2;negPhi};
        logZ = @(tau) -1/2*log(tau(:,3))+gammaln((tau(:,3)+alpha_0)/2)...
     +((tau(:,3)+alpha_0)/2).*(log(2)-log((tau(:,2)-tau(:,1).^2./tau(:,3)))) + 1/2*log(2*pi);
        % mu, sigma^2
        Theta = @(tau) [tau(:,1)./tau(:,3),((tau(:,2)-tau(:,1).^2./tau(:,3))./(tau(:,3)-1))];
        Predict = @(theta) theta(:,1); 
        cdf = @(x,theta) normcdf(x,theta(:,1),sqrt(theta(:,2)));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = LogNormal_Family(tau_0,alpha_0) 
%LOGNORMAL_FAMILY creates the Distr for the lognormal distribution
%Example:
%  Default          [Distr] = setup_distr('LogNormal',[0,0.001,0.001],1) 
%  Uniform Prior:   [Distr] = setup_distr('LogNormal',[0,0,0],1) 
%  Reference Prior: [Distr] = setup_distr('LogNormal',[0,0,0],3)
%               or: [Distr] = setup_distr('LogNormal',[0,0,0],5)
%  Empirical Bayes: [Distr] = setup_distr('LogNormal',...
%                      [mean(logData),var(logData)+mean(logData)^2,w],1)
%                    Where you pick w to be small.
%Syntax:
%  [Distr] = setup_distr('LogNormal',distr_prior,distr_para) 
%
%Input: 
%   distr_prior   - tau_0, prior parameters (Default: 0,0,0.001)
%   distr_para - alpha, additional prior parameter (Default: 1)
%   Note that prior distribution is 
%    logmean ~ Normal(tau(1)/tau(3),sigma^2/tau(3))
%    logvar ~ Gamma((tau(3)+alpha)/2,(tau(3)tau(2)-tau(1)^2)/2tau(3)))
%  
%Output: Distr struct.
%   Theta gives 'mean' and 'variance' (of logarithm of the data)
%   Predict is the estimate of mean of the lognormal: exp(mean+variance/2) 
%   distr's theta(1) is mean, theta(2) is variance 
%See Also setup_distr, setup_distr>Normal_Family, wsbm 
        if nargin < 2,alpha_0 = 1; end
        if nargin < 1, tau_0 = [0,0.001,0.001]; end
        if 3 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For LogNormal tau must be 1x3.');
        end
        name = 'LogNormal';
        distr = @(nrow,ncol,theta) exp(norminv(rand(nrow,ncol),theta(1),sqrt(theta(2))));
        logh = @(A)-1/2*log(2*pi)*ones(size(A));
        T_1 = @(x) log(x);
        T_2 = @(x) log(x).^2;
        N = @(x) ~isnan(x);
        T = {T_1; T_2; N};
        E_1 = @(tau) (tau(:,3)+alpha_0).*tau(:,1)./(tau(:,3).*tau(:,2)-tau(:,1).^2);
        E_2 = @(tau) -((tau(:,3)+alpha_0).*tau(:,3))./(2*(tau(:,3).*tau(:,2)-tau(:,1).^2));
        negPhi = @(tau) -(1/2*((tau(:,2)...
     +(1+(alpha_0-1)./tau(:,3)).*tau(:,1).^2)./(tau(:,3).*tau(:,2)-tau(:,1).^2)...
     - psi((tau(:,3)+alpha_0)/2) - log(2) + log(tau(:,2)-tau(:,1).^2./tau(:,3))));
        Eta = {E_1;E_2;negPhi};
        logZ = @(tau) -1/2*log(tau(:,3))+gammaln((tau(:,3)+alpha_0)/2)...
     +((tau(:,3)+alpha_0)/2).*(log(2)-log((tau(:,2)-tau(:,1).^2./tau(:,3)))) + 1/2*log(2*pi);
        % mu, sigma^2
        Theta = @(tau) [tau(:,1)./tau(:,3),((tau(:,2)-tau(:,1).^2./tau(:,3))./(tau(:,3)-1))];
        Predict = @(theta) exp(theta(:,1)+theta(:,2)/2);
        cdf = @(x,theta) logncdf(x,theta(:,1),sqrt(theta(:,2)));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = Exp_Family(tau_0,x_0) 
%EXP_FAMILY creates the Distr for the exponential distribution
%Example:
%  Default          [Distr] = setup_distr('Exp',[0.001,0],0) 
%  Uniform Prior:   [Distr] = setup_distr('Exp',[0,0],x_0) 
%  Reference Prior: [Distr] = setup_distr('Exp',[0,1],x_0) 
%
%Syntax:
%  [Distr] = setup_distr('Exp',distr_prior,distr_para) 
%
%Input: 
%   distr_prior - tau_0, prior parameters (Default: 0.001,0)
%   distr_para - x0, shift parameter (Default = 0)
%   Note that prior distribution is Gamma(tau(2)+1,tau(1))
%  
%Output: Distr struct.
%   Theta gives mean
%   Predict is the estimate of mean (approx)
%See Also setup_distr, wsbm 
        if nargin < 2,x_0 = 0; end
        if nargin < 1, tau_0 = [0.001,0]; end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Exp tau must be 1x2.');
        end
        name = 'Exp';
        distr = @(nrow,ncol,theta) expinv(rand(nrow,ncol),theta)+x_0;
        logh = @(x) 0;
        T_1 =@(x) x-x_0;
        N = @(x) ~isnan(x);
        T = {T_1; N};
        E_1 = @(tau) -(tau(:,2)+1)./tau(:,1);
        negPhi = @(tau) psi(tau(:,2)+1) - log(tau(:,1));
        Eta = {E_1; negPhi};
        logZ = @(tau) gammaln(tau(:,2)+1) - (tau(:,2)+1).*log(tau(:,1));
        Theta = @(tau) tau(:,1)./(tau(:,2)+1);
        Predict = @(theta) theta(:,1); 
        cdf = @(x,theta) expcdf(x,theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = Pareto_Family(tau_0,x_m) 
%PARETO_FAMILY creates the Distr for the Pareto (logexp) distribution,
%    shifted by 1 / centered at zero.
%Example:

%  Default          [Distr] = setup_distr('Exp',[0.001,0],1) 
%  Uniform Prior:   [Distr] = setup_distr('Exp',[0,0],x_m) 
%  Reference Prior: [Distr] = setup_distr('Exp',[0,1],x_m) 
%
%Syntax:
%  [Distr] = setup_distr('Exp',distr_prior,distr_para) 
%
%Input: 
%   distr_prior - tau_0, prior parameters (Default: 0.001,0)
%   distr_para - x_m, shift parameter (Default = 1)
%   Note that prior distribution is Gamma(tau(2)+1,tau(1))%Input: 
%         i.e. log((x+1)./x_m) follows exp distr,
%         or x_m*exp(x)-1 follows pareto distribution 
%  
%Output: Distr struct.
%   Theta gives "mean" of logarithm of data
%   Predict is the estimate of mean (approx/assumes rate is greater than 1)
%See Also setup_distr, setup_distr>Exp_Family, wsbm 
        if nargin < 2,x_m = 1; end
        if nargin < 1, tau_0 = [0.001,0]; end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Pareto tau must be 1x2.');
        end
        name = 'Pareto';
        distr = @(nrow,ncol,theta) exp(expinv(rand(nrow,ncol),theta))*x_m-1;
        logh = @(x) 0;
        T_1 =@(x) log((x+1)./x_m);
        N = @(x) ~isnan(x);
        T = {T_1; N};
        E_1 = @(tau) -(tau(:,2)+1)./tau(:,1);
        negPhi = @(tau) psi(tau(:,2)+1) - log(tau(:,1));
        Eta = {E_1; negPhi};
        logZ = @(tau) gammaln(tau(:,2)+1) - (tau(:,2)+1).*log(tau(:,1));
        Theta = @(tau) tau(:,1)./(tau(:,2)+1);
        Predict = @(theta) theta(:,1).*x_m./(theta(:,1)-1); %Not exactly accurate 
        cdf = @(x,theta) expcdf(log((x+1)./x_m),theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end

function [Distr] = DC_Family(tau_0) 
%DC_FAMILY creates the Distr for the Degree Corrected Poisson distribution
%   THIS IS FOR ADVANCED USERS ONLY. - There are some subtleties that are
%   easily overlooked.
%
%Example:
%  Default          [Distr] = setup_distr('dc',[0,0.001]) 
%  Uniform Prior:   [Distr] = setup_distr('dc',[0,0]) 
%  Reference Prior: [Distr] = setup_distr('dc',[-1/2,0])
%
%Syntax:
%  [Distr] = setup_distr('dc',distr_prior) 
%
%Input: 
%   distr_prior   - tau_0, prior parameters (Default: 1,0.001)
%   Note that the prior distribution is Gamma(tau(1)+1,tau(2))
%  
%Output: Distr struct.
%  !! Theta gives mean multiplied by the degree of the nodes 
%  !! Predict is the estimate of mean divided by the degree of the nodes
%  !!!!!!!!! These are off by a factor of each nodes DEGREE !!!!!!!!!!!!!!
%See Also setup_distr, wsbm 
        if nargin < 1, 
            tau_0 = [0,0.001]; 
            fprintf('It''s recommended to vary the prior for the Poisson distr.\n'); 
            fprintf('DC loglikelihood estimates are correct up to a constant!\n'); 
        end
        if 2 ~= size(tau_0,2) || 1 ~= size(tau_0,1),
            error('Invalid tau_0. For Poisson tau must be 1x2.');
        end
        name = 'DCPoisson';
        distr = @(nrow,ncol,theta) poissinv(rand(nrow,ncol),theta);
        logh = @(A)-gammaln(A+1); % Correct up to a constant
        T_1 = @(x) x;
        N = @(Edge_List,degree) (degree(2,Edge_List(:,1)).*degree(1,Edge_List(:,2)));
        T = {T_1; N};
        E_1 = @(tau) psi(tau(:,1)+1)-log(tau(:,2));
        negPhi = @(tau) -(tau(:,1)+1)./tau(:,2);
        Eta = {E_1; negPhi};
        logZ = @(tau) gammaln(tau(:,1)+1)-log(tau(:,2)).*(tau(:,1)+1);
        Theta = @(tau) (tau(:,1)+1)./(tau(:,2));
        Predict = @(theta) theta(:,1); 
        cdf = @(x,theta) poisscdf(x,theta(:,1));
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end
    
function [Distr] = No_Distr()
%NO_DISTR creates an empty Distr for the case where we assume no
% distribution
%
%Syntax:
%  [Distr] = setup_distr('None')
%
%Output: Distr struct.
%  !! Theta gives NaN 
%  !! Predict returns NaN
%See Also setup_distr, wsbm 
        name = 'None';
        distr = @(nrow,ncol,theta) NaN(nrow,ncol);
        tau_0 = [];
        logh = @(A) 0;
        T = {};
        Eta = {};
        logZ = @(tau) 0;
        Theta = @(tau) NaN;
        Predict = @(theta) NaN; 
        cdf = @(x,theta) NaN;
        Distr = struct('name',name,'distr',distr,'tau_0',tau_0,...
            'logh',logh,'T',{T},'Eta',{Eta},'logZ',logZ,'Theta',Theta,...
            'Predict',Predict,'cdf',cdf);
end
    
%%% End Of Distribution Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% End of SETUP_DISTR.M