function [] = plotCDF(Model, block, text)
%PLOTCDF - a tool for visualizing the fit of the Model. 
%
%   PLOTCDF plots the empirical cdf (blue) versus the fitted cdf 
%   (green) for each block. Empty blocks are blank. Large deviations 
%   between the distributions implies a lack of model fit. 
%   Although measuring the `fit' requires more theoretical analysis.
%
%   Examples:
%       plotCDF(Model,3)   % Plot the 3rd block without labels
%       plotCDF(Model,R,1) % Plot all blocks with labels
%       plotCDF(Model)
%   Syntax:
%       plotCDF(Model, block, text)
%   Inputs:
%       Model - Output of wsbm_driver or generateData
%       block - nxm mat of blocks/groups values to plot. 
%           Each element of block must be in {1,2,...,r}.
%       text  - boolean for whether to add text labels to each graph
%   Outputs:
%       nxm subplot of cdf comparison plots 
%   See Also WSBM, SETUP_DISTR 

% PLOTCDF
% Version 1.0 | Jan 2014 | Christopher Aicher
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
if nargin < 3, text = 0; end
if nargin < 2, block = Model.R_Struct.R; end

[n,m] = size(block);
for mm = 1:m,
    for nn = 1:n,
% For each block
b = block(nn,mm);
% Get Data for block b
[row,col] = find(Model.R_Struct.R == b);
k1 = Model.Para.mu(row,:) > .9; 
k2 = Model.Para.mu(col,:) > .9;
y = [];
A = Edg2Adj(Model.Data.Raw_Data);
for ii = 1:length(col),
    ytemp = A(k1(ii,:),k2(ii,:));
    y = [y;ytemp(:)];
end

y = y(~isnan(y));
% Plot Results
if ~isempty(y),
    subplot(n,m,mm+(nn-1)*n);
    [f,x] = ecdf(y);
    hold on;
    plot(x,f,'b');
    if ~isnan(Model.Para.theta_w(1)),
        ff = Model.W_Distr.cdf(x,Model.Para.theta_w(b,:));
        plot(x,ff,'g');
    end
    ylim([0 1.01]);
    if text,
        legend('Empirical CDF','Fitted CDF','Location','NorthWest');
        xlabel('x');
        ylabel('CDF');
    end
end
% End For Loops
    end
end

% Plot title
ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Weighted Empirical CDF (blue) vs Fitted CDF (green)','FontSize',14);


