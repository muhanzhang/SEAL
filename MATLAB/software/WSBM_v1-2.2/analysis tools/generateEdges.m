function [Edge_List,True_Model] = generateEdges(W_truth,E_truth,R,theta_w,theta_e,group_Size,degree_Para)
%GENERATEEDGES generates synthetic data sets from the WSBM.
%   This is for advanced users.
%
%   Examples:
%       Generate an default random example edge list
%   [Edge_List,True_Model] = generateEdges() 
%       Generic form:
%   [Edge_List,True_Model] =
%       generateEdges(W_Distr,E_Distr,R,theta_w,theta_e,group_Size)
%   [Edge_List,True_Model] =
%       generateEdges(W_Distr,E_Distr,R,theta_w,theta_e,group_Size,degree_Para)
%
%   Inputs:
%       W_truth - Distr Struct or Distr Name              (see setup_distr)
%       E_truth - a *discrete* Distr Struct or Distr Name (see setup_distr)
%       R       - a kxk matrix linking pairs of groups to edge bundles
%       theta_w   - rxd mat of parameters for each weight distr
%       theta_e   - rxd mat of parameters for each edge distr
%       group_Size - kx1 vec of group/block sizes (sum is n)
%       degree_Para - nx2 vec or mean in- and out- degrees (DC models only)
%
%   Outputs:
%       Edge_List - mx3 mat of Raw_Data values 
%   	True_Model - a struct containing: Data,W_Distr,E_Distr,R_Struct,Para     
%
%   Advanced Examples:
%   % Normal SBM4:
%       W = setup_distr('Normal');
%       E = setup_distr('Bernoulli');
%       R = reshape(1:16,4,4);
%       theta_w = [1,1; 2,1; 3,1; 4,1;
%                  5,1; 6,5; 7,1; 8,1; 
%                  9,1; 10,1; 11,10; 12,1;
%                  13,1; 14,1; 15,1; 16,5];
%       theta_e = [0.1; 0.2; 0.3; 0.4;
%                  0.2; 0.3; 0.4; 0.1;
%                  0.3; 0.4; 0.1; 0.2;
%                  0.4; 0.1; 0.2; 0.3];
%       group_Size = [30;30;30;30];
%       [Edge_List,True_Model] = generateData(W,E,R,...
%                                   theta_w,theta_e,group_Size);
%       plotWSBM(Edge_List);
%
%   % 2 Weighted Groups 2 Edge Groups:
%       R = [1,2,3,4; 
%            2,1,4,3; 
%            3,4,1,2;
%            4,3,2,1];
%       theta_w = [100,1;100,1; 0,1; 0,1];
%       theta_e = [0.1; 0.9; 0.1; 0.9];
%       group_sizes = [25;25;25;25];
%       [~,True_Model] = generateEdges('Normal','Bernoulli',...
%                           R,theta_w,theta_e,group_sizes);
%       plotWSBM(True_Model);
%
%   % Degree Corrected Groups:
%       R = [1,2; 3,4];
%       theta_w = [];
%       theta_e = [1; 4; 2; 3];
%       group_sizes = [50;50];
%       degree_Para = 2*rand(sum(group_sizes),2);
%       [~,True_Model] = generateEdges('None','DC',R,theta_w,theta_e,...
%           group_sizes,degree_Para);
%       plotWSBM(True_Model,'edges');
%
%   For more information, try 'type generateEdges.m'
% 
%   See Also WSBM, SETUP_DISTR, WSBMDEMO, SHUFFLE

% GENERATEEDGES
% Version 1.0 | December 2013 | Christopher Aicher
% Version 1.1 | April 2014 | Christopher Aicher | Added DC Options
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
% fprintf('Generating Edge List...\n');
if nargin == 0,
    % Defaults
    W_truth = 'Normal';
    E_truth = 'Bernoulli';
    R = [1,2; 
         3,4]; 
    theta_w = [-10,1;-5,10;5,10; 10,1];
    theta_e = [0.9;0.8;0.7;0.6];
    group_Size = [25;25];
end
% Check Inputs
checkInputs(); 
%-------------------------------------------------------------------------
% Setup mu_truth
n = sum(group_Size);
breaks = [0;cumsum(group_Size)];
mu_truth = zeros(n,R_truth.k);
for krow = 1:R_truth.k,
            mu_truth((breaks(krow)+1):breaks(krow+1),krow) = 1;
end

% Setup Edges
NumEdges = zeros(n,n);
for krow = 1:R_truth.k,
    for kcol = 1:R_truth.k,
rows = (breaks(krow)+1):breaks(krow+1);
cols = (breaks(kcol)+1):breaks(kcol+1);            
switch E_truth.name
    case {'DCPoisson'},
        NumEdges(rows,cols) = E_truth.distr(length(rows),length(cols),...
theta_e(R_truth.R(krow,kcol))*degree_Para(rows,:)*degree_Para(cols,:)');
    case {'None'}
        NumEdges(rows,cols) = 1;
    otherwise
        NumEdges(rows,cols) = E_truth.distr(length(rows),length(cols),...
            theta_e(R_truth.R(krow,kcol)));
end
    end
end
NumEdges(isnan(NumEdges(:))) = 0;

% Setup Weights
Edge_List = zeros(sum(NumEdges(:)),3);

ee = 1;
for ii = 1:n,
    for jj = 1:n,
        if NumEdges(ii,jj),
    Edge_List(ee:(ee+NumEdges(ii,jj)-1),1) = ii;
    Edge_List(ee:(ee+NumEdges(ii,jj)-1),2) = jj;
    ee = ee+NumEdges(ii,jj);
        end
    end
end
if ~strcmpi(W_truth.name,'None'),
    for krow = 1:R_truth.k,
        for kcol = 1:R_truth.k,
            EE = find( breaks(krow) < Edge_List(:,1) & ...
                       Edge_List(:,1) <= breaks(krow+1) & ...
                       breaks(kcol) < Edge_List(:,2) & ...
                       Edge_List(:,2) <= breaks(kcol+1));
            Edge_List(EE,3) = W_truth.distr(...
                length(EE),1,theta_w(R_truth.R(krow,kcol),:));
        end
    end
else
    Edge_List(:,3) = 1;
end          
    
            
    

% Create True_Model
Para = struct('mu',mu_truth','theta_w',theta_w,'theta_e',theta_e);
Data.Raw_Data = Edge_List;
True_Model = struct('name','Truth','Data',Data,'W_Distr',W_truth,'E_Distr',E_truth,...
                    'R_Struct',R_truth,'Para',Para);
% fprintf('...Done Generating Edge List\n');


%%% Nested Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function checkInputs()          
        % Check F_Distr
        if ~isstruct(W_truth),
            if ischar(W_truth),
                W_truth = setup_distr(W_truth);
            else
                error('W_truth needs to be an distr Struct or name');
            end
        end
        if ~isfield(W_truth,'distr'),
            error('W_truth is missing a distr function');
        end
        % Check E_Distr
        if ~isstruct(E_truth),
            if ischar(E_truth),
                E_truth = setup_distr(E_truth);
            else
                error('E_truth needs to be an distr Struct or name');
            end
        end
        if ~isfield(E_truth,'distr'),
            error('E_truth is missing a distr function');
        end
        if ~isstruct(R),
            [k1,k2] = size(R);
            if k1 ~= k2,
                error('R must be a square kxk matrix');
            end
            k = k1; 
            r = max(R(:));
            R_truth = struct('R',R,'r',r,'k',k);
        else
            R_truth = R;
        end
        % Check R_Struct
        if ~isstruct(R_truth),
            error('R_truth needs to be an R_Struct Struct');
        end
        if ~isfield(R_truth,'R'),
            error('R_truth is missing matrix R');
        end
        if ~isfield(R_truth,'k'),
            error('R_truth is missing value r');
        end
        if ~isfield(R_truth,'k'),
            error('R_truth is missing value k');
        end
        
        % Check theta_w
        if size(theta_w,1) ~= R_truth.r,
            if ~strcmpi(W_truth.name,'None'),
                error('theta_w must be an rxd matrix of parameters');
            end
        end
        % Check theta_e
        if size(theta_e,1) ~= R_truth.r,
            if ~strcmpi(E_truth.name,'None'),
                error('theta_e must be an rxd matrix of parameters');
            end
        end
        % Check group_Sizes
        if size(group_Size,1) ~= R_truth.k,
            group_Size = group_Size';
            if size(group_Size,1) ~= R_truth.k,
                error('group_Size must be a kx1 vector of group sizes');
            end
        end
        % Check degree_Para
        if exist('degree_Para','var'),
            [n1,n2] = size(degree_Para);
            if (n1 ~= sum(group_Size)) || (n2 ~= 2),
                error('degree_Para must be an nx2 vector of degree propensities');
            end
            if sum(degree_Para < 0),
                error('All elements of degree_Para must be nonnegative');
            end
        end
    end

%%% End of Nested Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


