function [] = plotWSBM(Raw_Data,Labels,style)
%PLOTWSBM plots the data (Raw_Data) sorted by the vertex-labels (Labels). 
% 
%   PLOTWSBM(Raw_Data,Labels) plots the Raw_Data as an adjacency matrix
%   sorted by the vertex-labels. The adjacency matrix's rows and columns
%   are permuted such that vertices in the same group are next to each 
%   other. This allows us to visualize the `block' structure. 
%   Alternative modes allow us to look at the edge count or weight on a
%   log scale.
%
%   Examples:
%   Regular Plot
%       [Label,Model] = wsbm(Raw_Data,k);
%       plotWSBM(Raw_Data); 
%       plotWSBM(Raw_Data,Labels);
%       plotWSBM(Raw_Data,Model.Para.mu);
%       plotWSBM(Model);
%   Log Plot
%       [Label,Model] = wsbm(Raw_Data,k);
%       plotWSBM(Raw_Data,Label,'log10');
%       plotWSBM(Model,'log10');
%   Edge Plot 
%       [Label,Model] = wsbm(Raw_Data,k);
%       plotWSBM(Raw_Data,Label,'edge');
%       plotWSBM(Model,'edge');
%
%
%   Inputs:
%       Raw_Data - mx3 (edge list) or nxn (adjacency matrix) of the network
%       Labels - 1xn vector of vertex-labels
%       mu  - kxn matrix ~ mu(g,v) is the prob of vertex v belonging to 
%                          group g
%       style - optional string for selecting 
%                   logplot - 'log10' (base10) or 'log' (base e)
%                               Note: values <= 0 are discarded
%                   edge plot - 'edge'
%   Outputs:
%       Plot of Raw_Data ordered by group labels.
%   See also WSBM, PLOTMU

% PLOTWSBM
% Version 1.0 | December 2013 | Christopher Aicher
% Version 1.1 | April 2014 | Christopher Aicher | added different options
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


% TODO: This is not written very intuitively due to cut and pasting old code.
if isstruct(Raw_Data),
    if nargin == 2 && ischar(Labels),
        style = Labels;
    end
    Model = Raw_Data;
    Raw_Data = Model.Data.Raw_Data;
    Labels = Model.Para.mu;
end
[n,m] = size(Raw_Data);
if n == m,
    fprintf('Treating Raw_Data as an Adj Matrix\n');
elseif m == 3,
    fprintf('Treating Raw_Data as an Edge List\n');
    n = max([Raw_Data(:,1);Raw_Data(:,2)]);
else
    error('Invalid Raw_Data Format');
end
if ~exist('Labels','var'),
    mu = ones(1,n);
elseif ~isnumeric(Labels),
    if ischar(Labels),
        style = Labels;
        mu = ones(1,n);
    else
        error('Labels is not a numeric array.');
    end
elseif size(Labels,1) == 1,
    mu = zeros(max(Labels),size(Labels,2));
    for kk = 1:max(Labels),
        mu(kk,Labels == kk) = 1;
    end

else
    mu = Labels;
end
if ~exist('style','var'),
    style = 'default';
end

opengl('software');

mu = mu';
[n,k] = size(mu);
A_sort = zeros(n);
list = zeros(1,n);
breaks = zeros(1,k);

if sum(sum(mu > 0)) > n
    fprintf('\nRounding mu, randomly picking uniform\n');
    mu = mu./(sum(mu,2)*ones(1,size(mu,2)));
    e = [zeros(size(mu,1),1) cumsum(mu,2)];
    e(:,end) = 1;
    e(e>1) = 1;
    mu = diff(e,1,2);
    mu = mnrnd(1,mu);
end

cur = 1;
for ii = 1:k
    indicies = find(mu(:,ii));
    list(cur:cur+length(indicies)-1) = indicies;
    cur = cur + length(indicies);
    breaks(ii) = cur-1;
end

switch lower(style),
    case {'edge','edges'}
        if m == 3,
            [~,Raw_Data] = Edg2Adj(Raw_Data);
        else
            Raw_Data = ones(size(Raw_Data));
        end
    case {'log10'}
        if m == 3,
            [Raw_Data,~] = Edg2Adj(Raw_Data);
        end
        Raw_Data(Raw_Data <= 0) = NaN; 
        Raw_Data = log10(Raw_Data);
    case {'log','loge'}
        if m == 3,
            [Raw_Data,~] = Edg2Adj(Raw_Data);
        end
        Raw_Data(Raw_Data <= 0) = NaN;
        Raw_Data = log(Raw_Data);
    otherwise
        if m == 3,
            [Raw_Data,~] = Edg2Adj(Raw_Data);
        end
end

for ii = 1:n
    A_sort(ii,:) = Raw_Data(list(ii),list);
end

    

%Plot the Matrix
h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
xlabel('Child'); 
ylabel('Parent');
set(h,'alphadata',~isnan(A_sort));
colorbar;
hold all; 

for ii = 1:k-1
    plot([breaks(ii)+.5,breaks(ii)+.5],[-.5,breaks(k)+.5],'-k','LineWidth',1.5);
    plot([-.5,breaks(k)+.5],[breaks(ii)+.5,breaks(ii)+.5],'-k','LineWidth',1.5);
end        
hold off;
    
end