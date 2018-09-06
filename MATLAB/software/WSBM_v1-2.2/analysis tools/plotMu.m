function [] = plotMu(Model)
%PLOTMU plots the posterior probability of vertex-label assignments. 
%
%   PLOTMU(Model) generates a IMAGESC plot of mu, the posterior probability
%   of vertex-label assigments. Each row is a vertex and each column is a
%   group. Ideally the posterior for each vertex is concentrated into one
%   group. If a vertex has a uniform or dispersed posterior, this indicates
%   a lack of fit or indecision in assigning a label.
%
%   Example:
%       [Label,Model] = wsbm(Raw_Data,k);
%       plotMu(Model)
%   Input:
%       mu   - kxn mat ~ mu(g,v) is the prob of vertex v belongs to group g
%   Output:
%       A IMAGESC plot. 
%
%   See also WSBM, PLOTWSBM

% PLOTMU
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
 

if isstruct(Model),
    mu = Model.Para.mu;
elseif isnumeric(Model),
    mu = Model;
else
    error('Unrecognized Input for Model');
end

opengl('software');
imagesc(mu');
title('Vertex-Label Posterior Plot');
xlabel('Probability in Group');
ylabel('Vertex');
colorbar
caxis([0 1]);
end