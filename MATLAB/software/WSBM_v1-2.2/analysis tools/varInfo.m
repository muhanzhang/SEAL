function [vi] = varInfo(MU1,MU2)
%VARINFO calculates the variation of information (VI) between two 
%   vertex-labelings. For more on VI see:
%   http://www.sciencedirect.com/science/article/pii/S0047259X06002016
%
%   Examples:
%       vi = varInfo(MU1,MU2);
%       vi = varInfo(Model1,Model2);
%       vi = varInfo(Model1.Para.mu,Model2.Para.mu);
%
%   Input:
%       MU1 - k1 x n ~ mu(g,v) is the prob of vertex v belongs to group g
%       MU2 - k2 x n ~ mu(g,v) is the prob of vertex v belongs to group g
%       Model1,Model2 - Model structs form WSBM
%   Output:
%       vi - a real valued metric of divergence between distribution 
%           Indentical Clusters <-> vi = 0
%
%       vi = H - I, where H is joint entropy and I is mutual information
%
%   See also WSBM, NMI

% varInfo comes with ABSOLUTELY NO WARRANTY
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

if isstruct(MU1),
    MU1 = MU1.Para.mu;
end
if isstruct(MU2),
    MU2 = MU2.Para.mu;
end
if size(MU1,2) ~= size(MU2,2),
    error('MU1 and MU2 must both have the same number of vertices (columns)');
end

[~,n] = size(MU1);
n_ij = MU1*MU2';
n_1 = sum(MU1,2);
n_2 = sum(MU2,2);
n_tot = n_1*n_2'; 
I = n_ij > 0.001;
vi = -1/n*sum(n_ij(I).*log(n_ij(I).^2./n_tot(I)));
end
