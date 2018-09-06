function [NMI] = nmi(MU1,MU2,type)
%NMI calculates the normalized mutual information (NMI) between two 
%   vertex-labelings. For more on NMI see:
%   http://jmlr.org/papers/volume11/vinh10a/vinh10a.pdf
%
%   Examples:
%       NMI = nmi(MU1,MU2);
%       NMI = nmi(Model1,Model2);
%       NMI = nmi(MU1,MU2,'sqrt');
%
%   Input:
%       MU1 - k1 x n ~ mu(g,v) is the prob of vertex v belongs to group g
%       MU2 - k2 x n ~ mu(g,v) is the prob of vertex v belongs to group g
%       Model1, Model2 - Model structs form WSBM
%       type - a string indicating the type of NMI calculated:
%           'sum' ~ 2I/(H1+H2)          *(Default)
%           'max' ~ I/max(H1,H2)
%           'joint' ~ I/H
%           'sqrt' ~ I/sqrt(H1*H2)
%               where I is mutual information, 
%                     H is joint entropy,
%                     H1 and H2 are the cluster entropies
%   Output:
%       nmi - a real valued [0,1] measure of information between clusters
%               Independent Clusters <-> NMI = 0
%               Identitcal Clusters <-> NMI = 1
%
%   See also WSBM, VARINFO

% nmi comes with ABSOLUTELY NO WARRANTY
% Version 1.0 | February 2014 | Christopher Aicher
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

if nargin < 3, type = 'sum'; end;
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
% Calculate Entropy
H_1 = -1/n*sum(n_1(n_1 > 0.001).*log(n_1(n_1 > 0.001)/n));
H_2 = -1/n*sum(n_2(n_2 > 0.001).*log(n_2(n_2 > 0.001)/n));
H = -1/n*sum(n_ij(n_ij > 0.001).*log(n_ij(n_ij > 0.001)/n));
% Calculate Mutual Information
I = 1/n*sum(n_ij(n_ij > 0.001).*log(n*n_ij(n_ij > 0.001)./n_tot(n_ij > 0.001)));

% Calculate NMI
NMI = NaN;
switch lower(type)
    case 'sum',
        NMI = 2*I/(H_1+H_2);
    case 'max',
        NMI = I/max([H_1,H_2]);
    case 'sqrt',
        NMI = I/sqrt(H_1*H_2);
    case 'joint',
        NMI = I/H;
    otherwise
        error('Unrecognized type: %s',type);
end

end
