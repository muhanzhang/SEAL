%% Installing WSBM MEX files
% Run INSTALLMEXFILES to create MEX files in the private folder 
%   for WSBM to use
% Q&A:
% Q: What are MEX files?
% A: MATLAB's way of mixing complied C++ code with MATLAB. 
%    See 'mex -help'. Also see:
%    http://www.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html
%
% Q: Why should I install MEX files?
% A1: Although MATLAB is fast at matrix computation, the inner computation
%     loops of WSBM can be replaced with C++ for a significant speed up.
% A2: The WSBM algorithm using MATLAB scales quadratically in the number of
%     vertices. The WSBM algorithm using MEX files scales linearly in the
%     number of edges. For sparse graphs, this can be significant.
%
% Q: How do I install MEX files?
% A: MATLAB is can be a bit picky about the specific complier and operating
%    system. Try 'mex -setup' to select a compatible complier. 
%    For more information see: 
%    http://www.mathworks.com/help/matlab/matlab_external/building-mex-file
%    s.html
%

% Version 1.0 | December 2013 | Christopher Aicher
% Version 1.1 | May 2014 |
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

% For Mac OS X 10.9 users it may be necessary to add `-Dchar16_t=uint16_t'
% to the end of each mex call Thanks to David Darmon for pointing this out. 
%
% Alternatively, For Mac OS X users it maybe necessary to edit the
% mexopts.sh file. See: http://www.mathworks.com/matlabcentral/answers/103258-mex-on-mavericks-with-r2012b#answer_112685


% Install VB MEX Files
fprintf('Installing create_T_w.cpp...\n');
mex('-outdir','private',['private',filesep,'create_T_w.cpp'])
fprintf('... Done\n');

fprintf('Installing create_T_e.cpp...\n');
mex('-outdir','private',['private',filesep,'create_T_e.cpp'])
fprintf('... Done\n');

fprintf('Installing calc_T_w_bra.cpp...\n');
mex('-outdir','private',['private',filesep,'calc_T_w_bra.cpp'])
fprintf('... Done\n');

fprintf('Installing calc_T_e_bra.cpp...\n');
mex('-outdir','private',['private',filesep,'calc_T_e_bra.cpp'])
fprintf('... Done\n');

fprintf('Installing vb_wsbm.cpp...\n');
mex('-outdir','private',['private',filesep,'vb_wsbm.cpp'])
fprintf('... Done\n');

% Install BP MEX files
fprintf('Installing create_T_bp.cpp...\n');
mex('-outdir','private',['private',filesep,'create_T_bp.cpp'])
fprintf('... Done\n');

fprintf('Installing bp_wsbm.cpp...\n');
mex('-outdir','private',['private',filesep,'bp_wsbm.cpp'])
fprintf('... Done\n');

fprintf('--Done Installing MEX files.--\n\n');

