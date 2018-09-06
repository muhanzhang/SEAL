function [T_w_in,T_w_out] = create_T_w(Edge_List,degrees_w)
%CREATE_T_W formats the Edge_List into an Adjancency List of 
% Weight Distribution Statistics 
%
%   CREATE_T_W is a helper function for wsbm.m
%
%   Syntax: 
%    function [T_w_in,T_w_out] = create_T_w(Edge_List,degrees_w)
% 
%   Inputs:
%       Edge_List - 2+size(tau_w) x m matrix of edges
%           First two entries are the edge vertices
%           The next size(tau_w) are the weight sufficient statistics
%       degrees_total - 2 x n matrix of weighted in and out degrees
%
%   Outputs:
%       T_w_in - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x dw_i+ matrix, first entry is parent 
%       T_w_out - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x dw_i- matrix, first entry is child
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/VB_WSBM, Private/CREATE_T_E,

% This is a documentation file for the MEX function.
% The code is in create_T_w.cpp

