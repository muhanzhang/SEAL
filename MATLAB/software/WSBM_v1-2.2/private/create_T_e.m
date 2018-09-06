function [T_e_in,T_e_out] = create_T_e(Edge_List,degrees_total)
%CREATE_T_E formats the Edge_List into an Adjancency List of 
% Edge Distribution Statistics 
%
%   CREATE_T_E is a helper function for wsbm.m
%
%   Syntax: 
%    function [T_e_in,T_e_out] = create_T_e(Edge_List,degrees_total)
% 
%   Inputs:
%       Edge_List - 2+size(tau_e) x m matrix of edges
%           First two entries are the edge vertices
%           The next size(tau_e) are the edge sufficient statistics
%       degrees_total - 2 x n matrix of total in and out degrees
%
%   Outputs:
%       T_e_in - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i+ matrix, first entry is parent 
%       T_e_out - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i- matrix, first entry is child
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/VB_WSBM, Private/CREATE_T_W,

% This is a documentation file for the MEX function.
% The code is in create_T_e.cpp

