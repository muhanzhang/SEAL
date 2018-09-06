function [newT_w_in,newT_w_out,newT_e_in,newT_e_out] = ...
    create_T_bp(T_w_in,T_w_out,T_e_in,T_e_out)
%CREATE_T_BP formats the Adjancency List to a BP friendly format
%
%   CREATE_T_BP is a helper function for wsbm.m
%       -Not exactly efficient, but faster than MATLAB
%
%   Syntax: 
%   function [newT_w_in,newT_w_out,newT_e_in,newT_e_out] = ...
%       create_T_bp(T_w_in,T_w_out,T_e_in,T_e_out)
% 
%   Inputs:
%       T_w_in - n x 1 cell array of weighted statistics
%           Each cell is a 1+t x d_i(in) matrix, first entry is parent
%           Sorted in increasing order
%       T_w_out - n x 1 cell array of weighted statistics
%           Each cell is a 1+t x d_i(out) matrix, first entry is child
%           Sorted in increasing order
%       T_e_in - n x 1 cell array of edge statistics
%           Each cell is a 1+s-1 x d_i(in) matrix, first entry is parent
%           Sorted in increasing order
%       T_e_out - n x 1 cell array of edge statistics
%           Each cell is a 1+s-1 x d_i(out) matrix, first entry is child
%           Sorted in increasing order
%
%   Outputs:
%     !!!Note: T_w and T_e are different in BP than in VB!!!
%       T_w_in - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x d_i*+ matrix, first entry is parent 
%       T_w_out - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x d_i*- matrix, first entry is child
%       T_e_in - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i*+ matrix, first entry is parent 
%       T_e_out - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i*- matrix, first entry is child
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/BP_WSBM,
%   Private/CREATE_T_W, Private/CREATE_T_E

% This is a documentation file for the MEX function.
% The code is in create_T_bp.cpp

