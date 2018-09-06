function [T_w_bra] = calc_T_w_bra(mu,T_w_out,R,r)
%CALC_T_W_BRA calculates the expected sufficient statistic for the weights 
% in each edge bundle.
%
%   CALC_T_W_BRA is a helper function for main_alg.m
%
%   Syntax: 
%       function [T_w_bra] = calc_T_w_bra(mu,T_w_out,R,r)
%
%   Inputs:
%       mu - k x n matrix of vertex label probabilities
%       T_w_out - n x 1 cell array of T_w lists
%           Each cell is a 1+t_w x d_i matrix, first entry is child
%       R - k x k group interaction bundling structure
%       r - the number of different edge bundles 
%
%   Outputs:
%       T_w_bra - t_w x r matrix of weighted sufficient statistics
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/BP_WSBM,
%   Private/CREATE_T_W, Private/CALC_T_E_BRA

% This is a documentation file for the MEX function.
% The code is in calc_T_w_bra.cpp

