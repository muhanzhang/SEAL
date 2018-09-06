function [T_e_bra] = calc_T_e_bra(mu,T_e_out,R,r,degrees_w,dc_flag)
%CALC_T_E_BRA calculates the expected sufficient statistic for the edges 
% in each edge bundle.
%
%   CALC_T_E_BRA is a helper function for main_alg.m
%
%   Syntax: 
%       function [T_e_bra] = calc_T_e_bra(mu,T_w_out,R,r,degrees_w,dc_flag)
%
%   Inputs:
%       mu - k x n matrix of vertex label probabilities
%       T_e_out - n x 1 cell array of T_e lists
%           Each cell is a 1+t_e-1 x d_i matrix, first entry is child
%       R - k x k group interaction bundling structure
%       r - the number of different edge bundles 
%       degrees_w - 2 x n matrix of weighted degree (ignoring missing)
%               First row is in-degree, Second row is out-degree
%       dc_flag - 0 = not degree corrected, 1 = degree-corrected.
%
%   Outputs:
%       T_e_bra - t_w x r matrix of weighted sufficient statistics
%           The last row needs to be updated to account for non-edges
%           See MAIN_ALG.m's code for more.
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/BP_WSBM,
%   Private/CREATE_T_E, Private/CALC_T_W_BRA

% This is a documentation file for the MEX function.
% The code is in calc_T_w_bra.cpp

