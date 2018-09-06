function [mu_new,flag,diff] = vb_wsbm(mu,T_w_in,T_w_out,T_e_in,T_e_out,...
    R,Eta_w_bra,Eta_e_bra,maxIter,Tol,verbosity,degrees_w,dc_flag,mu_0)
%VB_WSBM updates the vertex-label posteriors using a naive-VB 
%   VB_WSBM is a helper function to wsbm.m specifically used by main_alg.m
%
%   Syntax: 
%    function [mu_new,flag,diff] = vb_wsbm(mu,T_w_in,T_w_out,T_e_in,T_e_out,...
%       R,Eta_w_bra,Eta_e_bra,maxIter,Tol,verbosity,degrees_w,dc_flag)
% 
%   Inputs:
%       mu - k x n matrix of vertex label probabilities
%       T_w_in - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x dw_i+ matrix, first entry is parent 
%       T_w_out - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x dw_i- matrix, first entry is child
%       T_e_in - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i+ matrix, first entry is parent 
%       T_e_out - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i- matrix, first entry is child
%       R - k x k group interaction bundling structure
%       Eta_w_bra - r x t matrix of weight bundle expected natural parameter
%       Eta_e_bra - r x s matrix of edge bundle expected natural parameter
%       maxIter - the maximum number of iterations in the mu loop
%       Tol - min absolute tolerance for convergence 
%       verbosity - how verbose the function is
%         0 -silent, 1 -important only, 2 -more verbose, 3 -most verbose
%       degrees_w - 2 x n matrix of weighted degree (ignores NaNs)
%                   First row is in-degree, Second row is out-degree
%       dc_flag - 0 = not degree corrected, 1 = degree-corrected
%       mu_0 - k x n matrix of vertex label prior
%
%   Outputs:
%       mu_new - k x n matrix of new vertex label probabilities
%       flag - a flag indicating that mu converged within the tolerance
%       diff - the absolute max difference between iterations of mu
%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private/MAIN_ALG, Private/BP_WSBM,
%       Private/CREATE_T_W, Private/CREATE_T_E

% This is a documentation file for the MEX function.
% The code is in vb_wsbm.cpp
