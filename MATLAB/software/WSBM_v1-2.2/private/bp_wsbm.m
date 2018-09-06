function [mu_new,flag,diff] = bp_wsbm(mu,mes,T_w_in,T_w_out,T_e_in,T_e_out,...
    R,Eta_w_bra,Eta_e_bra,maxIter,Tol,verbosity,degrees_w,dc_flag,mu_0)
%BP_WSBM updates the vertex-label posteriors using belief propagation 
%   BP_WSBM is a helper function to wsbm.m specifically used by main_alg.m
%
%   Syntax: 
%     function [mu_new,flag,diff] = bp_wsbm(mu,mes,T_w_in,T_w_out,...
%       T_e_in,T_e_out,R,Eta_w_bra,Eta_e_bra,maxIter,Tol,verbosity,...
%       degrees_w,dc_flag)
% 
%   Inputs:
%       mu - k x n matrix of vertex label probabilities
%       mes - n x 1 cell array of messages
%          Each cell is a 1+k x d_i* matrix, first entry is neighbor 
%          mes{ii}[jj*d_i*] is the belief about jj in the absense of ii
%     !!!Note: T_w and T_e are different in BP than in VB!!!
%       T_w_in - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x d_i*+ matrix, first entry is parent 
%       T_w_out - n x 1 cell array of weighted statistics
%          Each cell is a 1+t x d_i*- matrix, first entry is child
%       T_e_in - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i*+ matrix, first entry is parent 
%       T_e_out - n x 1 cell array of edge statistics
%          Each cell is a 1+s-1 x d_i*- matrix, first entry is child
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
%       diff - the absolute max difference between iterations of mu%
%   See the additional notes folder for more information on how 
%    the algorithm works. The code for the MEX file is in the 
%    associate .cpp file in C++.     
%
%   See Also Private\MAIN_ALG, Private\VB_WSBM, Private\CREATE_T_BP

% This is a documentation file for the MEX function.
% The code is in bp_wsbm.cpp
