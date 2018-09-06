// Christopher Aicher
// Initialize Messages MeX File
// 12/20/2013
//
// Create T_w_in, T_w_out, T_e_in, T_e_out for BP
//
// Syntax:
// function [newT_w_in,newT_w_out,newT_e_in,newT_e_out] = 
//                              create_T_bp(T_w_in,T_w_out,T_e_in,T_e_out)
//
// Variables-
// Inputs:
//   T_w_in - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i(in) matrix, first entry is parent
//      Sorted in increasing order
//   T_w_out - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i(out) matrix, first entry is child
//      Sorted in increasing order
//   T_e_in - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i(in) matrix, first entry is parent
//      Sorted in increasing order
//   T_e_out - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i(out) matrix, first entry is child
//      Sorted in increasing order
//
// Outputs: 
//   mes - n x 1 cell array of messages
//      Each cell is a 1+k x d_i* matrix, first entry is neighbor 
//
//
//  Copyright 2013-2014 Christopher Aicher
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#include <math.h>
#include <matrix.h>
#include <mex.h>

// Prototypes


// Main Function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// Check for errors in inputs/outputs
    if (nrhs != 4) mexErrMsgTxt("Incorrect number of input arguments");
    if (nlhs != 4) mexErrMsgTxt("Incorrect number of output arguments");
    
// Variables
    mxArray *newT_w_in, *newT_w_out, *newT_e_in, *newT_e_out;
    int n, t, s;
    
// Figure Out Dimensions
    if(!mxIsEmpty(prhs[0])){
        t = (int) mxGetDimensions(mxGetCell(prhs[0],0))[0]-1;
    } else {
        t = 0;
    }
    if(!mxIsEmpty(prhs[2])){
        s = (int) mxGetDimensions(mxGetCell(prhs[2],0))[0];
    } else {
        s = 0;
    }
    if (s <= 0) mexErrMsgTxt("E_Distr statistics must not be empty");
    n = (int) mxGetDimensions(prhs[2])[0];
//     mexPrintf("n = %u, t = %u, s = %u\n",n,t,s);    
//     mexEvalString("drawnow;");
    
// Allocate Space for Outputs
    newT_w_in = plhs[0] = mxCreateCellMatrix(n,1);
    newT_w_out = plhs[1] = mxCreateCellMatrix(n,1);
    newT_e_in = plhs[2] = mxCreateCellMatrix(n,1);
    newT_e_out = plhs[3] = mxCreateCellMatrix(n,1);
    
// Function Stuff
    // Find d_i* = #unique neighbors
    for(int i = 0; i < n; i++){
        double* S_in = mxGetPr(mxGetCell(prhs[2],i));
        double* S_out = mxGetPr(mxGetCell(prhs[3],i));
        int d_in = mxGetDimensions(mxGetCell(prhs[2],i))[1];
        int d_out = mxGetDimensions(mxGetCell(prhs[3],i))[1];
        // Find d_i* = # unique neighbors
        int s_in = 0; 
        int s_out = 0;
        int d_i_star = 0;
        while((s_in < d_in) && (s_out < d_out)){
            if(S_in[s_in*(1+s-1)] > S_out[s_out*(1+s-1)]){
                s_out++; d_i_star++;
            }else if(S_in[s_in*(1+s-1)] < S_out[s_out*(1+s-1)]){
                s_in++; d_i_star++;
            }else{
                s_in++;
            }
        }
        d_i_star += (d_in-s_in)+(d_out-s_out);
//      Allocate Space for newT_cells
        mxSetCell(newT_w_in,i,mxCreateDoubleMatrix(1+t,d_i_star,mxREAL));
        mxSetCell(newT_w_out,i,mxCreateDoubleMatrix(1+t,d_i_star,mxREAL));
        mxSetCell(newT_e_in,i,mxCreateDoubleMatrix(1+s-1,d_i_star,mxREAL));
        mxSetCell(newT_e_out,i,mxCreateDoubleMatrix(1+s-1,d_i_star,mxREAL));      
        
// Update Each Cell        
        double* temp_t_in = mxGetPr(mxGetCell(newT_w_in,i));
        double* temp_t_out = mxGetPr(mxGetCell(newT_w_out,i));
        double* temp_s_in = mxGetPr(mxGetCell(newT_e_in,i));
        double* temp_s_out = mxGetPr(mxGetCell(newT_e_out,i));
        // Label temp 
        int d_cur = 0;
        s_in = 0; 
        s_out = 0;
        while((s_in < d_in) && (s_out < d_out)){
            if(S_in[s_in*(1+s-1)] > S_out[s_out*(1+s-1)]){
                temp_t_in[d_cur*(1+t)] = S_out[s_out*(1+s-1)];
                temp_t_out[d_cur*(1+t)] = S_out[s_out*(1+s-1)];
                temp_s_in[d_cur*(1+s-1)] = S_out[s_out*(1+s-1)];
                temp_s_out[d_cur*(1+s-1)] = S_out[s_out*(1+s-1)];
                s_out++; d_cur++;
            }else if(S_in[s_in*(1+s-1)] < S_out[s_out*(1+s-1)]){
                temp_t_in[d_cur*(1+t)] = S_in[s_in*(1+s-1)];
                temp_t_out[d_cur*(1+t)] = S_in[s_in*(1+s-1)];
                temp_s_in[d_cur*(1+s-1)] = S_in[s_in*(1+s-1)];
                temp_s_out[d_cur*(1+s-1)] = S_in[s_in*(1+s-1)];
                s_in++; d_cur++;
            }else{
                s_in++;
            }           
        } 
        while(s_in < d_in){
            temp_t_in[d_cur*(1+t)] = S_in[s_in*(1+s-1)];
            temp_t_out[d_cur*(1+t)] = S_in[s_in*(1+s-1)];
            temp_s_in[d_cur*(1+s-1)] = S_in[s_in*(1+s-1)];
            temp_s_out[d_cur*(1+s-1)] = S_in[s_in*(1+s-1)];        
            s_in++; d_cur++;
        }
        while(s_out < d_out){
            temp_t_in[d_cur*(1+t)] = S_out[s_out*(1+s-1)];
            temp_t_out[d_cur*(1+t)] = S_out[s_out*(1+s-1)];
            temp_s_in[d_cur*(1+s-1)] = S_out[s_out*(1+s-1)];
            temp_s_out[d_cur*(1+s-1)] = S_out[s_out*(1+s-1)];
            s_out++; d_cur++;
        }
        
        // Fill in temp values
        s_in = 0;
        s_out = 0;
        int t_in = 0;
        int t_out = 0;
        for(d_cur = 0;  d_cur < d_i_star; d_cur++){
            // Update temp_t_in
        if(t > 0){
            double* T_in = mxGetPr(mxGetCell(prhs[0],i));
            double* T_out = mxGetPr(mxGetCell(prhs[1],i));
            if(temp_t_in[d_cur*(1+t)] == T_in[t_in*(1+t)]){
                for(int tt = 0; tt < t; tt++){
        temp_t_in[d_cur*(1+t)+tt+1] = T_in[t_in*(1+t)+tt+1]; 
                }
                t_in++;
            } else {
                for(int tt = 0; tt < t; tt++){
        temp_t_in[d_cur*(1+t)+tt+1] = mxGetNaN();  
                }
            }
            // Update temp_t_out
            if(temp_t_out[d_cur*(1+t)] == T_out[t_out*(1+t)]){
                for(int tt = 0; tt < t; tt++){
        temp_t_out[d_cur*(1+t)+tt+1] = T_out[t_out*(1+t)+tt+1]; 
                }
                t_out++;
            } else {
                for(int tt = 0; tt < t; tt++){
        temp_t_out[d_cur*(1+t)+tt+1] = mxGetNaN();  
                }
            }
        }
            // Update temp_s_in
            if(temp_s_in[d_cur*(1+s-1)] == S_in[s_in*(1+s-1)]){
                for(int ss = 0; ss < s-1; ss++){
        temp_s_in[d_cur*(1+s-1)+ss+1] = S_in[s_in*(1+s-1)+ss+1]; 
                }
                s_in++;
            } else {
                for(int ss = 0; ss < s-1; ss++){
        temp_s_in[d_cur*(1+s-1)+ss+1] = 0;  
                }
            }
            // Update temp_s_out
            if(temp_s_out[d_cur*(1+s-1)] == S_out[s_out*(1+s-1)]){
                for(int ss = 0; ss < s-1; ss++){
        temp_s_out[d_cur*(1+s-1)+ss+1] = S_out[s_out*(1+s-1)+ss+1]; 
                }
                s_out++;
            } else {
                for(int ss = 0; ss < s-1; ss++){
        temp_s_out[d_cur*(1+s-1)+ss+1] = 0;  
                }
            }
            
        }
    }
 return;   
}