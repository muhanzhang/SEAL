// Christopher Aicher
// calc_T_e_bra MeX File
// 7/14/2013
// Updated: 12/1/2013
//
// Updates mu using VB for the WSBM O(m*k^2*t) version
//
// Syntax:
// function [S_bra] = calc_T_e_bra(mu,S_out,R,r,degrees_w,dc_flag)
//
// Variables-
// Inputs:
//   mu - k x n matrix of vertex label probabilities
//   S_out - n x 1 cell array of S lists
//      Each cell is a 1+s-1 x d_i matrix, first entry is child
//   R - k x k group interaction bundling structure
//   r - the number of different edge bundles 
//   degrees_w - 2 x n matrix of weighted degree (ignores NaNs)
//               First row is in-degree, Second row is out-degree
//   dc_flag - 0 = not degree corrected, 1 = degree-corrected
//   nanType - 0 = NaN is Missing, 1 = NaN is Non-Edge
//
// Outputs: 
//   S_bra - s x r matrix of the sufficient statistic sum of each edge bundle
// Other:
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
    if (nrhs != 6) mexErrMsgTxt("Incorrect number of input arguments");
    //if (nlhs != 1) mexErrMsgTxt("Incorrect number of output arguments");
    if (!mxIsCell(prhs[1])) mexErrMsgTxt("Second Arg: S_out needs to be a cell array");
    //mexPrintf("Hello from start\n");
    //mexEvalString("drawnow;");
  
// Variables
    const mwSize *dim_mu;
    double *mu, *R, *S_bra, *S_out, *degrees_w;
    int k, n, r, s;
    bool convergenceFlag = 0;

// Access Matlab Objects
    mu = mxGetPr(prhs[0]);
    R = mxGetPr(prhs[2]);
    r = (int) mxGetScalar(prhs[3]); 
    degrees_w = mxGetPr(prhs[4]);
     
// Figure Out Dimensions
    dim_mu = mxGetDimensions(prhs[0]);
    k = (int) dim_mu[0]; 
    n = (int) dim_mu[1];
    s = (int) mxGetDimensions(mxGetCell(prhs[1],0))[0]-1+1;
//     mexPrintf("k = %u, n = %u, t = %u, r = %u\n",k,n,t,r);    
//     mexEvalString("drawnow;");
//     mexEvalString("keyboard");
    
// Allocate Space
    plhs[0] = mxCreateDoubleMatrix(s,r,mxREAL);        
    S_bra = mxGetPr(plhs[0]);
    //mexPrintf("Stuff was allocated\n");
    //mexEvalString("drawnow;");
    
// Function Stuff
    for(int i = 0; i < n; i++){
        int d_out = (int) mxGetDimensions(mxGetCell(prhs[1],i))[1];
        S_out = mxGetPr(mxGetCell(prhs[1],i));
        //mexPrintf("Got Cell i = %u\n",i);
        //mexEvalString("drawnow;");
        for(int e = 0; e < d_out; e++){
            int j = (int) S_out[e*(s)]-1;
            if(!mxIsNaN(S_out[e*(s)+1])){ // Ignore NaN = Missing Entries
            for(int k1 = 0; k1 < k; k1++){
                for(int k2 = 0; k2 < k; k2++){
                    if( R[k1+k2*k] > 0 ){
                        if( R[k1+k2*k] > r){
                            mexErrMsgTxt("Misspecification: R[k1,k2] > r\n");
                        }
                        for(int ss = 0; ss < s-1; ss++){
// S_bra(R(k1,k2)) = mu(i,k1)mu(j,k2)S(i,j)
//mexPrintf("Calculating S_bra[%u,%u] = mu[%u,%u] mu[%u,%u] S_out[%u,%u]\n",(int)(R[k1+k2*k]-1),tt,k1,i,k2,j,e,tt+1);
//mexEvalString("drawnow;");
S_bra[((int)(R[k1+k2*k]-1))*s+ss] += mu[k1+i*k]*mu[k2+j*k]*S_out[e*(s)+ss+1];
                        }                   
                    }
                }
            }
            } else {
            for(int k1 = 0; k1 < k; k1++){
                for(int k2 = 0; k2 < k; k2++){
                    if( R[k1+k2*k] > 0 ){
                        if( R[k1+k2*k] > r){
                            mexErrMsgTxt("Misspecification: R[k1,k2] > r\n");
                        }    
                        if(mxGetScalar(prhs[5])){
// Degree Corrected Case
S_bra[((int)(R[k1+k2*k]-1))*s+s-1] -= degrees_w[2*i+1]*mu[k1+i*k]*
                                        degrees_w[2*j]*mu[k2+j*k];
                        } else { 
// Non-Degree Corrected Case
S_bra[((int)(R[k1+k2*k]-1))*s+s-1] -= mu[k1+i*k]*mu[k2+j*k];
                        }
                    }
                }
            }
            }
        }
    }
 return;   
}







