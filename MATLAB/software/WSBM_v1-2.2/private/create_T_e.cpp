// Christopher Aicher
// MeX File
// 7/14/2013
//
// Create Edge Lists S_in S_out
//
// Syntax:
// function [S_in,S_out] = create_T_e(EdgeList,degrees)
//
// Variables-
// Inputs:
//    EdgeList - 2+t x m matrix
//         First two entries are the edge vertices 
//         the next t are the sufficient statistic of the Weight
//    degrees - 2 x n matrix of in and out degrees
//
// Outputs: 
//   S_in - n x 1 cell array 
//      Each cell is a 1+t x d_i(in) matrix, first entry is parent 
//   S_out - n x 1 cell array
//      Each cell is a 1+t x d_i(out) matrix, first entry is child 
//   (i.e. S_in{i}[s,1] = the parent of the s-th edge into i
//         S_in{i}[s,t+1] = the t-th statistic of the s-th edge of i)
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
    if (nrhs != 2) mexErrMsgTxt("Incorrect number of input arguments");
    if (nlhs != 2) mexErrMsgTxt("Incorrect number of output arguments");
    
    //mexPrintf("Hello from start\n");
    //mexEvalString("drawnow;");
  
// Variables
    mxArray *S_in_cell, *S_out_cell,*S_in_temp, *S_out_temp, *degree_mat;
    const mwSize *dim_Edge_list;
    double *Edge_list, *degree, *S_in, *S_out;
    double Sol;
    int n, m, t, *degree_count;

// Access Matlab Objects
    Edge_list = mxGetPr(prhs[0]);
    degree = mxGetPr(prhs[1]);
    
// Figure Out Dimensions
    n = (int) mxGetDimensions(prhs[1])[1];
    dim_Edge_list = mxGetDimensions(prhs[0]);
    t = (int) dim_Edge_list[0] - 2; 
    m = (int) dim_Edge_list[1];
    //mexPrintf("n = %u, m = %u, t = %u,\n",n,m,t);    
    //mexEvalString("drawnow;");
    
// Allocate Space
    S_in_cell = plhs[0] = mxCreateCellMatrix(n,1);
    S_out_cell = plhs[1] = mxCreateCellMatrix(n,1);
    S_in_temp = mxCreateCellMatrix(n,1);
    S_out_temp = mxCreateCellMatrix(n,1);
    for(int i = 0; i < n; i++){
        //mexPrintf("For i = %u, S_in size %u by %u, S_out size %u by %u\n",i,1+t,(int)degree[i*2],1+t,(int)degree[i*2+1]);
        //mexEvalString("drawnow;");
        mxSetCell(S_in_temp,i,mxCreateDoubleMatrix(1+t,(int) degree[i*2],mxREAL));
        mxSetCell(S_out_temp,i,mxCreateDoubleMatrix(1+t,(int) degree[i*2+1],mxREAL));
    }
    degree_mat = mxCreateNumericMatrix(2,n,mxINT64_CLASS,mxREAL);
    degree_count = (int*) mxGetPr(degree_mat);
    //mexPrintf("Stuff was allocated\n",n,m,t);    
    //mexEvalString("drawnow;");
    
// Function Stuff
    // Fill in S_temp
    for(int e = 0; e < m; e++){
        //mexPrintf("Current Edge %u\n",e);    
        //mexEvalString("drawnow;");
        int i = (int) Edge_list[e*(2+t)]-1;
        int j = (int) Edge_list[e*(2+t)+1]-1;
        //mexPrintf("Current i = %u, j = %u \n",i,j);    
        //mexEvalString("drawnow;");
        S_in = mxGetPr(mxGetCell(S_in_temp,j));
        //mexPrintf("Current S_in[%u]\n",degree_count[j*2]*(1+t));
        //mexEvalString("drawnow;");
        S_in[degree_count[j*2]*(1+t)] = (double) i+1;
        S_out = mxGetPr(mxGetCell(S_out_temp,i));
        //mexPrintf("Current S_out[%u]\n",degree_count[i*2+1]*(1+t));
        //mexEvalString("drawnow;");
        S_out[degree_count[i*2+1]*(1+t)] = (double) j+1;
        for(int tt = 0; tt < t; tt++){
            //mexPrintf("Updating e = %u for t = %u\n",e,tt);    
            //mexEvalString("drawnow;");
            S_in[degree_count[j*2]*(1+t)+tt+1] = Edge_list[e*(2+t)+tt+2];
            S_out[degree_count[i*2+1]*(1+t)+tt+1] = Edge_list[e*(2+t)+tt+2];
        }
        degree_count[j*2]++;
        degree_count[i*2+1]++;
    }
    // Compress Space
    // Fill in S_in_cell
    for(int i = 0; i < n; i++){
        int d_tot = 0; // Number of Neighbors
        double cur = 0;
        double* S_temp = mxGetPr(mxGetCell(S_in_temp,i));
        int d = (int) mxGetDimensions(mxGetCell(S_in_temp,i))[1];
        // Count the number of Neighbors
        for(int dd = 0; dd < d; dd++){
            if(S_temp[dd*(1+t)] > cur){
                d_tot++;
                cur = S_temp[dd*(1+t)];
            }            
        }
        // Allocate Space for Output
        mxSetCell(S_in_cell,i,mxCreateDoubleMatrix(1+t,d_tot,mxREAL));
        double* S_in = mxGetPr(mxGetCell(S_in_cell,i));
        cur = 0;
        int d_cur = -1;
        for(int dd = 0; dd < d; dd++){
            if(S_temp[dd*(1+t)] > cur){
                d_cur++;
                cur = S_temp[dd*(1+t)];
                S_in[d_cur*(1+t)] = cur;
            }
//             mexPrintf("Current d %u of %u\n",d_cur,d_tot);
//             mexEvalString("drawnow;");
            for(int tt = 0; tt < t; tt++){
                S_in[d_cur*(1+t)+tt+1] += S_temp[dd*(1+t)+tt+1];
            }
        }
    }
    // Fill in S_out_cell
    for(int i = 0; i < n; i++){
        int d_tot = 0; // Number of Neighbors
        double cur = 0;
        double* S_temp = mxGetPr(mxGetCell(S_out_temp,i));
        int d = (int) mxGetDimensions(mxGetCell(S_out_temp,i))[1];
        // Count the number of Neighbors
        for(int dd = 0; dd < d; dd++){
            if(S_temp[dd*(1+t)] > cur){
                d_tot++;
                cur = S_temp[dd*(1+t)];
            }            
        }
        // Allocate Space for Output
        mxSetCell(S_out_cell,i,mxCreateDoubleMatrix(1+t,d_tot,mxREAL));
        double* S_out = mxGetPr(mxGetCell(S_out_cell,i));
        cur = 0;
        int d_cur = -1;
        for(int dd = 0; dd < d; dd++){
            if(S_temp[dd*(1+t)] > cur){
                d_cur++;
                cur = S_temp[dd*(1+t)];
                S_out[d_cur*(1+t)] = cur;
            }
            for(int tt = 0; tt < t; tt++){
                S_out[d_cur*(1+t)+tt+1] += S_temp[dd*(1+t)+tt+1];
            }
        }
    }          
//  Deallocate Memory
//     mexPrintf("Freeing Allocated Memory\n");
//     mexEvalString("drawnow;");
    mxDestroyArray(S_in_temp);
    mxDestroyArray(S_out_temp);
    mxDestroyArray(degree_mat);
 return;   
}



