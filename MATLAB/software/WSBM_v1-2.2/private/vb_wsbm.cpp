// Christopher Aicher
// vb_wsbm MeX File
//
// 7/25/2013 v1 - It works!
// 8/16/2013 v2 - Added NaN + Edge Distr Values
// 12/1/2013 v3 - Updated for new calling function
//
// Updates mu using VB for the WSBM O(n+m*k^2*t) version
//
// Syntax:
// function [mu_new,flag,diff] = vb_wsbm(mu,T_in,T_out,S_in,S_out,...
//        R,Eta_bra,SEta_bra,maxIter,Tol,verbosity,...
//        degrees_w,dc_flag,mu_0)
//
// Variables-
// Inputs:
//   mu - k x n matrix of vertex label probabilities
//   T_in - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i matrix, first entry is parent 
//   T_out - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i matrix, first entry is child
//   S_in - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i matrix, first entry is parent 
//   S_out - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i matrix, first entry is child
//   R - k x k group interaction bundling structure
//   Eta_bra - r x t matrix of weight bundle expected natural parameter
//   SEta_bra - r x s matrix of edge bundle expected natural parameter
//   maxIter - the maximum number of iterations in the mu loop
//   Tol - min absolute tolerance for convergence 
//   verbosity - how verbose the function is
//         0 -silent, 1 -important only, 2 -more verbose, 3 -most verbose
//   degrees_w - 2 x n matrix of weighted degree (ignores NaNs)
//               First row is in-degree, Second row is out-degree
//   dc_flag - 0 = not degree corrected, 1 = degree-corrected
//   mu_0 - k x n matrix of vertex label prior
//   nanType - 0 = NaN is Missing, 1 = NaN is Non-Edge
//
// Outputs:
//   mu_new - k x n matrix of new vertex label probabilities
//   flag - a flag indicating that mu converged within the tolerance Tol
//   diff - the absolute max difference between iterations of mu
//   
// Other:
//   dT_bra - t x r x k matrix of <T> derivatives
//   dS_bra - s x r x k matrix of <S> derivatives
//   S_last - k x 2 for dc_flag matrix
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
void update_mu(double* mu, const mxArray *prhs[], double* R, double* Eta_bra, 
        double* SEta_bra, int &r, int &n, int &k, int &t, int &s, double* dT_bra, 
        double* dS_bra, double* S_last, double* H, double* mu_new, double verbosity,
        double* mu_0);

void calc_dT_bra(int i, double* mu, const mxArray *prhs[], double* R, int &r, int &n, int &t, int &k, double* dT_bra);
void calc_dS_bra(int i, double* mu, const mxArray *prhs[], double* R, int &r, int &n, int &s, int &k, double* dS_bra, double* S_last);
void calc_S_last(double* mu, const mxArray *prhs[], int &n, int &k, double* S_last);


// Main Function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// Check for errors in inputs/outputs
    if (nrhs != 14) mexErrMsgTxt("Incorrect number of input arguments");
    //if (nlhs != 1) mexErrMsgTxt("Incorrect number of output arguments");
    if (!mxIsCell(prhs[1])) mexErrMsgTxt("Second Arg: T_in needs to be a cell array");
    if (!mxIsCell(prhs[2])) mexErrMsgTxt("Second Arg: T_out needs to be a cell array");
//     mexPrintf("Hello from start\n");
//     mexEvalString("drawnow;");
//   
// Variables
    mxArray *dT_mat, *dS_mat, *S_mat, *H_mat;
    const mwSize *dim_mu, *dim_Eta_cell;
    double *mu, *mu_new, *H, *R, *dT_bra, *dS_bra,
            *S_last, *Eta_bra, *SEta_bra, *mu_0;
    double Tol, verbosity, maxDiff;
    int k, n, r, t, s, maxIter;
    bool convergenceFlag = 0;

// Access Matlab Objects
    mu = mxGetPr(prhs[0]);
    R = mxGetPr(prhs[5]);
    Eta_bra = mxGetPr(prhs[6]); 
    SEta_bra = mxGetPr(prhs[7]);
    maxIter = (int) mxGetScalar(prhs[8]);
    Tol = mxGetScalar(prhs[9]);
    verbosity = mxGetScalar(prhs[10]);
    mu_0 = mxGetPr(prhs[13]);
    
// Figure Out Dimensions
    dim_mu = mxGetDimensions(prhs[0]);
    k = (int) dim_mu[0]; 
    n = (int) dim_mu[1];
    dim_Eta_cell = mxGetDimensions(prhs[6]);
    r = (int) dim_Eta_cell[0];
    t = (int) dim_Eta_cell[1];
    s = (int) mxGetDimensions(prhs[7])[1];
    if(verbosity > 4){
        mexPrintf("k = %u, n = %u, t = %u, s = %u, r = %u, maxIter = %u, Tol = %2.2e\n",k,n,t,s,r,maxIter,Tol);    
        mexEvalString("drawnow;");
    } 
// Allocate Space
    dT_mat = mxCreateDoubleMatrix(r*k*t,1,mxREAL);
    dT_bra = mxGetPr(dT_mat);
    dS_mat = mxCreateDoubleMatrix(r*k*s,1,mxREAL);
    dS_bra = mxGetPr(dS_mat);
    S_mat = mxCreateDoubleMatrix(k*2,1,mxREAL);
    S_last = mxGetPr(S_mat);
    H_mat = mxCreateDoubleMatrix(k,1,mxREAL);
    H = mxGetPr(H_mat);
    plhs[0] = mxCreateDoubleMatrix(k,n,mxREAL);
    mu_new = mxGetPr(plhs[0]);
    
// Function Stuff
    for(int iter = 0; iter < maxIter; iter++){
        if(verbosity > 2){
            mexPrintf("Iteration %u of %u\n",iter, maxIter);    
            mexEvalString("drawnow;");
        }
        update_mu(mu, prhs, R, Eta_bra, SEta_bra, r, n, k, t, s,
                dT_bra, dS_bra, S_last, H, mu_new, verbosity, mu_0);
        maxDiff = 0;
        for(int i = 0; i < n; i++){
            for(int kk = 0; kk < k; kk++){
                double diff = abs(mu[kk+i*k]-mu_new[kk+i*k]);
                mu[kk+i*k] = mu_new[kk+i*k];
                if(diff > maxDiff) maxDiff = diff;
            }
        }
        if(verbosity > 2){
            mexPrintf("Max Difference = %2.2e\n",maxDiff);
            mexEvalString("drawnow;");
        }
        if(maxDiff < Tol){
            convergenceFlag = 1;
            break;
        }
    }
    // Return convergenceFlag
    plhs[1] = mxCreateLogicalScalar(convergenceFlag);
    plhs[2] = mxCreateDoubleScalar(maxDiff);
    
    // Free Memory
    mxDestroyArray(dT_mat);
    mxDestroyArray(dS_mat);
    mxDestroyArray(S_mat);
    mxDestroyArray(H_mat);    
 return;   
}

///////////// Update Mu
void update_mu(double* mu, const mxArray *prhs[], double* R, double* Eta_bra, 
        double* SEta_bra, int &r, int &n, int &k, int &t, int &s, double* dT_bra, 
        double* dS_bra, double* S_last, double* H, double* mu_new, double verbosity,
        double* mu_0){
 // Calculate S_Last
 calc_S_last(mu,prhs,n,k,S_last);  
 for(int i = 0; i < n; i++){
    // Set H to zero
    for(int kk = 0; kk<k; kk++){
        H[kk] = 0;
    }
    // Calculate dT_bra
    if(t > 0){
        calc_dT_bra(i,mu,prhs,R,r,n,t,k,dT_bra);
    }
    for(int kk = 0; kk < k; kk++){
        for(int rr = 0; rr < r; rr++){        
            for(int tt = 0; tt < t; tt++){
                H[kk] += dT_bra[tt+rr*t+kk*r*t]*Eta_bra[rr+tt*r]; 
//                 mexPrintf("Updating kk %u, rr %u, tt = %u\n",kk,rr,tt);
//                 mexEvalString("drawnow;");
            }
        }
    }
    // Calculate dS_bra
    if(s > 0){
        calc_dS_bra(i,mu,prhs,R,r,n,s,k,dS_bra,S_last);
    }
    for(int kk = 0; kk < k; kk++){
        for(int rr = 0; rr < r; rr++){
            for(int ss = 0; ss < s; ss++){
                H[kk] += dS_bra[ss+rr*s+kk*r*s]*SEta_bra[rr+ss*r]; 
//                  if(verbosity > 4){
//                      mexPrintf("dSbra = %2.2f, SEta_bra = %2.2f\n",dS_bra[ss+rr*s+kk*r*s],SEta_bra[rr+ss*r]);
//                      mexPrintf("Updating kk %u, rr %u, ss = %u, H[%u] = %2.2f\n",kk,rr,ss,kk,H[kk]);
//                      mexEvalString("drawnow;");
//                  }
            }
        }
    }
    double maxH = -1000000000000; // This should be -Inf
    // Find max of H
    for(int kk = 0; kk < k; kk++){
        if ((maxH < H[kk]) & (mu_0[kk+i*k] > 0)){
            maxH = H[kk];
        }
        if(verbosity > 4){
            mexPrintf("H[%u] is %f \n",kk,H[kk]);
//             mexPrintf("mu_0[%u] is %f \n",kk,mu_0[kk+i*k]);
            mexEvalString("drawnow;");
        }
    }
    
    if(verbosity > 4){
        mexPrintf("maxH is %f \n",maxH);
        mexEvalString("drawnow;");
    }
    double magH = 0;
    // H = mu_0*exp(H - maxH)
    for(int kk = 0; kk < k; kk++){
        if(mu_0[kk+i*k] > 0){
            H[kk] -= maxH;
            H[kk] = mu_0[kk+i*k]*exp(H[kk]);
            magH += H[kk];
        } else { // mu_0 = 0
            H[kk] = 0;
        }
    }
    if(verbosity > 4){
        mexPrintf("magH of H is %f \n",magH);
        mexEvalString("drawnow;");
    }
    
    // Normalize H and replace mu
    for(int kk = 0; kk < k; kk++){
        mu_new[kk+i*k] = H[kk]/magH;
        if(verbosity > 4){
            mexPrintf("mu_new[%u,%u] is %f \n",kk,i,H[kk]/magH);
            mexEvalString("drawnow;");
//            mexEvalString("keyboard");
        }
    }   
 }
 return; 
}

////////////// Helper Function for calculating d<T>/d\mu
void calc_dT_bra(int i, double* mu, const mxArray *prhs[], double* R, int &r, int &n, int &t, int &k, double* dT_bra){
    // Set dT_bra to zero
    for(int kk = 0; kk < k; kk++){
        for(int rr = 0; rr < r; rr++){    
            for(int tt = 0; tt < t; tt++){
                dT_bra[tt+rr*t+kk*r*t] = 0;
            }
        }
    }
//     mexPrintf("Loading T_out\n");
//     mexEvalString("drawnow;");
    double* T_out = mxGetPr(mxGetCell(prhs[2],i));
    // Sum over d_iW^+
    int d_out = (int) mxGetDimensions(mxGetCell(prhs[2],i))[1];
//     mexPrintf("d_out[%u] = %u\n",i,d_out);
//     mexEvalString("drawnow;");
    for(int e = 0; e < d_out; e++){
        int j = (int) T_out[e*(t+1)]-1;
        for(int k1 = 0; k1 < k; k1++){
            for(int k2 = 0; k2 < k; k2++){
                if( R[k1+k2*k] > 0 ){//R = zeros implies ignore box    
                    for(int tt = 0; tt < t; tt++){
// dT_bra[R(k1,k2),k1] += T[i,j]*mu[j,k2]; 
dT_bra[tt+((int) R[k1+k2*k]-1)*t+k1*r*t] += T_out[e*(t+1)+tt+1]*mu[k2+j*k];
                    }
                }
            }
        }
    }
//     mexPrintf("Loading T_in\n");
//     mexEvalString("drawnow;");
    double* T_in = mxGetPr(mxGetCell(prhs[1],i));
    // Sum over d_iW^-
    int d_in = (int) mxGetDimensions(mxGetCell(prhs[1],i))[1];
//     mexPrintf("d_in[%u] = %u\n",i,d_in);
//     mexEvalString("drawnow;");
    for(int e = 0; e < d_in; e++){
        int j = (int) T_in[e*(t+1)]-1;
        for(int k1 = 0; k1 < k; k1++){
            for(int k2 = 0; k2 < k; k2++){
                if(!mxIsNaN(T_in[e*(t+1)+1])){//NaN implies missing
                    if( R[k2+k1*k] > 0 ){//R = zeros implies ignore box
                        for(int tt = 0; tt < t; tt++){
// dT_bra[R(k2,k1),k1] += T[j,i]*mu[j,k2];
dT_bra[tt+((int) R[k2+k1*k]-1)*t+k1*r*t] += T_in[e*(t+1)+tt+1]*mu[k2+j*k];
                        }
                    }
                }
            }
        }
    }
    return;
}

///////////////// Helper Function for calculating d<S>/d\mu
void calc_dS_bra(int i, double* mu, const mxArray *prhs[], double* R, int &r, int &n, int &s, int &k, double* dS_bra, double* S_last){
    double *degrees_w = mxGetPr(prhs[11]); 
    double *temp_dc;
    temp_dc = new double [2*k];
    
     // Set dS_bra to zero (also temp_dc to zero)
    for(int kk = 0; kk < k; kk++){
        for(int rr = 0; rr < r; rr++){    
            for(int ss = 0; ss < s; ss++){
                dS_bra[ss+rr*s+kk*r*s] = 0;
            }
        }
        temp_dc[kk*2] = 0;
        temp_dc[kk*2+1] = 0;
    }
//     mexPrintf("Loading S_out\n");
//     mexEvalString("drawnow;");
// UPDATE THE DS_BRA via when mu(z) is parent
    double* S_out = mxGetPr(mxGetCell(prhs[4],i));
    int d_out = (int) mxGetDimensions(mxGetCell(prhs[4],i))[1];
//     mexPrintf("d_out[%u] = %u\n",i,d_out);
//     mexEvalString("drawnow;");
    for(int e = 0; e < d_out; e++){
        int j = (int) S_out[e*(s)]-1;
        if(!mxIsNaN(S_out[e*(s)+1])){//NaN implies missing 
            for(int k1 = 0; k1 < k; k1++){
                for(int k2 = 0; k2 < k; k2++){
                    if( R[k1+k2*k] > 0 ){//R = zeros implies ignore box
                        for(int ss = 0; ss < s-1; ss++){
// dS_bra[R(k1,k2),k1] += S[i,j]*mu[j,k2];            
dS_bra[((int) R[k1+k2*k]-1)*s+k1*r*s+ss] += S_out[e*(s)+ss+1]*mu[k2+j*k];
                        }
                    }
                }
            }
        } else { // If missing
if(mxGetScalar(prhs[12])){
    // Degree Corrected
    for(int k2 = 0; k2 < k; k2++){
        temp_dc[k2*2+1] += degrees_w[j*2]*mu[k2+j*k]; // d_in
    }
} else {                       
    // Not Degree Corrected
    for(int k1 = 0; k1 < k; k1++){
        for(int k2 = 0; k2 < k; k2++){
            if( R[k1+k2*k] > 0 ){//R = zeros implies ignore box
dS_bra[((int) R[k1+k2*k]-1)*s+k1*r*s+(s-1)] -= mu[k2+j*k];
            }
        }
    }
}
        }
            
    }
//     mexPrintf("Loading S_in\n");
//     mexEvalString("drawnow;");
// UPDATE THE DS_BRA via when mu(z) is child
    double* S_in = mxGetPr(mxGetCell(prhs[3],i));
    int d_in = (int) mxGetDimensions(mxGetCell(prhs[3],i))[1];
//     mexPrintf("d_in[%u] = %u\n",i,d_in);
//     mexEvalString("drawnow;");
    for(int e = 0; e < d_in; e++){
        int j = (int) S_in[e*(s)]-1;
        if(!mxIsNaN(S_in[e*(s)+1])){//NaN implies missing            
            for(int k1 = 0; k1 < k; k1++){
                for(int k2 = 0; k2 < k; k2++){
                    if( R[k2+k1*k] > 0 ){//R = zeros implies ignore box
                        for(int ss = 0; ss < s-1; ss++){
// dS_bra[R(k2,k1),k1] += mu[j,k2]*S[j,i];
dS_bra[((int) R[k2+k1*k]-1)*s+k1*r*s+ss] += S_in[e*(s)+ss+1]*mu[k2+j*k];
                        }
                    }
                }
            }
        } else { // If missing
if(mxGetScalar(prhs[12])){
    // Degree Corrected
    for(int k2 = 0; k2 < k; k2++){
        temp_dc[k2*2] += degrees_w[j*2+1]*mu[k2+j*k]; // d_out
    }
} else {                       
    // Not Degree Corrected
    for(int k1 = 0; k1 < k; k1++){
        for(int k2 = 0; k2 < k; k2++){
            if( R[k2+k1*k] > 0 ){//R = zeros implies ignore box
dS_bra[((int) R[k2+k1*k]-1)*s+k1*r*s+(s-1)] -= mu[k2+j*k];
            }
        }
    }
}
        }
    }
//Update the final statistic of dS, S = 1;
    if(s > 1){
    if(mxGetScalar(prhs[12])){
        // Degree Corrected
        for(int k1 = 0; k1 < k; k1++){
            for(int k2 = 0; k2 < k; k2++){
                if( R[k2+k1*k] > 0 ){
    //dS_bra[R(k1,k2),k1] += d_wi^+ * (sum_in - missing_dj_in)
    dS_bra[((int) R[k1+k2*k]-1)*s+k1*r*s+(s-1)] += 
            degrees_w[i*2+1]*(S_last[k2+k]-temp_dc[k2*2+1]);
    //dS_bra[R(k2,k1),k1] += d_wi^- * (sum_out - missing_dj_out) 
    dS_bra[((int) R[k2+k1*k]-1)*s+k1*r*s+(s-1)] += 
            degrees_w[i*2]*(S_last[k2]-temp_dc[k2*2]);
                }
            }
        } 
    } else {
       // Not Degree Corrected
       for(int k1 = 0; k1 < k; k1++){
            for(int k2 = 0; k2 < k; k2++){
                if( R[k2+k1*k] > 0 ){
    dS_bra[((int) R[k1+k2*k]-1)*s+k1*r*s+(s-1)] += S_last[k2+k];
    dS_bra[((int) R[k2+k1*k]-1)*s+k1*r*s+(s-1)] += S_last[k2];   
                }
            }
        }  
    }
    }
    delete [] temp_dc;
    return;
}

//////////////// Helper Function for Calculating S_last
void calc_S_last(double* mu, const mxArray *prhs[], int &n, int &k, double* S_last){
    double *degrees_w = mxGetPr(prhs[11]);
    // Set S_last to zero
    for(int k1 = 0; k1 < k; k1++){
        S_last[k1] = 0;
        S_last[k1+k] = 0;
    }
    if(mxGetScalar(prhs[12])){
        // Degree Corrected
        for(int jj = 0; jj < n; jj++){
           for(int k1 = 0; k1 < k; k1++){        
    S_last[k1] += degrees_w[jj*2+1]*mu[k1+jj*k]; // d_out(jj)*mu(k1,jj)
    S_last[k1+k] += degrees_w[jj*2]*mu[k1+jj*k]; // d_in(jj)*mu(k1,jj)   
            }
        }
//         if(mxGetScalar(prhs[10]) > 4){
//             for(int kk = 0; kk < k*2; kk++){
//                 mexPrintf("S_last[%u] = %f \n",kk,S_last[kk]);
//             }
//             mexEvalString("drawnow;");
//         }
    } else {
        // Not Degree Corrected
        for(int jj = 0; jj < n; jj++){
            for(int k1 = 0; k1 < k; k1++){
    S_last[k1] += mu[k1+jj*k]; // mu(k1,jj)
    S_last[k1+k] += mu[k1+jj*k]; // mu(k1,jj)
            }
        } 
    }
    
    return;
}

// The End