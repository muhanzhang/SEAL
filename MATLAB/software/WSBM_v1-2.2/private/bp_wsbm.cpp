// Christopher Aicher
// BP_Update MeX File
// 12/20/2013
//
// Updates mu using BP for the WSBM O(d_w*(n+m)*k^2*t) version
//
// Syntax:
// function [mu_new,flag,diff] = bp_wsbm(mu,mes,T_in,T_out,S_in,S_out,R,...
//                      Eta_bra,SEta_bra,maxIter,Tol,verbosity,mu_0)
//
// Variables-
// Inputs:
//   mu - k x n matrix of vertex label probabilities
//   mes - n x 1 cell array of messages
//      Each cell is a 1+k x d_i* matrix, first entry is neighbor 
//      mes{ii}[jj*d_i*] is the belief about jj in the absense of ii
//   T_in - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i* matrix, first entry is parent
//      Sorted in increasing order
//   T_out - n x 1 cell array of weighted statistics
//      Each cell is a 1+t x d_i* matrix, first entry is child
//      Sorted in increasing order
//   S_in - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i* matrix, first entry is parent
//      Sorted in increasing order
//   S_out - n x 1 cell array of edge statistics
//      Each cell is a 1+s-1 x d_i* matrix, first entry is child
//      Sorted in increasing order
//   R - k x k group interaction bundling structure
//   Eta_bra - r x t matrix of weight bundle expected natural parameter
//   SEta_bra - r x s matrix of edge bundle expected natural parameter
//   maxIter - the maximum number of iterations in the mu loop
//   Tol - min absolute tolerance for convergence 
//   verbosity - how verbose the function is
//         0 -silent, 1 -important only, 2 -more verbose, 3 -most verbose
//   mu_0 - k x n vector of vertex label prior probabilities
//
// Outputs:  
//   mu_new - k x n matrix of new vertex label probabilities
//   flag - boolean indicating proper convergence (1) or (0) if maxIter exceeded
//   diff - abs max difference between current and prev iteration of mu
//   
// Other:
//   M_cell - n x 1 cell array of pseudo-evidence values
//      Each cell is a k^2+1 x d_i* matrix, last entry is neighbor
//   logM_null - k x n matrix of non-edge approximate pseudo-evidence values
//
// PseudoCode Outline:
//    Initialize mes
//    Calculate M_cell for each weighted edge
//    Loop:
//       Approximate logM_null for all non-weighted edges
//       Update mes (for weighted edges)
//       Update mu  (for all vertices)
//       Check Convergence
//    End Loop
//    Return Result
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
void calc_M_cell(const mxArray *prhs[], mxArray *M_cell, 
        double* R, double* Eta_bra, double* SEta_bra, 
        int &r, int &n, int &k, int &t, int &s);
// void update_M_ij(double* M_ij,double* T, double* R, double* Eta, 
//         int &k, int &r, int &t, int &d_cur, int &t_cur, bool IsOut);
// void update_M_ijS(double* M_ij, double* S, double* R, double* SEta, 
//         int &k, int &r, int &s, int &d_cur, int &s_cur, bool IsOut);
void update_mu(double* mu, const mxArray *prhs[], mxArray *M_cell,
        int &n, int &k, double* mu_new, double verbosity, 
        double *logM_null, double * mu_0);
void update_mes(double* mu, const mxArray *prhs[], mxArray *M_cell,
        int &r, int &n, int &k, int &t, int &s, double *logM_null, 
        double * mu_0);
void update_logM_null(double *mu, const mxArray *prhs[], double *logM_null, 
        double *R, double* SEta_bra, int &r, int &n, int &k, int &s);

// Main Function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// Check for errors in inputs/outputs
    if (nrhs != 13) mexErrMsgTxt("Incorrect number of input arguments");
    if (nlhs != 3) mexErrMsgTxt("Incorrect number of output arguments");
    if (!mxIsCell(prhs[1])) mexErrMsgTxt("Arg: Mes needs to be a cell array");
    if (!mxIsCell(prhs[2])) mexErrMsgTxt("Arg: T_in needs to be a cell array");
    if (!mxIsCell(prhs[3])) mexErrMsgTxt("Arg: T_out needs to be a cell array");
    if (!mxIsCell(prhs[4])) mexErrMsgTxt("Arg: S_in needs to be a cell array");
    if (!mxIsCell(prhs[5])) mexErrMsgTxt("Arg: S_out needs to be a cell array");
//     mexPrintf("Hello from start\n");
//     mexEvalString("drawnow;");
//   
// Variables
    mxArray *M_cell, *logM_null_mat;
    double *mu, *mu_new, *R, *Eta_bra, *SEta_bra, *logM_null, *mu_0;
    double Tol, verbosity, maxDiff;
    int k, n, r, t, s, maxIter;
    bool convergenceFlag = 0;

// Access Matlab Objects
    mu = mxGetPr(prhs[0]);
    R = mxGetPr(prhs[6]);
    Eta_bra = mxGetPr(prhs[7]); 
    SEta_bra = mxGetPr(prhs[8]);
    maxIter = (int) mxGetScalar(prhs[9]);
    Tol = mxGetScalar(prhs[10]);
    verbosity = mxGetScalar(prhs[11]);
    mu_0 = mxGetPr(prhs[12]);
    
// Figure Out Dimensions
    k = (int) mxGetDimensions(prhs[0])[0]; // Mu
    n = (int) mxGetDimensions(prhs[0])[1];
    r = (int) mxGetDimensions(prhs[7])[0]; // Eta_Bra
    t = (int) mxGetDimensions(prhs[7])[1];
    s = (int) mxGetDimensions(prhs[8])[1]; // SEta_Bra
    if (s <= 0) mexErrMsgTxt("E_Distr statistics must not be empty");
    if(verbosity > 4){
        mexPrintf("k = %u, n = %u, t = %u, s = %u, r = %u, maxIter = %u, Tol = %2.2e\n",k,n,t,s,r,maxIter,Tol);    
        mexEvalString("drawnow;");
    } 
// Allocate Space
    M_cell = mxCreateCellMatrix(n,1);
    plhs[0] = mxCreateDoubleMatrix(k,n,mxREAL);
    mu_new = mxGetPr(plhs[0]);
    logM_null_mat = mxCreateDoubleMatrix(k,n,mxREAL);
    logM_null = mxGetPr(logM_null_mat);    
    if(verbosity > 4){
        mexPrintf("Stuff was allocated\n");
        mexEvalString("drawnow;");
    }
    
// Function Stuff
    
    // Calculate M_ij
    if(verbosity > 4){
        mexPrintf("Calculating M_cell\n");
        mexEvalString("drawnow;");
    }
    calc_M_cell(prhs, M_cell, R, Eta_bra, SEta_bra, r, n, k, t, s);
    
    // Mu Loop
    for(int iter = 0; iter < maxIter; iter++)
    {
        if(verbosity > 2)
        {
            mexPrintf("Iteration %u of %u\n",iter, maxIter);    
            mexEvalString("drawnow;");
        }
        
        // Update logM_null
        if(verbosity > 4){
            mexPrintf("Updating logM_null\n");
            mexEvalString("drawnow;");
        }  
        update_logM_null(mu, prhs, logM_null, R, SEta_bra, r, n, k, s);
        
        // Update Mes
        if(verbosity > 4){
            mexPrintf("Updating Mes\n");
            mexEvalString("drawnow;");
        } 
        update_mes(mu, prhs, M_cell, r, n, k, t, s, logM_null,mu_0);
//         mexEvalString("keyboard;");
        
        // Update Mu
        if(verbosity > 4){
            mexPrintf("Updating Mu\n");
            mexEvalString("drawnow;");
        }
        update_mu(mu, prhs, M_cell, n, k, mu_new, verbosity, logM_null,mu_0);
//         mexEvalString("keyboard;");
        
        // Calculate the Difference
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
    mxDestroyArray(M_cell);
    mxDestroyArray(logM_null_mat);
 return;   
}


////////////// Helper Function for calculating M_ij
void calc_M_cell(const mxArray *prhs[], mxArray *M_cell, 
        double* R, double* Eta_bra, double* SEta_bra, 
        int &r, int &n, int &k, int &t, int &s){
    // For each vertex i
    for(int i = 0; i < n; i++){
        // Allocate Space for M_cell
//         mexPrintf("Allocating Space for M_cell %u\n",i);
//         mexEvalString("drawnow;");
        int d_i_star = mxGetDimensions(mxGetCell(prhs[1],i))[1];
        mxSetCell(M_cell,i,mxCreateDoubleMatrix(1+k*k,d_i_star,mxREAL));
        double* M_ij = mxGetPr(mxGetCell(M_cell,i));
        
//         mexPrintf("Size of M_Cell, %u by %u ",
//                 mxGetDimensions(mxGetCell(M_cell,i))[0],
//                 mxGetDimensions(mxGetCell(M_cell,i))[1]);
        
//         mexPrintf("d_i_star, %u ",d_i_star);
//         mexEvalString("drawnow;");
        
        // Initialize M_ij to ones
        for(int dd = 0; dd < d_i_star; dd++){
            for(int k2 = 0; k2 < k; k2++){
                for(int k1 = 0; k1 < k; k1++){
                    // M_ij[k1,k2] = exp(0)
                    M_ij[k1+k2*k+dd*(1+k*k)] = 1;
                }
            }
        }
        
//         mexPrintf("...Done\n");
//         mexEvalString("drawnow;");
        // Calculate M_ij
        double* T_in = mxGetPr(mxGetCell(prhs[2],i));
        double* T_out = mxGetPr(mxGetCell(prhs[3],i));
        double* S_in = mxGetPr(mxGetCell(prhs[4],i));
        double* S_out = mxGetPr(mxGetCell(prhs[5],i));
        for(int d_cur = 0; d_cur < d_i_star; d_cur++){
            // Update Incoming S,T Terms
            if(S_in[d_cur*(1+s-1)+1] == 1){ 
                // Weighted Edge In
                for(int k2 = 0; k2 < k; k2++){
                    for(int k1 = 0; k1 < k; k1++){
                        double sum = 0;
                        int rr = (int) R[k2+k1*k]-1;
                        // Update Weight
                        for(int tt = 0; tt < t; tt++){
                            sum += T_in[d_cur*(1+t)+tt+1]*Eta_bra[rr+tt*r];
                        }
                        // Update Edge
                        for(int ss = 0; ss < s-1; ss++){
                            sum += S_in[d_cur*(1+s-1)+ss+1]*SEta_bra[rr+ss*r];
                        }
                        sum += SEta_bra[rr+(s-1)*r];
                        M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(sum);
                    }
                } 
            } else if(S_in[d_cur*(1+s-1)+1] == 0){
                // Non Edge
                for(int k2 = 0; k2 < k; k2++){
                    for(int k1 = 0; k1 < k; k1++){
                        int rr = (int) R[k2+k1*k]-1;
                        M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(SEta_bra[rr+(s-1)*r]);
                    }
                }
            } // Otherwise Missing Edge -> Do Nothing
            
            // Update Outgoing S,T Terms
            if(S_out[d_cur*(1+s-1)+1] == 1){ 
                // Weighted Edge out
                for(int k2 = 0; k2 < k; k2++){
                    for(int k1 = 0; k1 < k; k1++){
                        double sum = 0;
                        int rr = (int) R[k1+k2*k]-1;
                        // Update Weight
                        for(int tt = 0; tt < t; tt++){
                            sum += T_out[d_cur*(1+t)+tt+1]*Eta_bra[rr+tt*r];
                        }
                        // Update Edge
                        for(int ss = 0; ss < s-1; ss++){
                            sum += S_out[d_cur*(1+s-1)+ss+1]*SEta_bra[rr+ss*r];
                        }
                        sum += SEta_bra[rr+(s-1)*r];
                        M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(sum);
                    }
                } 
            } else if(S_out[d_cur*(1+s-1)+1] == 0){
                // Non Edge
                for(int k2 = 0; k2 < k; k2++){
                    for(int k1 = 0; k1 < k; k1++){
                        int rr = (int) R[k1+k2*k]-1;
                        M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(SEta_bra[rr+(s-1)*r]);
                    }
                }
            } // Otherwise Missing Edge -> Do Nothing        
        }
//         int d_cur = 0;
//         int t_in = 0; 
//         int t_out = 0;
//         while((t_in < d_in) && (t_out < d_out)){
// //             mexPrintf("Update: t_in %u of %u, t_out %u of %u, d_cur %u\n",
// //                     t_in,d_in,t_out,d_out,d_cur);
// //             mexEvalString("drawnow;");
//             if(T_in[t_in*(1+t)] > T_out[t_out*(1+t)]){
//                 // Update T_out
//                 update_M_ij(M_ij,T_out,R,Eta_bra,k,r,t,d_cur,t_out,1);
//                 t_out++; d_cur++;
//             }else if(T_in[t_in*(1+t)] < T_out[t_out*(1+t)]){
//                 // Update T_in
//                 update_M_ij(M_ij,T_in,R,Eta_bra,k,r,t,d_cur,t_in,0);
//                 t_in++; d_cur++;
//             }else{
//                 // Update T_in                
//                 update_M_ij(M_ij,T_in,R,Eta_bra,k,r,t,d_cur,t_in,0);
//                 t_in++;
//             }           
//         } 
//         while(t_in < d_in){
// //             mexPrintf("Update t_in %u of %u, d_cur %u\n",t_in,d_in,d_cur);
// //             mexEvalString("drawnow;");    
//             // Update T_in
//             update_M_ij(M_ij,T_in,R,Eta_bra,k,r,t,d_cur,t_in,0);
//             t_in++; d_cur++;
//         }
//         while(t_out < d_out){
//             // Update T_out
// //             mexPrintf("Update t_out %u of %u, d_cur %u\n",t_out,d_out,d_cur);
// //             mexEvalString("drawnow;");    
//             update_M_ij(M_ij,T_out,R,Eta_bra,k,r,t,d_cur,t_out,1);
//             t_out++; d_cur++;
//         }
//         
//         // Update S terms
//         d_cur = 0;
//         int s_in = 0;
//         int s_out = 0;
//         while((s_in < d_in) && (s_out < d_out)){
// //             mexPrintf("Update: s_in %u of %u, s_out %u of %u, d_cur %u\n",
// //                     s_in,d_in,s_out,d_out,d_cur);
// //             mexEvalString("drawnow;");
//             if(S_in[s_in*(1+s-1)] > S_out[s_out*(1+s-1)]){
//                 // Update S_out
//                 update_M_ijS(M_ij,S_out,R,SEta_bra,k,r,s,d_cur,s_out,1);
//                 s_out++; d_cur++;
//             }else if(S_in[s_in*(1+s-1)] < S_out[s_out*(1+s-1)]){
//                 // Update S_in
//                 update_M_ijS(M_ij,S_in,R,SEta_bra,k,r,s,d_cur,s_in,0);
//                 s_in++; d_cur++;
//             }else{
//                 // Update S_in                
//                 update_M_ijS(M_ij,S_in,R,SEta_bra,k,r,s,d_cur,s_in,0);
//                 s_in++;
//             }           
//         } 
//         while(s_in < d_in){
// //             mexPrintf("Update s_in %u of %u, d_cur %u\n",s_in,d_in,d_cur);
// //             mexEvalString("drawnow;");    
//             // Update S_in
//             update_M_ijS(M_ij,S_in,R,SEta_bra,k,r,s,d_cur,s_in,0);
//             s_in++; d_cur++;
//         }
//         while(s_out < d_out){
//             // Update S_out
// //             mexPrintf("Update s_out %u of %u, d_cur %u\n",s_out,d_out,d_cur);
// //             mexEvalString("drawnow;");    
//             update_M_ijS(M_ij,S_out,R,SEta_bra,k,r,s,d_cur,s_out,1);
//             s_out++; d_cur++;
//         }
// //         mexEvalString("keyboard;");
    }
    return;
}
// 
// ////////////// Update an element of M_ij for T
// void update_M_ij(double* M_ij, double* T, double* R, double* Eta, 
//         int &k, int &r, int &t, int &d_cur, int &t_cur, bool IsOut){
//     for(int k2 = 0; k2 < k; k2++){
//         for(int k1 = 0; k1 < k; k1++){
//             double sum = 0;
//             int rr;
//             if(IsOut){
//                 rr = (int) R[k1+k2*k]-1;
//             } else {
//                 rr = (int) R[k2+k1*k]-1;
//             }
//             for(int tt = 0; tt < t; tt++){
//                 sum += T[t_cur*(1+t)+tt+1]*Eta[rr+tt*r];
//             }             
//             M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(sum);
//         }
//     }
//     M_ij[k*k+d_cur*(1+k*k)] = T[t_cur*(1+t)];
//     return;
// }
// ////////////// Update an element of M_ij for S
// void update_M_ijS(double* M_ij, double* S, double* R, double* SEta, 
//         int &k, int &r, int &s, int &d_cur, int &s_cur, bool IsOut){
//     for(int k2 = 0; k2 < k; k2++){
//         for(int k1 = 0; k1 < k; k1++){
//             double sum = 0;
//             int rr;
//             if(IsOut){
//                 rr = (int) R[k1+k2*k]-1;
//             } else {
//                 rr = (int) R[k2+k1*k]-1;
//             }
//             for(int ss = 0; ss < s-1; ss++){
//                 sum += S[s_cur*(1+s-1)+ss+1]*SEta[rr+ss*r];
//             }
//             sum += SEta[rr+(s-1)*r];
//             M_ij[k1+k2*k+d_cur*(1+k*k)] *= exp(sum);
//         }
//     }
//     return;
// }

///////////// Update Mu
void update_mu(double* mu, const mxArray *prhs[], mxArray *M_cell,
        int &n, int &k, double* mu_new, double verbosity, 
        double *logM_null, double *mu_0){
    double* logmu;
    logmu = new double[k];
    double* sum_null; // Sum of logM_null parts 
    sum_null = new double[k];
    for(int k1 = 0; k1 < k; k1++){
        sum_null[k1] = 0;
    }
    for(int i = 0; i < n; i++){
        for(int k1 = 0; k1 < k; k1++){
            sum_null[k1] += logM_null[k1+k*i];
        }
    }   
    
    for(int i = 0; i < n; i++)
    {
        // Initialize mu_i 
        for(int k1 = 0; k1 < k; k1++){
            logmu[k1] = log(mu_0[k1+i*k])+sum_null[k1];
        }
        // Update mu_i with weighted messages and correcting for non-edges
        double* mes_i = mxGetPr(mxGetCell(prhs[1],i));
        int d_i_star = (int) mxGetDimensions(mxGetCell(prhs[1],i))[1];
        double* M_ij = mxGetPr(mxGetCell(M_cell,i));
        for(int dd = 0; dd < d_i_star; dd++)
        {
            int j = (int) mes_i[dd*(1+k)]-1;
            for(int k1 = 0; k1 < k; k1++){
                double sum = 0;
                for(int k2 = 0; k2 < k; k2++){
                    // sum_k2 M_jl[k1,k2]*mu_l[k2]
                    sum += M_ij[k1+k2*k+dd*(k*k+1)]*mes_i[dd*(k+1)+k2+1];
                }
                logmu[k1] += log(sum)-logM_null[k1+j*k]; 
            }
        }
        // Find maxlogmu
        double maxlogmu = logmu[0];
        for(int k1 = 1; k1 < k; k1++){
            if(logmu[k1] > maxlogmu){
                maxlogmu = logmu[k1];
            }
        }
        double maglogmu = 0;
        // logmu = exp(logmu - maxlogmu)
        for(int k1 = 0; k1 < k; k1++){ 
            logmu[k1] -= maxlogmu;
            logmu[k1] = exp(logmu[k1]);
            maglogmu += logmu[k1];
        }
        // Normalize and save mu_new
        for(int k1 = 0; k1 < k; k1++){
            mu_new[k1+i*k] = logmu[k1]/maglogmu;        
        }    
    }
    
    
    delete [] logmu;
    delete [] sum_null;
    return; 
}


////////////// Helper Function for calculating mes
void update_mes(double* mu, const mxArray *prhs[], mxArray *M_cell,
         int &r, int &n, int &k, int &t, int &s, double *logM_null,
         double *mu_0)
{
    
    double* logmes;
    logmes = new double[k];
    
    double* sum_null; // Sum of logM_null parts 
    sum_null = new double[k];
    for(int k1 = 0; k1 < k; k1++){
        sum_null[k1] = 0;
    }
    for(int i = 0; i < n; i++){
        for(int k1 = 0; k1 < k; k1++){
            sum_null[k1] += logM_null[k1+k*i];
        }
    }
    
    for(int i = 0; i < n; i++)
    {
        double* mes_i = mxGetPr(mxGetCell(prhs[1],i));
        int d_i_star = (int) mxGetDimensions(mxGetCell(prhs[1],i))[1];
        for(int dd = 0; dd < d_i_star; dd++)
        {
            // Calculate mes i -> j
            int j = (int) mes_i[dd*(1+k)]-1;
            if(i == j){ // Skip if i == j
                continue;
            }
            double* mes_j = mxGetPr(mxGetCell(prhs[1],j));
            int d_j_star = (int) mxGetDimensions(mxGetCell(prhs[1],j))[1];
            double* M_jl = mxGetPr(mxGetCell(M_cell,j));
//             mexPrintf("i %u, j %u, dd %u, djstar %u, M_jlsize %u \n",i,j,dd, d_j_star, mxGetDimensions(mxGetCell(M_cell,j))[1]);
//             mexEvalString("drawnow;");
            // Initialize mes i -> j 
            for(int k1 = 0; k1 < k; k1++){
                logmes[k1] = log(mu_0[k1+i*k])+sum_null[k1];
            }
            
            // Update mes_j with weighted messages and correcting for non-edges
            for(int d2 = 0; d2 < d_j_star; d2++){
                // Set Sum to Zero
                int ell = (int) mes_j[d2*(1+k)]-1;
//                 mexPrintf("ell %u, i %u, j%u, d2 %u of %u\n",ell,i,j,d2,d_j_star);
//                 mexEvalString("drawnow;");
                if(ell != i && ell != j){
                    for(int k1 = 0; k1 < k; k1++){
                        double sum = 0;
                        for(int k2 = 0; k2 < k; k2++){
                            // sum_k2 M_jl[k1,k2]*mu_l[k2]
                            sum += M_jl[k1+k2*k+d2*(k*k+1)]*mes_j[k2+1+d2*(k+1)]; 
//                             mexPrintf("M_jell k1 %u, k2 %u = %2.2e, mes_j %2.2e\n",k1,k2,M_jl[k1+k2*k+d2*(k*k+1)],mes_j[d2*(k+1)+k2+1]);
//                             mexEvalString("drawnow;");
                        }
//                         mexPrintf("k1 %u, sum %2.2e, log(sum) %2.2e\n",k1,sum,log(sum));
//                         mexEvalString("drawnow;");
                        logmes[k1] += log(sum)-logM_null[k1+k*ell]; 
                    }
                }
            }
//             mexPrintf("Find Max");
//             mexEvalString("drawnow;");
            // Find maxlogmes
            double maxlogmes = logmes[0];
            for(int k1 = 1; k1 < k; k1++){
                if(logmes[k1] > maxlogmes){
                    maxlogmes = logmes[k1];
                }
            }
//             mexPrintf(" ... %2.2e\nFind Mag",maxlogmes);
//             mexEvalString("drawnow;");
            double maglogmes = 0;
            // logmes = exp(logmes - maxlogmes)
            for(int k1 = 0; k1 < k; k1++){ 
                logmes[k1] -= maxlogmes;
                logmes[k1] = exp(logmes[k1]);
                maglogmes += logmes[k1];
            }
//             mexPrintf(" ... %2.2e\nSave Stuff \n",maglogmes);
//             mexEvalString("drawnow;");
            // Normalize and save mes i->j
            for(int k1 = 0; k1 < k; k1++){
                mes_i[k1+1+dd*(k+1)] = logmes[k1]/maglogmes;        
            }
//             mexEvalString("keyboard;");
        }
    }   
    delete [] logmes;
    delete [] sum_null;
    return;
}

//////////////// Update logM_null the messages for Non-Edges
void update_logM_null(double *mu, const mxArray *prhs[], double *logM_null, 
        double *R, double *SEta_bra, int &r, int &n, int &k, int &s)
{
    // Set logM_null to zero
    for(int i = 0; i < n; i++){
        for(int k1 = 0; k1 < k; k1++){
            logM_null[k1+k*i] = 0;
        }
    }
if(s > 0){
    for(int i = 0; i < n; i++){
        for(int k1 = 0; k1 < k; k1++){
            double sum = 0;
            for(int k2 = 0; k2 < k; k2++){
                // Non DC Case
                sum += exp(SEta_bra[((int)R[k2+k1*k]-1)+(s-1)*r]+
                       SEta_bra[((int)R[k1+k2*k]-1)+(s-1)*r])*mu[k2+k*i]; 
            }
            logM_null[k1+k*i] = log(sum);
        }
    }
}          
    return;
}

// The End