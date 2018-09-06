% DegreeCorrectedDEMO Demo Script for the WSBM
%   Press Ctrl + Enter in each (%%) cell to execute the code
% See also WSBM, WSBMDEMO
% Version 1.0 | April 2014 | Christopher Aicher

%% Add Analysis Tool Dir
addpath('analysis tools'); 

%% Install MEX Files (Optional)
InstallMEXFiles

%% Simple Example ---------------------------------------------------------
%% Generate Data
R = [1,2; 3,4];
theta_w = [];
theta_e = [1; 4; 2; 3];
group_sizes = [50;50];
degree_Para = 2*rand(sum(group_sizes),2);
[~,True_Model] = generateEdges('None','DC',R,theta_w,theta_e,...
    group_sizes,degree_Para);
plotWSBM(True_Model,'edges');
title('Synthetic Data');
%% Fit Degree Corrected Model
[~,DC_Model] = wsbm(True_Model.Data.Raw_Data,2,'W_Distr','None','E_Distr','DC');
subplot(1,2,1);
plotMu(DC_Model);
subplot(1,2,2);
plotWSBM(DC_Model,'edge');
title('DC Model Permuted Adjacency Matrix')
%% Fit Naive Poisson Model
[~,Naive_Model] = wsbm(True_Model.Data.Raw_Data,2,'W_Distr','None','E_Distr','Poisson');
subplot(1,2,1);
plotMu(Naive_Model);
subplot(1,2,2);
plotWSBM(Naive_Model,'edge');
title('Naive Poisson Model Permuted Adjacency Matrix')
% Note that the Naive Poisson Model sorts by degree instead of connectiv

%% Advanced Example -------------------------------------------------------
%% Generate Data
% Advanced User Stuff -> It creates an interesting dataset
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1];
theta_w = [100,1;100,1; 0,1; 0,1];
theta_e = [0.1; 0.9; 0.1; 0.9];
group_sizes = [25;25;25;25];
degree_Para = rand(sum(group_sizes),2)*5;
[~,True_Model] = generateEdges('Normal','DC',R,theta_w,theta_e,...
    group_sizes,degree_Para);
subplot(1,2,1);
plotWSBM(True_Model);
title('Synthetic Data: Average Weight');
subplot(1,2,2);
plotWSBM(True_Model,'edge');
title('Synthetic Data: Number of Edges');

%% Fit Degree Corrected WSBM
[~,DC_Model] = wsbm(True_Model.Data.Raw_Data,4,'W_Distr','Normal','E_Distr','DC');
subplot(1,2,1);
plotMu(DC_Model);
subplot(1,2,2);
plotWSBM(DC_Model,'edge');
title('DC Model Permuted Adjacency Matrix')
%% Fit Naive WSBM
[~,Naive_Model] = wsbm(True_Model.Data.Raw_Data,4,'W_Distr','Normal','E_Distr','Poisson');
subplot(1,3,1);
plotMu(Naive_Model);
subplot(1,3,2);
plotWSBM(Naive_Model,'edge');
title('Naive Poisson with Weights Permuted Adjacency Matrix')
subplot(1,3,3);
plotWSBM(Naive_Model);
%% Fit Naive SBM
[~,Naive_Model] = wsbm(True_Model.Data.Raw_Data,4,'W_Distr','None','E_Distr','Poisson');
subplot(1,3,1);
plotMu(Naive_Model);
subplot(1,3,2);
plotWSBM(Naive_Model,'edge');
title('Naive Poisson ignoring Weights Permuted Adjacency Matrix')
subplot(1,3,3);
plotWSBM(Naive_Model);


