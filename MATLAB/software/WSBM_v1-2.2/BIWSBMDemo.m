% BIWSBMDEMO Demo Script for the BIWSBM
%   Press Ctrl + Enter in each (%%) cell to execute the code
% See also BIWSBM, WSBM, WSBMDEMO
% Version 1.0 | April 2014 | Christopher Aicher

%% Add Analysis Tool Dir
addpath('analysis tools'); 

%% Install MEX Files (Optional)
InstallMEXFiles

%% Simple Example----------------------------------------------------------
%% Make Adj_Matrix
% A is a 40x40 Adjacency Matrix for two groups of 10 vertices for each type
A = [ones(10,10),zeros(10,10);
     zeros(10,10),ones(10,10)];
A = [zeros(20),A; A, zeros(20)];
% Visualization of the Adjacency Matrix A
plotWSBM(A);
title('Synthetic Data Adjacency Matrix');
%% Infer Truth Using BIWSBM
% Run the inference on A for (2,2) groups
K_0 = 2; 
K_1 = 2;
types = zeros(40,1);
types(21:40) = 1;
[~,Model] = biwsbm(A,K_0,K_1,types,'W_Distr','Bernoulli','E_Distr','None');
% Display the inferred labels
subplot(1,2,1);
plotWSBM(Model);
subplot(1,2,2);
plotMu(Model);

%% Advanced Example 1------------------------------------------------------
%% Generate Bipartite Graph
R = [1,1,1,2,3; 
     1,1,1,3,2;
     1,1,1,4,4;
     2,3,4,1,1;
     3,2,4,1,1];
theta_w = [0; 2; 10; 1];
theta_e = [0; 1; 1; 1];
group_sizes = [20;20;20;20;20];
[~,True_Model] = generateEdges('Poisson','Bernoulli',R,theta_w,theta_e,group_sizes);
plotWSBM(True_Model);
title('Synthetic Data');

%% Infer Truth Using BIWSBM
K_0 = 3; % Number of Groups of Type 1 
K_1 = 2; % Number of Groups of Type 2
types = [zeros(60,1);ones(40,1)]; % Vertex Type Indicator
[~,Bipartite_Model] = biwsbm(True_Model.Data.Raw_Data,K_0,K_1,types,...
    'W_Distr','Poisson','E_Distr','None'); 
subplot(1,2,1);
plotMu(Bipartite_Model);
subplot(1,2,2);
plotWSBM(Bipartite_Model);
title('Bipartite Model Permuted Adjacency Matrix');

%% Compare with Using a Naive WSBM
K = 5;
[~,Naive_Model] = wsbm(True_Model.Data.Raw_Data,K,'W_Distr','Normal',...
    'E_Distr','None');
subplot(1,2,1);
plotMu(Naive_Model);
subplot(1,2,2);
plotWSBM(Naive_Model);
title('Naive Model Permuted Adjacency Matrix')


%% Advanced Example 2------------------------------------------------------
%% Generate Bipartite Graph
R = [1,1,1,2,3; 
     1,1,1,3,2;
     1,1,1,4,4;
     2,3,4,1,1;
     3,2,4,1,1];
theta_w = [];
theta_e = [0; 4; 2; 1]/10;
group_sizes = [20;20;20;20;20];
degree_Para = 5*rand(100,2); 
[~,True_Model] = generateEdges('None','DC',R,theta_w,theta_e,group_sizes,degree_Para);
plotWSBM(True_Model,'edges');
title('Synthetic Data');

%% Infer Truth Using BIWSBM
K_0 = 3; 
K_1 = 2;
types = [zeros(60,1);ones(40,1)];
[~,Bipartite_Model] = biwsbm(True_Model.Data.Raw_Data,K_0,K_1,types,...
    'W_Distr','None','E_Distr','DC'); 
subplot(1,2,1);
plotMu(Bipartite_Model);
subplot(1,2,2);
plotWSBM(Bipartite_Model,'edges');
title('Bipartite Model Permuted Adjacency Matrix')

%% Infer Truth Using Naive Non-DC Model
K_0 = 3; 
K_1 = 2;
types = [zeros(60,1);ones(40,1)];
[~,Bipartite_Model] = biwsbm(True_Model.Data.Raw_Data,K_0,K_1,types,...
    'W_Distr','None','E_Distr','Poisson'); 
subplot(1,2,1);
plotMu(Bipartite_Model);
subplot(1,2,2);
plotWSBM(Bipartite_Model,'edges');
title('Bipartite Model Permuted Adjacency Matrix')



