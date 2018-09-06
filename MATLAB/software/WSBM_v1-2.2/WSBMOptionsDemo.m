% WSBMOPTIONDEMO Demo Script for the WSBM
%   This demo highlights how to change some of the default inference 
%   options including the prior distribution of edge distributions, 
%   vertex-labels and tolerance/iterations values.
%
%   Press Ctrl + Enter in each (%%) cell to execute the code
% See also WSBM, WSBMDEMO
% Version 1.0 | May 2014 | Christopher Aicher
%

%% Add Analysis Tool Dir
addpath('analysis tools'); 

%% Install MEX Files (Optional)
InstallMEXFiles

%% Generate Data 
[Edge_List] = generateEdges();
Edge_List(:,3) = Edge_List(:,3)+20;
plotWSBM(Edge_List);

%% Fit Default WSBM
[~,Model] = wsbm(Edge_List);
% Number of Groups = 4
% Weight Distribution is Normal/Gaussian
% Edge Distribution is Poisson
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Vary Number of Groups
[~,Model] = wsbm(Edge_List,2)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian
% Edge Distribution is Poisson
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Vary Weight Distribution
[~,Model] = wsbm(Edge_List,2,'W_Distr','Exp')
% Number of Groups = 2
% Weight Distribution is Exponential
% Edge Distribution is Poisson
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Vary Weight Distribution Prior
% Very Strong Prior Distribution (Don't do this) 
normal_distr = setup_distr('Normal',[0,100,1000],1);
% See `help setup_distr' for more details
[~,Model] = wsbm(Edge_List,2,'W_Distr',normal_distr)
% Number of Groups = 2
% Weight Distribution is Normal with very strong prior
% Edge Distribution is Poisson
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Ignore the Weight Distribution
[~,Model] = wsbm(Edge_List,2,'W_Distr','None')
% Number of Groups = 2
% Weight Distribution is Ignored (None).
% Edge Distribution is Poisson
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Vary Edge Distribution
[~,Model] = wsbm(Edge_List,2,'E_Distr','Bernoulli')
% Number of Groups = 2
% Weight Distribution is Normal
% Edge Distribution is Bernoulli
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model)

%% Vary Edge Distribution Prior
% Biased Prior Distribution 
bernoulli_distr = setup_distr('Binomial',[0,20],1);
% See `help setup_distr' for more details
[~,Model] = wsbm(Edge_List,2,'E_Distr',bernoulli_distr,'W_Distr','None')
% Number of Groups = 2
% Weight Distribution is None
% Edge Distribution is Bernoulli with a biased prior
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Ignore the Edge Distribution
[~,Model] = wsbm(Edge_List,2,'E_Distr','None')
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Ignored (None).
% Number of Random Initializations is 50
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Vary Prior over Vertex Labels
mu_0 = zeros(2,max([Edge_List(:,1);Edge_List(:,2)]));
mu_0(1,:) = 1/4;
mu_0(2,:) = 3/4;
[~,Model] = wsbm(Edge_List,2,'mu_0',mu_0)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Ignored (None).
% Number of Random Initializations is 50
% Prior favors group 2 three times more than group 1
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Run Inference Silently
[~,Model] = wsbm(Edge_List,2,'verbosity',0)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson
% Number of Random Initializations is 50
% No Messages

%% Display Extra Messages
[~,Model] = wsbm(Edge_List,2,'verbosity',2)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson
% Number of Random Initializations is 50
% Displays Additional Convergence Messages

%% Run Inference in Parallel
[~,Model] = wsbm(Edge_List,2,'parallel',1)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson
% Number of Random Initializations is 50
% Uses MATLABPOOL (if you have it the toolbox) to run wsbm.m in parallel

%% Run Additional Random Initialized Trials
[~,Model] = wsbm(Edge_List,2,'numTrials',200)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson
% Number of Random Initializations is 50
% Displays Additional Convergence Messages

%% Change Convergence Tolerance
[~,Model] = wsbm(Edge_List,2,'muTol',0.1,'mainTol',0.1)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson.
% Number of Random Initializations is 50
% Changed Convergence Tolerance from default (1e-3) to (1e-1)
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Change Convergence Max Iterations
[~,Model] = wsbm(Edge_List,2,'muMaxIter',5,'mainMaxIter',5)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson.
% Number of Random Initializations is 50
% Changed Convergence Max Iterations from default 80 and 50 to 5 and 5
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);

%% Run From Seed
seed = ones(2,max([Edge_List(:,1);Edge_List(:,2)]))/2;%(This is a bad seed)
[~,Model] = wsbm(Edge_List,2,'seed',seed)
% Number of Groups = 2
% Weight Distribution is Normal/Gaussian.
% Edge Distribution is Poisson
% Number of Random Initializations is 50
% Runs Algorithm Once From Seed
subplot(1,2,1);
plotMu(Model);
subplot(1,2,2);
plotWSBM(Model);


