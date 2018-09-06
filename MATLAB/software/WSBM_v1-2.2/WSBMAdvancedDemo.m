% WSBMDEMO Demo Script for the WSBM
%   Press Ctrl + Enter in each (%%) cell to execute the code
% See also WSBM
% Version 1.0 | January 2014 | Christopher Aicher

%% Add Analysis Tool Dir
addpath('analysis tools'); 

%% Install MEX Files (Optional)
InstallMEXFiles

%% Fit Multiple Models at Once    -----------------------------------------
%% Generate Data + Setup Models + Setup Score
% Generate Data
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1];
theta_w = [100,1;100,1; 0,1; 0,1];
theta_e = [0.1; 0.9; 0.1; 0.9];
group_sizes = [25;25;25;25];
[E,True_Model] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
%plotWSBM(True_Model);
title('Synthetic Data');

% Create a List of Model Inputs to test
Model1 = {2,'W_Distr','None','E_Distr','Bernoulli'}; % A pure SBM
Model2 = {2,'W_Distr','Normal','E_Distr','None'};    % A pure WSBM
Model3 = {4,'W_Distr','Normal','E_Distr','Bernoulli','numTrials',100}; % A mixed Model
ModelInputs = {Model1; Model2; Model3};

% Create a List of model score functions
scoreNMI = @(Model) nmi(True_Model,Model);
scoreVI = @(Model) varInfo(True_Model,Model);
scorefuncs = {scoreNMI;scoreVI};

%% Fit + Score Multiple Models
[Best_Model,Scores,Models] = wsbmLooper(E,ModelInputs,scorefuncs);

%% Plot Models
% Note that the mixed model (row 3) recovers the correct structure)
% subplot(3,2,1);
% plotMu(Models{1});
% subplot(3,2,2);
% plotWSBM(Models{1});
% subplot(3,2,3);
% plotMu(Models{2});
% subplot(3,2,4);
% plotWSBM(Models{2});
% subplot(3,2,5);
% plotMu(Models{3});
% subplot(3,2,6);
% plotWSBM(Models{3});

%% Learn K Example---------------------------------------------------------
%% Generate Data 
% Generate Data
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1];
theta_w = [100,1;100,1; 0,1; 0,1];
theta_e = [0.1; 0.9; 0.1; 0.9];
group_sizes = [25;25;25;25];
[E,True_Model] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
%plotWSBM(True_Model);
title('Synthetic Data');
%% Model Selection for the number of groups K
% Create a List of Model Inputs to test
R_Structs = 2:1:8;      %   Test K = 2,3,...,8
W_Distr = 'Normal';     %   Use Normal Weight Distribution:
E_Distr = 'Bernoulli';  %   Use Bernoulli Edge Distribution:
ModelInputs = cell(numel(R_Structs),1); 
for rr = 1:numel(R_Structs),
    ModelInputs{rr} = {R_Structs(rr),'W_Distr',W_Distr,'E_Distr',E_Distr};
end

% Fit the list of models and select the Best according to LogEvidence
[Best_Model,Scores,Models] = wsbmLooper(E,ModelInputs);

%% Alternatively, Fit and Score Models according to various criteria
% Create a List of Model Inputs to test
R_Structs = 2:1:8;      %   Test K = 2,3,...,8
W_Distr = 'Normal';     %   Use Normal Weight Distribution:
E_Distr = 'Bernoulli';  %   Use Bernoulli Edge Distribution:
ModelInputs = cell(numel(R_Structs),1); 
for rr = 1:numel(R_Structs),
    ModelInputs{rr} = {R_Structs(rr),'W_Distr',W_Distr,'E_Distr',E_Distr};
end

% Create a List of model score functions
LogEvidence = @(Model) Model.Para.LogEvidence;      % LogEvidence
scoreNMI = @(Model) nmi(True_Model,Model);          % NMI
scoreVI = @(Model) -varInfo(True_Model,Model);      % VI
Scorefuncs = {LogEvidence;scoreNMI;scoreVI};

[~,Scores,Models] = wsbmLooper(E,ModelInputs,Scorefuncs);

%% Plot Resulting Model Scores
plot(R_Structs',Scores(:,1));
xlabel('K');
ylabel('LogEvidence');
subplot(1,3,2);
plot(R_Structs',Scores(:,2));
xlabel('K');
ylabel('NMI');
subplot(1,3,3);
plot(R_Structs',Scores(:,3));
xlabel('K');
ylabel('VI');

%% Cross Validation Example -----------------------------------------------
%% Generate Data
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1];
theta_w = [8,1;10,1; 2,1; 0,1];
theta_e = [0.25; 0.9; 0.1; 0.75];
group_sizes = [30;30;30;30];
[E,True_Model] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
plotWSBM(True_Model);
title('Synthetic Data');
% Split into Training and Test Datasets 
[E_Train,E_Test] = crossValidEdg(E,0.1);

%% Fit + Score Various Models
% Create a List of Model Inputs to test
Model1 = {4,'W_Distr','None','E_Distr','Bernoulli','parallel',1,'numTrials',100};
Model2 = {4,'W_Distr','Normal','E_Distr','None','parallel',1,'numTrials',100};
Model3 = {4,'W_Distr','Normal','E_Distr','Bernoulli','parallel',1,'numTrials',100};
ModelInputs = {Model1; Model2; Model3};

% Create a List of model score functions
Error = @(x,y) (x-y).^2; % MSE
score1 = @(Model) -predictW_Error(Model,E_Test,Error);
score2 = @(Model) -predictE_Error(Model,E_Test,Error);
scorefuncs = {score1,score2};

% Fit
[Best_Model,Scores,Models] = wsbmLooper(E_Train,ModelInputs,scorefuncs);
disp(Scores) % Less Negative == Better
% Note that by combining Edge + Weight Information, it is possible to
% outperform both the Pure SBM and the Pure WSBM in this (specific) case.

%% Plot Models
subplot(3,2,1);
plotMu(Models{1});
subplot(3,2,2);
plotWSBM(Models{1});
subplot(3,2,3);
plotMu(Models{2});
subplot(3,2,4);
plotWSBM(Models{2});
subplot(3,2,5);
plotMu(Models{3});
subplot(3,2,6);
plotWSBM(Models{3});
