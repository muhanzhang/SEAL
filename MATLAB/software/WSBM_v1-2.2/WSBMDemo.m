% WSBMDEMO Demo Script for the WSBM
%   Press Ctrl + Enter in each (%%) cell to execute the code
% See also WSBM
% Version 1.0 | January 2014 | Christopher Aicher
% Version 1.1 | May 2014 | Added NFL2009 Dataset Example

%% Add Analysis Tool Dir
addpath('analysis tools'); 

%% Install MEX Files (Optional)
InstallMEXFiles

%% Simple Example----------------------------------------------------------
%% Make Adj_Matrix
% A is a 20x20 Adjacency Matrix for two groups of 10 vertices each,
% Edges within groups have weight 1 and edges between groups have weight 0
A = [ones(10,10),zeros(10,10);
     zeros(10,10),ones(10,10)];
% Visualization of the Adjacency Matrix A
%plotWSBM(A);
%title('Synthetic Data Adjacency Matrix');
%% Infer Truth
% Run the inference on A for 2 groups
Labels = wsbm(A,2);
% Display the inferred labels
disp('Labels:');
disp(Labels');
%% Shuffled Adj_Matrix
% Shuffle the numbering of the vertices of A 
%  (i.e. permute the rows and columns)
Ashuffle = shuffle(A);
%plotWSBM(Ashuffle);
%title('Shuffled Synthetic Data Adjacency Matrix');
%% Infer Truth
% Run the inference on Ashuffle for 2 groups
[Labels, model] = wsbm(Ashuffle,2);
% Display the inferred labels
disp('Labels:');
disp(Labels');

%% 2 Group Weighted Network Example----------------------------------------
%% Generate Random Synthetic Data
[Edge_List,True_Model] = generateEdges();
%plotWSBM(True_Model);
%title('Synthetic Data Adjacency Matrix');
%% Fit WSBM
num_groups = 2; 
[Labels, Model] =  wsbm(Edge_List,num_groups);
plotWSBM(Model);
title('Inferred Permuted Adjacency Matrix');

%% 4 Group Mixed Network Example-------------------------------------------
%% Generate Data
% Advanced User Stuff -> It creates an interesting dataset
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1];
theta_w = [100,1;100,1; 0,1; 0,1];
theta_e = [0.1; 0.9; 0.1; 0.9];
group_sizes = [25;25;25;25];
[~,True_Model] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
plotWSBM(True_Model);
title('Synthetic Data');
%% Fit Edge Model
[~,Edge_Model] = wsbm(True_Model.Data.Raw_Data,2,'W_Distr','None','E_Distr','Bernoulli'); 
subplot(1,2,1);
plotMu(Edge_Model);
subplot(1,2,2);
plotWSBM(Edge_Model);
title('Edge Model Permuted Adjacency Matrix')
%% Fit Weight Model
[~,Weight_Model] = wsbm(True_Model.Data.Raw_Data,2,'W_Distr','Normal','E_Distr','None'); 
subplot(1,2,1);
plotMu(Weight_Model);
subplot(1,2,2);
plotWSBM(Weight_Model);
title('Weight Model Permuted Adjacency Matrix')
%% Fit Mixed Model
[~,Mixed_Model] = wsbm(True_Model.Data.Raw_Data,4);
subplot(1,2,1);
plotMu(Mixed_Model);
subplot(1,2,2);
plotWSBM(Mixed_Model);
title('Mixed Model Permuted Adjacency Matrix')

%% NFL2009 Example --------------------------------------------------------
%% Load Edge List Data
data = importdata(['NFL2009_network',filesep,'NFL2009_EdgeList.txt']);
E = zeros(2*size(data.data,1),3);
E(1:size(data.data,1),1:2) = data.data(:,1:2);
E(size(data.data,1)+1:end,1:2) = [data.data(:,2),data.data(:,1)];
% Weights = Score Difference
E(1:size(data.data,1),3) = data.data(:,3)-data.data(:,4);
E(size(data.data,1)+1:end,3) = data.data(:,4)-data.data(:,3); 

%% Load Vertex Metadata
data = importdata(['NFL2009_network',filesep,'NFL2009_VertexMetadata.txt']);
mu_conf = zeros(2,size(data.data,1));
for ii = 1:2,
    mu_conf(ii,data.data(:,1) == ii) = 1;
end
mu_division = zeros(8,size(data.data,1));
for ii = 1:8,
    mu_division(ii,data.data(:,2) == ii) = 1;
end
%% Plot Raw Edge List
subplot(1,2,1)
plotWSBM(E);
title('Plot of Avg. Weight');
subplot(1,2,2)
plotWSBM(E,'edges'); 
title('Plot of Number of Edges');

%% Plot Edge List Sorted by Conference
subplot(1,2,1)
plotWSBM(E,mu_conf);
title('Plot of Avg. Weight');
subplot(1,2,2)
plotWSBM(E,mu_conf,'edges'); 
title('Plot of Number of Edges');

%% Plot Edge List Sorted by Division
subplot(1,2,1)
plotWSBM(E,mu_division);
title('Plot of Avg. Weight');
subplot(1,2,2)
plotWSBM(E,mu_division,'edges'); 
title('Plot of Number of Edges');

%% Infer Models
[~,Weight_Model] = wsbm(E,4,'W_Distr','Normal','E_Distr','None','numTrials',100);
[~,Edge_Model] = wsbm(E,4,'W_Distr','None','E_Distr','Poisson','numTrials',100);

%% Plot Inferred Model
subplot(1,2,1)
plotWSBM(Weight_Model);
title('Plot of Avg. Weight Sorted By Weight WSBM');
subplot(1,2,2)
plotWSBM(Edge_Model,'edges'); 
title('Plot of Edge Counts Sorted By Edge SBM');


