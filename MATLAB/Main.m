%  Main evaluation code for "Link Prediction Based on Graph Neural Networks".
%
%  *author: Muhan Zhang, Washington University in St. Louis

delete(gcp('nocreate'))
addpath(genpath('utils'));
datapath = 'data/';

% change your general exp setting here
setting = 1;
switch setting
case 1  % traditional link prediction benchmarks
    numOfExperiment = 1;        
    %workers = 5;  % number of workers running parallelly
    workers = 0;  % change workers to 0 to disable parallel loop of multiple exps
    ratioTrain = 0.9; % train split ratio
    connected = false; % whether to sample test links while ensuring the remaining net is connected
    %dataname = strvcat('USAir','NS','PB','Yeast','Celegans','Power','Router','Ecoli');
    %dataname = strvcat('USAir','NS','Yeast','Celegans','Power','Router'); % 
    %dataname = strvcat('PB', 'Ecoli');  % set workers 5,  h=1 for SEAL due to memory issues
    %dataname = strvcat('PB', 'Ecoli'); % set workers 2, h=1 for WL alone due to memory issues
    dataname = strvcat('USAir');
    %method = [1, 2, 3, 4, 5, 6, 7, 8, 9];  % 1: SEAL,  2: Heuristic methods 3: Traditional latent feature methods,  4: WLNM 5: WL graph kernel, 6: Embedding methods
    %method =[1, 2, 3, 4, 5, 6];
    method =[1];
    h = 'auto';  % the maximum hop to extract enclosing subgraphs, h = 'auto' means to automatically select h from 1, 2
    %h = 1;
    include_embedding = 0;  % whether to include node embeddings in node information matrix of SEAL, needs node2vec software
    include_attribute = 0;
    portion = 1;  % portion of observed links selected as training data
case 2  % network embedding benchmark datasets without node attributes
    numOfExperiment = 5;        
    workers = 0;  % disable parallel computing for large networks
    ratioTrain = 0.5; 
    dataname = strvcat('facebook', 'arxiv');  
    connected = true;
    method = [1];
    h = 1;
    include_embedding = 1;
    include_attribute = 0;
    portion = 1;
case 3  % network embedding benchmark datasets with node attributes
    numOfExperiment = 5;        
    workers = 0;  % disable parallel computing for large networks
    ratioTrain = 0.5; 
    dataname = strvcat('PPI_subgraph', 'Wikipedia', 'BlogCatalog');  % networks with node attributes
    connected = true;
    method = [1];
    h = 1;
    include_embedding = 1;
    include_attribute = 1;
    portion = 10000;  % randomly select 10000 observed links as positive training (since selecting all is out of memory)
case 4  % cora, citeseer, pubmed
    numOfExperiment = 5;        
    workers = 5;  % disable parallel computing for large networks
    ratioTrain = 0.9; 
    dataname = strvcat('cora', 'citeseer', 'pubmed');  % networks with node attributes
    connected = false;
    method = [1];
    h = 1;
    include_embedding = 1;
    include_attribute = 1;
    portion = 1; 
end

tic;
num_in_each_method = [1, 9, 2, 1, 1, 3];  % how many algorithms in each type of method
num_of_methods = sum(num_in_each_method(method));  % the total number of algorithms

auc_for_dataset = [];
precision_for_dataset = [];
for ith_data = 1:size(dataname, 1)                          
    tempcont = ['processing the ', int2str(ith_data), 'th dataset...', dataname(ith_data,:)];
    disp(tempcont);
    thisdatapath = strcat(datapath,dataname(ith_data,:),'.mat');  % load the net
    load(thisdatapath);                                 
    aucOfallPredictor = zeros(numOfExperiment, num_of_methods); 
    precisionOfallPredictor = zeros(numOfExperiment, num_of_methods); 
    PredictorsName = [];
    
    % parallelize the repeated experiments
    %for ith_experiment = 1:numOfExperiment
    parfor (ith_experiment = 1:numOfExperiment, workers)
        ith_experiment
        if mod(ith_experiment, 10) == 0
                tempcont = strcat(int2str(ith_experiment),'%... ');
                disp(tempcont);
        end

        rng(ith_experiment);  % generate fixed network splits for different methods and exp numbers

        % divide net into train/test
        [train, test] = DivideNet(net, ratioTrain, connected); % train test are now symmetric adjacency matrices without self loops 
        train = sparse(train); test = sparse(test);  % convert to sparse matrices

        % sample negative links 
        htrain = triu(train, 1);  % half train adjacency matrix
        htest = triu(test, 1);
        [train_pos, train_neg, test_pos, test_neg] = sample_neg(htrain, htest, 1, portion, false);
        test = {};
        test.pos = test_pos; test.neg = test_neg;  % evaluate performance on sampled test links
        train_mix = {};
        train_mix.train = train;  % the observed network to extract enclosing subgraphs
        train_mix.data_name = strip(dataname(ith_data, :));  % store the name of the current data set
        train_mix.pos = train_pos; train_mix.neg = train_neg;  % the train pos and neg links (used by learning based methods, not by heuristic methods)

        ithAUCvector = []; ithprecisionvector = []; Predictors = []; % for recording results

        % Run link prediction methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SEAL
        if ismember(1, method)
        disp('SEAL...');
        [auc, precision] = SEAL(train_mix, test, h, include_embedding, include_attribute, ith_experiment);
            Predictors = [Predictors 'SEAL	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Heuristic methods
        if ismember(2, method)

        all_sims = {}

        disp('CN...');
        [auc, precision, sim] = CN(train, test);                  % Common Neighbor
            Predictors = [Predictors 'CN	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
        
        disp('Jaccard...');
        [auc, precision, sim] = Jaccard(train, test);             % Jaccard Index
            Predictors = [Predictors 'Jaccard	'];  ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim

        disp('PA...');
        [auc, precision, sim] = PA(train, test);                  % Preferential Attachment
            Predictors = [Predictors 'PA	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
        
        disp('AA...');
        [auc, precision, sim] = AA(train, test);                  % Adar-Adamic Index
            Predictors = [Predictors 'AA	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
        
        disp('RA...');
        [auc, precision, sim] = RA(train, test);                  % Resourse Allocation
            Predictors = [Predictors 'RA	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
       
        disp('Katz 0.001...');
        [auc, precision, sim] = Katz(train, test, 0.001);         % Katz Index, beta=0.001
            Predictors = [Predictors '~.001	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim

        disp('RWR 0.85...');
        [auc, precision, sim] = RWR(train, test, 0.85);           % Random walk with restart (rooted PageRank), d=0.85
            Predictors = [Predictors 'RWR.85	'];   ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
        
        disp('SimRank 0.8...');
        [auc, precision, sim] = SimRank(train, test, 0.8);        % SimRank
            Predictors = [Predictors 'SimR	'];      ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            all_sims{length(all_sims)+1} = sim
        
        disp('Ensemble heuristics...');
        [auc, precision] = ensemble_heuristics(train_mix, test, all_sims);      % Logistic regression on all heuristics
            Predictors = [Predictors 'Ensemble  '];      ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Traditional latent feature methods
        if ismember(3, method)
        % matrix factorization
        disp('MF...');
        [auc, precision] = MF(train, test, 5, ith_experiment);                 % matrix factorization
            Predictors = [Predictors 'MF	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
            
        % stochastic block model
        disp('SBM...');
        [auc, precision] = SBM(train, test, 12);                 % stochastic block models
            Predictors = [Predictors 'SBM	'];       ithAUCvector = [ithAUCvector auc];  
            ithprecisionvector = [ithprecisionvector precision];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Weisfeiler-Lehman Neural Machine (WLNM)
        if ismember(4, method)
        disp('WLNM...');
        [auc, precision] = WLNM(train_mix, test, 10, ith_experiment);                  % WLNM
            Predictors = [Predictors 'WLNM	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Weisfeiler-Lehman graph kernel for Link Prediction
        if ismember(5, method)
        disp('WLK...');
        [auc, precision] = WLK(train_mix, test, h, ith_experiment);
            Predictors = [Predictors 'WLK	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Embedding + LR for Link Prediction
        if ismember(6, method)
        disp('Embedding...');
        [auc, precision] = embedding_lp(train_mix, test, 'node2vec', ith_experiment);
            Predictors = [Predictors 'node2vec	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
        
        disp('Embedding...');
        [auc, precision] = embedding_lp(train_mix, test, 'LINE', ith_experiment);
            Predictors = [Predictors 'LINE	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];

        disp('Embedding...');
        [auc, precision] = embedding_lp(train_mix, test, 'SPC', ith_experiment);
            Predictors = [Predictors 'SPC	'];      ithAUCvector = [ithAUCvector auc];
            ithprecisionvector = [ithprecisionvector precision];
        end

        aucOfallPredictor(ith_experiment, :) = ithAUCvector; PredictorsName = Predictors;
        precisionOfallPredictor(ith_experiment, :) = ithprecisionvector;
    end
    if exist('poolobj')
        delete(poolobj)
    end

    %% write the results for this dataset
    avg_auc = mean(aucOfallPredictor,1)
    avg_precision = mean(precisionOfallPredictor,1)
    auc_for_dataset = [auc_for_dataset, avg_auc];
    precision_for_dataset = [precision_for_dataset, avg_precision];
    std_auc = std(aucOfallPredictor, 0, 1)
    std_precision = std(precisionOfallPredictor, 0, 1)
    respath = strcat(datapath,'result/',dataname(ith_data,:),'_res.txt');         
    dlmwrite(respath,{PredictorsName}, '');
    dlmwrite(respath,[avg_auc; std_auc; avg_precision; std_precision], '-append','delimiter', '	','precision', 4);
    
end 
toc;
auc_for_dataset'
precision_for_dataset'


