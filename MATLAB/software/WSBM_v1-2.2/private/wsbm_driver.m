function [Model] = wsbm_driver(Raw_Data,R_Struct,varargin)
% See 'help wsbm.m' or 'type wsbm.m' for information
%-------------------------------------------------------------------------%

% WSBM_Driver
% Version 1.0 | December 2013 | Christopher Aicher
%
%   Copyright 2013-2014 Christopher Aicher
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>

%-------------------------------------------------------------------------%
% WSBM_DRIVER CODE
%-------------------------------------------------------------------------%
W_Distr = [];
E_Distr = [];
Data = [];
Options = [];
% Parse and Setup Input 
format_input()
% 3 Cases:
% -Case 1: Once From Seed
if ~isempty(Options.seed),
    check_seed(Options.seed);
    % Run Inference
    [Para,Flags] = main_alg(Data,W_Distr,E_Distr,R_Struct,Options.seed,Options);
    % Score Inference
    Para.LogEvidence = calc_logEvidence(Data,W_Distr,E_Distr,R_Struct,Para,Options);
    
    % Print Info
    if Options.verbosity > 0, 
        fprintf('Best Result LogEvidence: %2.4e  \n',Para.LogEvidence);
    end
    
% -Case 2: Trials in Serial
elseif ~Options.parallel,
    if Options.verbosity > 0,
        fprintf('Running %u Trials in Series \n', Options.numTrials);
    end
    LogEvidenceBest = -Inf;
    %Begin Trial Loop
    for trial_num = 1:Options.numTrials,     
        % Set Seed
        seed = random_seed();
        % Run Inference
        [Para,Flags] = main_alg(Data,W_Distr,E_Distr,R_Struct,seed,Options);
        % Score Inference
        Para.LogEvidence = calc_logEvidence(Data,W_Distr,E_Distr,R_Struct,Para,Options);
        % Print Info
        if Options.verbosity > 0, 
            fprintf('Trial %u of %u : %2.4e | Best: %2.4e \n',...
                trial_num,Options.numTrials,Para.LogEvidence,LogEvidenceBest);
        end
        % Keep Track of Best Trial
        if Para.LogEvidence > LogEvidenceBest,
            ParaBest = Para;
            FlagsBest = Flags;
            seedBest = seed;
            LogEvidenceBest = Para.LogEvidence;
            if Options.verbosity > 0, % Indicate Better Model Found
                fprintf(' ** \n');
            end
        end
    end % End Trial Loop
    % Save Result
    Para = ParaBest;
    Flags = FlagsBest;
    Options.seed = seedBest;
    
% -Case 3: Trials in Parallel
else
    if Options.verbosity > 0,
        fprintf('Running %u Trials in Parallel \n', Options.numTrials);
    end
    % Setup
    seeds = cell(Options.numTrials,1);
    Paras = cell(Options.numTrials,1);
    Flagss = cell(Options.numTrials,1);
    for trial_num = 1:Options.numTrials,
        % Set Seeds
        seeds{trial_num} = random_seed();
    end
    % Begin Para Trial Loop
    parfor trial_num = 1:Options.numTrials,
        % Run Inference
        [Para,Flags] = main_alg(Data,W_Distr,E_Distr,R_Struct,seeds{trial_num},Options);
        % Score Inference
        Para.LogEvidence = calc_logEvidence(Data,W_Distr,E_Distr,R_Struct,Para,Options);
        % Save Temp Results
        Paras{trial_num} = Para;
        Flagss{trial_num} = Flags;
        % Print Info
        if Options.verbosity > 0,
            fprintf('Trial %u complete \n',trial_num);
        end
    end % End Para Trial Loop
    % Find Best Trial
    LogEvidenceBest = zeros(Options.numTrials,1);
    for trial_num = 1:Options.numTrials,
        LogEvidenceBest(trial_num) = Paras{trial_num}.LogEvidence;
    end
    [~,bestTrial] = max(LogEvidenceBest);
    if Options.verbosity > 0,
        fprintf('Best Result LogEvidence: %2.4e  \n',LogEvidenceBest(bestTrial));
    end
    % Save Result
    Para = Paras{bestTrial};
    Flags = Flagss{bestTrial};
    Options.seed = seeds{bestTrial};
end % End Cases
% Save Model Struct
Model = struct('name',[W_Distr.name,'-',E_Distr.name,'-',R_Struct.name],...
               'Data',Data,...
               'W_Distr',W_Distr,...
               'E_Distr',E_Distr,...
               'R_Struct',R_Struct,...
               'Para',Para,...
               'Flags',Flags,...
               'Options',Options);
if Options.verbosity > 0,
    fprintf('... wsbm.m Done\n');
end
           
%-------------------------------------------------------------------------%
% NESTED FUNCTIONS
%-------------------------------------------------------------------------%

% FORMAT_INPUT FUNC -------------------------------------------------------
function format_input()
    % Parse Input
    % Setup + Check Options
    parse_varagin();
    % Setup + Check W_Distr
    setup_W_Distr();
    % Setup + Check E_Distr
    setup_E_Distr();
    % Check for BP cases
    if strcmpi(Options.algType,'bp'),
        if Options.mexFile == 0,
            error('BP only works with MEX files. Consider setting algType = vb');
        end
        if isempty(E_Distr.tau_0),
            error('BP only works with E_Distr nonempty. Consider setting alpha = 0');
        end
        if strncmpi(E_Distr.name,'dc',2),
            error('BP + Degree Correction has not been implemented.');
        end
    end

    % Setup + Check R_Struct
    setup_R_Struct();
    % Setup Data
    if Options.mexFile,
        setup_Data(); % With Mex files
    else
        fprintf('Running Code without MEX files.\n');
        fprintf('Installing MEX files leads to significantly\n');
        fprintf('faster performance and better scaling...\n');
        if strncmpi(W_Distr.name,'DC',2) || strncmpi(E_Distr.name,'DC',2)
            error('Degree Corrected Models have only been implemented with MEX files.');
        end
        setup_Data_NoMex(); % Without Mex files
    end
    % Setup Mu_0
    if isempty(Options.mu_0),
        Options.mu_0 = ones(R_Struct.k,Data.n)/R_Struct.k;
    elseif size(Options.mu_0,1) ~= R_Struct.k ||...
            size(Options.mu_0,2) ~= Data.n,
        error('User specified mu_0 is not a k by n matrix');
    end     
    % Setup Save Directory
    if Options.save,
        if ~exist(Options.outputPath,'dir'),
            disp('Making Save Dir');
            mkdir(Options.outputPath);
        else
            disp('Adding/Overwriting Save Dir');
        end
    end
    % Setup Parallel Computing
    if Options.parallel && matlabpool('size') == 0,
        try
            matlabpool;
        catch err,
            fprintf('Not Running in Parallel: %s\n', err.message);
            Options.parallel = 0;
        end
    end    
end

% PARSE_INPUT FUNC --------------------------------------------------------
function parse_varagin()
    % Setup Default Options Struct
    Options = struct('algType','vb',...
                     'alpha',0.5,...
                     'networkType','directed',...
                     'nanType','missing',...
                     'parallel',0,...
                     'save',0,...
                     'outputPath',[cd,filesep,'WSBM_Temp_Output_MAT'],...
                     'verbosity',1,...
                     'numTrials',50,...
                     'mainMaxIter',80,...
                     'mainTol',0.001,...
                     'muMaxIter',50,...
                     'muTol',0.001,...
                     'mexFile',1,...
                     'seed',[],...
                     'mu_0',[]);
    % Parse Varagin
    for ii = 1:2:length(varargin)-1,
        argok = 1;
        if ischar(varargin{ii}),
            switch lower(varargin{ii}),
                case 'w_distr',
                    W_Distr = varargin{ii+1};
                case 'e_distr',
                    E_Distr = varargin{ii+1};
                case 'algtype',
                    Options.algType = varargin{ii+1};
                case 'alpha',
                    Options.alpha = varargin{ii+1};
                case 'networktype',
                    if ischar(varargin{ii+1}),
                        Options.networkType = varargin{ii+1};
                        if ~strcmpi(Options.networkType,'directed'),
                            fprintf('!WSBM.m currently only runs directed models!\n');
                            error('Use networkType == directed');
                        end
                    else
                        error('networkType must be a string (in varargin %u)',ii+1);
                    end
                case 'nantype',
                    if ischar(varargin{ii+1}),
                        Options.nanType = varargin{ii+1};
                        if ~strcmpi(Options.nanType,'missing'),
                            fprintf('!WSBM.m currently only allows nans to be missing!\n');
                            error('Use nanType == missing');
                        end
                    else
                        error('nanType must be a string (in varargin %u)',ii+1);
                    end
                case 'parallel',
                    Options.parallel = varargin{ii+1};
                case 'save',
                    Options.save = varargin{ii+1};
                case 'outputpath',
                    if ischar(varargin{ii+1}),
                        Options.outputpath = varargin{ii+1};
                    else
                        error('outputpath must be a string (in varargin %u)',ii+1);
                    end
                case 'verbosity',
                    Options.verbosity = varargin{ii+1};
                case 'seed',
                    Options.seed = varargin{ii+1};
                case 'numtrials',
                    Options.numTrials = varargin{ii+1};
                case 'mainmaxiter',
                    if isnumeric(varargin{ii+1}),
                        Options.mainMaxIter = varargin{ii+1};
                    else
                        error('Invalid mainMaxIter parameter in varargin %u',ii+1);
                    end
                case 'maintol',
                    if isnumeric(varargin{ii+1}),
                        Options.mainTol = varargin{ii+1};
                    else
                        error('Invalid mainTol parameter in varargin %u',ii+1);
                    end
                case 'mumaxiter',
                    if isnumeric(varargin{ii+1}),
                        Options.muMaxIter = varargin{ii+1};
                    else
                        error('Invalid muMaxIter parameter in varargin %u',ii+1);
                    end
                case 'mutol',
                    if isnumeric(varargin{ii+1}),
                        Options.muTol = varargin{ii+1};
                    else
                        error('Invalid muTol parameter in varargin %u',ii+1);
                    end
                case 'mexfile',
                    if isnumeric(varargin{ii+1}),
                        Options.mexFile = varargin{ii+1};
                    else
                        error('Invalid mexFile parameter in varargin %u',ii+1);
                    end
                case 'mu_0',
                    Options.mu_0 = varargin{ii+1};
                otherwise,
                    argok = 0;
            end
        else
            error(['Invalid argument #',num2str(ii)]);
        end
        if ~argok,
            error('Unknown argument #%u: %s', ii,varargin{ii});
        end  
    end % End Varagin For Loop
    
    % Check if mexFiles exist + compiled
    if Options.mexFile == 1,
        if ~(exist('calc_T_e_bra','file') == 3 &&...
             exist('calc_T_w_bra','file') == 3 &&...
             exist('vb_wsbm','file') == 3 &&...
             exist('create_T_e','file') == 3 && ...
             exist('create_T_w','file') == 3 && ...
             exist('create_T_bp','file') == 3 && ...
             exist('bp_wsbm','file') == 3),
            fprintf('Required MEX files not found... see InstallMEXFiles.m\n');
            fprintf('Running with slower O(n^2) MATLAB code instead of O(m)...\n');
            Options.mexFile = 0;
        end
    end
    
    % Check verbosity
    if Options.verbosity > 4,
        % Comment this out if you know what you are doing.
        error('Note: Options.verbosity > 4 is only used for debugging purposes.\n');
    end
    
end

% SETUP_W_DISTR FUNC ------------------------------------------------------
function setup_W_Distr()
    % Setup W_Distr Struct
    if isempty(W_Distr),
        W_Distr = setup_distr('normal');
        if Options.verbosity > 0,
            fprintf('W_Distr set to Normal (default)\n');
        end
    elseif ischar(W_Distr),
        W_Distr = setup_distr(W_Distr);
    elseif ~isstruct(W_Distr)
        error('Unrecognized W_Distr');
    end        
    % Check W_Distr Struct
    check_distr(W_Distr);
end

% SETUP_E_DISTR FUNC ------------------------------------------------------
function setup_E_Distr()
    % Setup E_Distr Struct
    if isempty(E_Distr),
        E_Distr = setup_distr('bernoulli');
        if Options.verbosity > 0,
            fprintf('E_Distr set to Bernoulli (default)\n');
        end
    elseif ischar(E_Distr),
        E_Distr = setup_distr(E_Distr);
    elseif ~isstruct(E_Distr)
        error('Unrecognized E_Distr');
    end        
    % Check E_Distr Struct
    check_distr(E_Distr);
end

% SETUP_R_Struct FUNC -----------------------------------------------------
function setup_R_Struct()
% Setup R_Struct 
    if ~exist('R_Struct','var'),
        % Default
        R_Struct = SBM_Struct(4);
        if Options.verbosity > 0,
            fprintf('R_Struct set to 4 blocks (default)     \n');
        end
    elseif isnumeric(R_Struct),
        if numel(R_Struct) == 1,
            % Block Model
            R_Struct = SBM_Struct(R_Struct);
        else
            % Custom Matrix
            newR_Struct.R = R_Struct;
            newR_Struct.k = length(newR_Struct.R);
            newR_Struct.r = max(newR_Struct.R(:));
            newR_Struct.name = sprintf('Custom%u',newR_Struct.k);
            R_Struct = newR_Struct;
        end
    elseif ~isstruct(R_Struct),
        error('Unrecognized R_Struct');
    end
    % Check R_Struct
    check_r_struct(R_Struct);
    % Nested Functions:
    % SBM_Struct
    function [out] = SBM_Struct(k)
        if strcmpi(Options.networkType,'directed')
            r = 0;
            R = zeros(k);
            for ii = 1:k
                for jj = 1:k
                    r = r+1;
                    R(ii,jj) = r;
                end
            end
        else
            error('networkType = %s has not been implemented yet',Options.networkType);
        end
        out.R = R;
        out.k = k;
        out.r = r;
        out.name = sprintf('SBM%u',k);
    end
end
% SETUP_DATA_NOMEX FUNC ---------------------------------------------------
function setup_Data_NoMex()
%SETUP_DATA_NOMEX creates/formats the Data struct used in WSBM when MeX
%files are not loaded/complied.
    
    % Convert E into an Adj_Matrix
    [n,m] = size(Raw_Data);
    if (n ~= m),
        if m == 3 && exist('Edg2Adj','file') == 2,
            fprintf('Assuming Raw Data is an Edge List\n');
            fprintf('Converting Raw Data to an Adjacency Matrix\n');
            Adj_Mat = Edg2Adj(Raw_Data);
            n = size(Adj_Mat,1);
        else
            error('Raw Data is not a square n by n Adjancency Matrix');
        end
    else
        fprintf('Treating Raw Data as an Adjacency Matrix\n');
        Adj_Mat = Raw_Data;
    end
    
    % Create Sufficient Statistics
    T_w = cell(numel(W_Distr.T),1);
    for tw = 1:numel(W_Distr.T),
        temp = W_Distr.T{tw}(Adj_Mat);
        temp(isnan(temp)) = 0;
        T_w{tw} = temp;
    end
    T_e = cell(numel(E_Distr.T),1);
    for te = 1:numel(E_Distr.T),
        T_e{te} = E_Distr.T{te}(~isnan(Adj_Mat)*1);
    end
            
    % Calculate Additive LogLikelihood Constants
    logHw = sum(W_Distr.logh(Adj_Mat(~isnan(Adj_Mat))));
    logHe = sum(E_Distr.logh(isnan(Adj_Mat(:))));
    % COMMENT: This would need to be changed for symmetric / non-Missing cases 
    
    % Create the Struct
    Data = struct('n',n,'Raw_Data',Raw_Data,'T_w',{T_w},'T_e',{T_e},...
    'logHw',logHw,'logHe',logHe);
end
% SETUP_DATA FUNC ---------------------------------------------------------
function setup_Data()
%SETUP_DATA creates/formats the Data struct used in WSBM 
    
    % Convert E into an Edge_List
    [m,s1] = size(Raw_Data);
    if (s1 ~= 3),
        if s1 == m && exist('Adj2Edg','file') == 2,
            fprintf('Assuming Raw Data is an Adjacency Matrix\n');
            fprintf('Converting Raw Data to an Edge List\n');
            Edge_List = Adj2Edg(Raw_Data);
            [m,~] = size(Edge_List);
        else
            error('Raw Data is not an m by 3 Edge List\n');
        end
    else
        fprintf('Treating Raw Data as an Edge List\n');
        Edge_List = Raw_Data;
    end
    Edge_List = sortrows(Edge_List,[1,2]);
    n = max(max(Edge_List (:,1:2)));
    
    % Create Sufficient Statistics 
    T_w = zeros(m,numel(W_Distr.T)+2);
    T_w(:,1:2) = Edge_List(:,1:2);
    if strncmpi(W_Distr.name,'DC',2) 
        % Handle Degree Corrected Case
        for tw = 1:numel(W_Distr.T)-1,
            T_w(:,tw+2) = W_Distr.T{tw}(Edge_List(:,3));
        end
        A = Edg2Adj(Edge_List);
        T_w(:,end) = W_Distr.T{end}(Edge_List,[nansum(A,2)';nansum(A,1)]);
    else
        % Handle General Corrected Case
        for tw = 1:numel(W_Distr.T),
            T_w(:,tw+2) = W_Distr.T{tw}(Edge_List(:,3));
        end
    end
    T_w = T_w'; % Convert to t_w+2 by m_E
    T_e = zeros(m,numel(E_Distr.T)+1);
    T_e(:,1:2) = Edge_List(:,1:2);
    edge_temp = ones(m,1);
    edge_temp(isnan(Edge_List(:,3))) = NaN;
    for te = 1:numel(E_Distr.T)-1,
        T_e(:,te+2) = E_Distr.T{te}(edge_temp);
    end
    T_e = T_e'; % Convert to t_e+1 by m_E
    
    % Calculate In- Out- Degrees
    % degrees(1,:) are the In Degrees
    % degrees(2,:) are the Out Degrees
    degrees_total = zeros(2,n); % Counts missing / NaN edges
    degrees_w = zeros(2,n); % Ignores missing edges
    for ee = 1:m,
        degrees_total(2,Edge_List(ee,1)) = degrees_total(2,Edge_List(ee,1))+1;
        degrees_total(1,Edge_List(ee,2)) = degrees_total(1,Edge_List(ee,2))+1;
        if ~isnan(Edge_List(ee,3)),
            degrees_w(2,Edge_List(ee,1)) = degrees_w(2,Edge_List(ee,1))+1;
            degrees_w(1,Edge_List(ee,2)) = degrees_w(1,Edge_List(ee,2))+1;
        end
    end
    
    % Create Sufficient Statistics Cell Arrays
    if Options.verbosity > 0,
        fprintf('Setting Up Weighted Statistics...\n');
    end
    if numel(W_Distr.T) > 0,
        [T_w_in,T_w_out] = create_T_w(T_w,degrees_w);
    else
        T_w_in = {};
        T_w_out = {};
    end
    if Options.verbosity > 0,
        fprintf('Setting Up Edge Statistics...\n');
    end
    if numel(E_Distr.T) > 0,
        [T_e_in,T_e_out] = create_T_e(T_e,degrees_total);
    else
        T_e_in = {};
        T_e_out = {};
    end
    if strcmpi(Options.algType,'bp'),
        if Options.verbosity > 0,
            fprintf('Setting Up BP Statistics...\n');
        end
        [T_w_in,T_w_out,T_e_in,T_e_out] = ...
            create_T_bp(T_w_in,T_w_out,T_e_in,T_e_out);
    end
    
    % Calculate Additive LogLikelihood Constants
    logHw = sum(W_Distr.logh(Edge_List(:,3)));
    logHe = m*E_Distr.logh(1); % Not exactly correct for Multigraphs
    logHe = logHe+(n*(n-1)-m)*E_Distr.logh(0);
    % COMMENT: This would need to be changed for symmetric / non-Missing cases 
    
    % Create the Struct
Data = struct('n',n,'Raw_Data',Raw_Data,'T_w',T_w,'T_e',T_e,'degrees_w',degrees_w,...
    'degrees_total',degrees_total,'T_w_in',{T_w_in},'T_w_out',{T_w_out},...
    'T_e_in',{T_e_in},'T_e_out',{T_e_out},'logHw',logHw,'logHe',logHe);
end

% RANDOM_SEED  ------------------------------------------------------------
function theSeed = random_seed(tuning)
%RANDOM_SEED returns an appropriate random seed for the WSBM.m inference
% tuning is a tuning parameter which biases the seed to concentrate mass
%    tuning = 1 -> random/flat, 5 -> default (0 is a bad idea).
if nargin < 1, tuning = 5; end;
mu_seed = rand(R_Struct.k,Data.n).*Options.mu_0;
mu_seed = (mu_seed./(ones(R_Struct.k,1)*max(mu_seed,[],1))).^tuning;
theSeed = mu_seed./(ones(R_Struct.k,1)*sum(mu_seed,1));            
end

% CHECK_SEED --------------------------------------------------------------
function check_seed(theSeed)
%CHECK_SEED checks to make sure seed is properly formatted. 
if strncmpi(Options.algType,'vb',2),
    [k,n] = size(theSeed);
    %Check VB Seed Format (mu)
    if n ~= Data.n || k ~= R_Struct.k,
        error('Invalid Mu Seed Format:\nSeed is %u by %u needs to be k by n',n,k);
    end
elseif strncmpi(Options.algType,'bp',2),
    error('BP SEED CHECKER needs to be fixed');
    %Check BP Seed Format (mu,mes)
    if  ~isstruct(theSeed),
        error('For BP, Seed needs to be struct of mu and mes initial values');
    end
    [n,m] = size(theSeed.mu);
    if n ~= Data.n || m ~= R_Struct.k,
        error('Invalid Mu Seed Format:\nSeed is %u by %u needs to be n by k',n,m);
    end
    [n,m] = size(theSeed.mes);
    if n ~= Data.m || m~= R_Struct.k,
        error('Invalide Mes Seed Format:\nSeed is %u by %u needs to be m by k',n,m);
    end
else
    error('Unrecognized algType: %s',Options.algType);
end    
end

%-------------------------------------------------------------------------%
end % End of WSBM_DRIVER.M
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% HELPER FUNCTIONS
%-------------------------------------------------------------------------%

% CHECK DISTR STRUCT ------------------------------------------------------
function [] = check_distr(Distr)
%CHECK_DISTR checks to make sure a distr struct has the necessary fields
%for WSBM.m
%See `help SETUP_DISTR.m' for more details
    if ~isstruct(Distr),
        error('Distr needs to be a Distr Struct');
    end
    if ~isfield(Distr,'tau_0'),
        error('Distr is missing the field tau_0');
    end
    if ~isfield(Distr,'logh'),
        error('Distr is missing the field logh');
    end
    if ~isfield(Distr,'T'),
        error('Distr is missing the field T');
    end
    if size(Distr.tau_0,2) ~= size(Distr.T,1),
        error('size(tau_0,2) is not dimT in Distr');
    end
    if ~isfield(Distr,'Eta'),
        error('Distr is missing the field Eta');
    end
    if size(Distr.Eta,1) ~= size(Distr.T,1),
        error('size(Eta,1) is not size(T,1) in Distr');
    end
    if ~isfield(Distr,'logZ'),
        error('Distr is missing the field logZ');
    end
    if ~isfield(Distr,'Theta'),
        error('Distr is missing the field Theta');
    end
    if ~isfield(Distr,'Predict'),
        error('Distr is missing the field Predict');
    end
    if ~isfield(Distr,'name'),
        error('Distr is missing the field name');
    end
end

% CHECK R_Struct ----------------------------------------------------------
function [] = check_r_struct(R_Struct)
%CHECK_R_STRUCT checks to make sure R_list has proper format. Throws errors 
% if R_Struct does not have proper format.
    if ~isstruct(R_Struct),
        error('R_Struct needs to be a struct');
    end
    if ~isfield(R_Struct,'R')
        error('R_Struct is missing the field R');
    end
    if size(R_Struct.R,1) ~= size(R_Struct.R,2),
        error('R is not square in R_Struct');
    end
    if ~isfield(R_Struct,'r')
        error('R_Struct is missing the field r');
    end
    if R_Struct.r ~= max(R_Struct.R(:)),
        error('R and r do not agree on the number of edge bundles\n r ~= max(R) in R_Struct');
    end
    if ~isfield(R_Struct,'k')
        error('R_Struct is missing the field k');
    end
    if size(R_Struct.R,1) ~= R_Struct.k,
        error('R is not a kxk matrix in R_Struct');
    end
    if ~isfield(R_Struct,'name')
        error('R_Struct is missing the field name');
    end
end 

% EOF