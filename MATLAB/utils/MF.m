function [ auc, precision ] = MF(train, test, k, ith_experiment, SGD)
%  Usage: link prediction using matrix factorization, need libFM
%  --Input--
%  -train: the training matrix, positive links are 1, otherwise 0
%  -test: the testing matrix, positive links are 1, otherwise 0
%  -k: number of latent factors in matrix factorization
%  -SGD: if SGD = 1, will use adaptive sgd for optimization (which can 
%        output parameters w, v). Otherwise, will use MCMC (default).
%  --Output--
%  -auc: the AUC score
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin < 3
    k = 8;
end

if nargin < 4
    ith_experiment = 1;
end

if nargin < 5
    SGD = 0;
end

if exist('tempdata/libFM') == 0
    !cp software/libfm-1.42.src/bin/libFM tempdata/
end

FMtrain = sprintf('FMtrain_exp%d', ith_experiment);
FMtest = sprintf('FMtest_exp%d', ith_experiment);
FMconvert(train, FMtrain, 1);  % call FMconvert.m to convert train, all to libfm format
delete(strcat('tempdata/', FMtrain));
all = ones(size(train));  % to predict all the possible links' scores
RC = FMconvert(all, FMtest);
delete(strcat('tempdata/', FMtest));

switch SGD
    case 0
        % to use MCMC optimization
        cd tempdata;
        cmd = sprintf('./libFM -task r -train %s.libfm -test %s.libfm -dim "1,1,%d" -out out_%d.txt -iter 100 -ith_experiment "%d"', FMtrain, FMtest, k, ith_experiment, ith_experiment);  % change to -task r for a regression loss
        system(cmd);    %run libFM
        pred = dlmread(sprintf('out_%d.txt', ith_experiment));    %load the output file of libFM
        delete(sprintf('%s.libfm', FMtrain));
        delete(sprintf('%s.libfm', FMtest));
        delete(sprintf('out_%d.txt', ith_experiment));
        cd ..;
    case 1
        % to use adaptive sgd with separate validation set
        [r,c,v] = find(train);
        rperm = randperm(length(r));
        val = rperm(1: floor(0.1 * length(r)));    % split the train and validation data
        tra = rperm(floor(0.1 * length(r)) + 1: end);
        validation = zeros(size(train));        % reassemble the matrices
        validation(sub2ind(size(validation),r(val),c(val))) = v(val);
        train = zeros(size(train));
        train(sub2ind(size(train),r(tra),c(tra))) = v(tra);
        FMconvert(train, FMtrain);
        FMvalidation = sprintf('FMvalidation_exp%d', ith_experiment);
        FMconvert(validation, FMvalidation);
        cd tempdata;
        cmd = sprintf('./libFM -task r -train %s.libfm -test %s.libfm -validation %s.libfm -dim "1,1,%d" -out out_%d.txt -iter 100 -method sgda -learn_rate 0.001 -init_stdev 0.01', FMtrain, FMtest, FMvalidation, k, ith_experiment);
        system(cmd);                            % run libFM
        pred = dlmread(sprintf('out_%d.txt', ith_experiment));    % load the output file of libFM
        delete(sprintf('out_%d.txt', ith_experiment));
        delete(sprintf('%s.libfm', FMtrain));
        delete(sprintf('%s.libfm', FMtest));
        delete(sprintf('%s.libfm', FMvalidation));
        delete(FMtrain);
        delete(FMtest);
        delete(FMvalidation);
        cd ..;
end

%%
sim = zeros(size(train));
sim(sub2ind(size(sim), RC(:,1), RC(:,2))) = pred;    % assign FM predictions to sim
dsim = diag(diag(sim));     % diagonal elements of sim
sim = sparse((sim + sim'- dsim));  % NOTE: (sim + sim' - dsim) keeps diagonal elements

[auc, precision] = CalcAUC(train, test, sim);
end





function [RC] = FMconvert(train, libfmname, neg)
%  A script to convert matrix to libFM formats.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin < 3
    neg = 0;  % do not also sample negatives
end

datapath = strcat(pwd, '/tempdata/');       % path of the data folder
train = triu(train);
[r,c,v] = find(train);  % store row, col, value of train>0 data to libfmname, for further processing using python
RC = [r, c];
a = [v, r-1, c-1];

if neg == 1  % if neg = 1, sample the same number of negative links
    [r2, c2, v2] = find(train == 0);
    perms = randperm(length(r2));
    r2 = r2(perms(1: length(r)));
    c2 = c2(perms(1: length(c)));
    v2 = v2(perms(1: length(v)));
    RC2 = [r2, c2];
    a2 = [v2 - 1, r2 - 1, c2 - 1];
    a = [a; a2];
    RC = [RC; RC2];
end

process_method = 2;
switch process_method
case 1
% directly process text in MATLAB, slower
fid = fopen(sprintf('tempdata/%s.libfm', libfmname), 'w+');
for i = 1:size(a, 1)
    b = a(i, :);
    fprintf(fid, strcat(num2str(b(1)), 32, num2str(b(2)), ':1', 32, num2str(b(3)), ':1', '\r\n'));
end
fclose(fid);

case 2
%% call python to process text, faster
dlmwrite(strcat(datapath, libfmname), a, 'delimiter', ' ');
cd tempdata;
cmd = sprintf('python ../utils/processFM.py %s', libfmname);     % call python, convert to [v, r, c] to [v, r:1, c:1]
system(cmd);
cd ..;
end

end
