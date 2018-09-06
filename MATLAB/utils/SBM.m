% wrapper for SBM, you need to install WSBM software into ./software

function [auc, precision] = SBM(train, test, group)
    %addpath('software/WSBM_v1-2.2');
    cd software/WSBM_v1-2.2;
    %InstallMEXFiles
    [i, j] = find(train);
    Train_List = [i, j, ones(length(i), 1)];
    maxi = size(train, 1);
    Train_List = [Train_List; maxi, maxi, 0];  % avoid missing some nodes in Model
    [I, J] = find(ones(size(train)));  % give scores to all possible links
    Test_List = [I, J, ones(length(I), 1)];
    [Labels, Model] = wsbm(Train_List, group, 'alpha', 1);
    Mu_I = Model.Para.mu(:,I);
    Mu_J = Model.Para.mu(:,J);
    edge_Predict = reshape(Model.E_Distr.Predict(...
        Model.Para.theta_e(Model.R_Struct.R(:),:)),...
        size(Model.R_Struct.R,1),size(Model.R_Struct.R,2));
    edge_Predict(isnan(edge_Predict)) = 0;
    edge_list = sum(Mu_I'*edge_Predict.*Mu_J',2);
    cd ../..;
    sim = zeros(size(train));
    sim(sub2ind(size(sim), I, J)) = edge_list;
    [auc, precision] = CalcAUC(train,test,sim);   
   
end
