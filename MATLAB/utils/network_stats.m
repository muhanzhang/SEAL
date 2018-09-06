dataname = strvcat('USAir','NS','PB','Yeast','Celegans','Power','Router'); data_num = 7;
%dataname = strvcat('Power','Router'); data_num = 2;
%dataname = strvcat('USAir'); data_num = 1;
datapath = strcat(pwd,'/data/');       %数据集所在的路径
stats = zeros(data_num+1,5);
for ith_data = 1:data_num+1              %consider the meta data +1                           
    % 遍历每一个数据
    tempcont = strcat('正在处理第 ', int2str(ith_data), '个数据...');
    disp(tempcont);
    tic;
    if ith_data <= data_num
        thisdatapath = strcat(datapath,dataname(ith_data,:),'.txt');    % 第ith个数据的路径
        linklist = load(thisdatapath);                                  % 导入数据（边的list）
        net = FormNet(linklist); clear linklist;                       % 根据边的list构成邻接矩阵
    else
        load(strcat(datapath,'EcoliModel.mat'));
        S = iAF1260.S;
        net = S*S';   %Adjacency matrix
        net = net - diag(diag(net));
        net = logical(net);
    end
    [largest,components] = largestcomponent(net);
    
    net1 = triu(net);
    N = size(net1,1);    %nodes number
    M = nnz(net1);     %links number
    average_degree = nnz(net)/N;
    stats(ith_data,:) = [N,M,largest,components,average_degree];    %nodes number, links number, ..., ..., 
    
    subplot(4,2,ith_data);
    vet = sum(net,2);
    histogram(vet);
    
   
    
end
