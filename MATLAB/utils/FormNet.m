function [ net ] = FormNet( linklist )
%% 读入连边列表linklist，构建网络邻接矩阵net
    %---- 如果节点编号从0开始，将所有节点编号加1（matlab的下标从1开始）
    if ~all(all(linklist(:,1:2)))
        linklist(:,1:2) = linklist(:,1:2)+1;
    end
    
    %----对无向图，将第三列元素置为1
    linklist(:,3) = 1;
    net = spconvert(linklist);    %
    nodenum = length(net);
    net(nodenum,nodenum) = 0;                               
    % 此处删除自环，对角元为0以保证为方阵
    net = net-diag(diag(net));
    net = spones(net + net'); 
    % 确保邻接矩阵为对称矩阵，即对应于无向网络  %in case of dataset noise
end 
% 转换过程结束，得到网络的邻接矩阵
