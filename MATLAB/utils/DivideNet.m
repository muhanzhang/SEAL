function [train test] = DivideNet(net, ratioTrain, connected)
%%划分训练集和测试集，保证训练集连通
    net = triu(net) - diag(diag(net));  % convert to upper triangular matrix
    num_testlinks = floor((1-ratioTrain) * nnz(net));      
    % 确定测试集的边数目
    [xindex, yindex] = find(net);  linklist = [xindex yindex];    
    % 将网络（邻接矩阵）中所有的边找出来，存入linklist  
    clear xindex yindex;  
    % 为每条边设置标志位，判断是否能删除
    test = sparse(size(net,1),size(net,2));                 
    while (nnz(test) < num_testlinks)               %For power dataset, maximum 636 test links <660 expected. 
        if length(linklist) <= 2
            break;
        end
        %---- 随机选择一条边
        index_link = ceil(rand(1) * length(linklist));
        
        uid1 = linklist(index_link,1); 
        uid2 = linklist(index_link,2);    
        net(uid1,uid2) = 0;  
        
        
          %% 
        %---- 判断所选边两端节点uid1和uid2是否可达，若可达则可放入测试集，否则重新挑选一条边
         
        % 将这条边从网络中挖去用以判断挖掉后的网络是否还连通
        tempvector = net(uid1,:);
        % 取出uid1一步可达的点，构建成一维向量
        sign = 0;  
        % 标记此边是否可以被移除，sign=0表示不可； sign=1表示可以
        uid1TOuid2 = tempvector * net + tempvector;        
        % uid1TOuid2表示二步内可达的点
        if uid1TOuid2(uid2) > 0
            sign = 1;               
            % 二步即可达
        else
            while (nnz(spones(uid1TOuid2) - tempvector) ~=0)   
            % 直到可达的点到达稳定状态，仍然不能到达uid2，此边就不能被删除
                tempvector = spones(uid1TOuid2);
                uid1TOuid2 = tempvector * net + tempvector;    
                % 此步的uid1TOuid2表示K步内可达的点
                if uid1TOuid2(uid2) > 0
                    sign = 1;      
                     % 某步内可达
                    break;
                end
            end
        end 
        % 结束-判断uid1是否可达uid2
        
        if connected == false
            sign = 1;  % overwrite, keep all selected links in test, no matter whether the remaining net is connected
        end

        %% 
        
        %----若此边可删除，则将之放入测试集中，并将此边从linklist中移除
        if sign == 1 %此边可以删除
            linklist(index_link,:) = []; 
            test(uid1,uid2) = 1;
        else
            linklist(index_link,:) = [];
            net(uid1,uid2) = 1;   
        end   
        % 结束-判断此边是否可以删除并作相应处理
    end   
    % 结束（while）-测试集中的边选取完毕
    train = net + net';  test = test + test';
    % 返回为训练集和测试集
end
