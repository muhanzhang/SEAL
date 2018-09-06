dataname = strvcat('USAir','NS','PB','Yeast','Celegans','FWFB','Power','Router'); 
datapath = 'D:\current works\writting book\data\';      %数据集所在的路径
thisdatapath = strcat(datapath,dataname(7,:),'.txt');
data = load(thisdatapath);
a = [data(:,1); data(:,2)];
b = unique(a);
m = max(a);
n = min(a);
[length(b) m n]
