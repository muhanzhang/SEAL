% A script used to transfrom txt nets to mat formats.

dataname = strvcat('nettest', 'USAir','NS','PB','Yeast','Celegans','Power','Router');

for ith_data = 1:size(dataname, 1)                           
    tempcont = ['processing the ', int2str(ith_data), 'th dataset...', dataname(ith_data,:)];
    disp(tempcont);
    thisdatapath = strcat('data/raw_data/', dataname(ith_data,:),'.txt');    
    linklist = load(thisdatapath);                                 
    net = FormNet(linklist); clear linklist; 
    name = strcat('data/', dataname(ith_data, :), '.mat');
    save(name, 'net');
end
