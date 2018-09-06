a=[0,9,5,1;9,0,2,0;5,2,0,1;1,0,1,0];

[xx,yy]=ind2sub(size(sim),i);
s=full(sim);


%for test the optimal k in SPC
K=[9,20,50,200,500,800,900,1000,1100,1200,1300,1400,1500,1668];
res=[];
for i = 1:length(K)
    k=K(i)
    tempauc = SPC2(wtrain.^alpha, test,k)
    res=[res,tempauc];
end
K
res






% to test the random precision
n=1324;
for k=500
t=[];

    for i=1:100000
    K=randperm(n);
    t = [t,nnz(K(1:k)<=k)];
    end
mean(t)
k^2/n
end
