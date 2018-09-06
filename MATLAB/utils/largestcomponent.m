function [largest,components] = largestcomponent(A)

%A = triu(A);
n=length(A);

for i=1:n
    mz{i}=find(A(i,:));
end

x(1:n)=0;
z=0;
k=0;
for i=1:n
    if x(i)==0;
        z=z+1;
        clear v
        v(1)=i;
        while nnz(v)>0
                    x(v(1))=z;
                    k=union(k,v(1));
                    b=setdiff(mz{v(1)},k);
                    v=setdiff(v,v(1));
                    v=union(v,b);         
        end
    end
end

c(1:max(x))=0;
for i=1:max(x)
    c(i)=length(find(x==i)); 
end

largest = max(c);
components = size(c,2);

%cm=find(c==max(c));

%B=find(x==cm(1));


end
                    
        