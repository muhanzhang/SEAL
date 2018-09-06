function [n] = card_3inter(L1, L2, L3, l1, l2, l3)

% count cardinalities of subsets L1&L2&L3 (n(7)), L1&L2\L3 (n(4)), L1&L3\L2 (n(5)), L2&L3\L1 (n(6)), 
% L1\(L2|L3) (n(1)), L2\(L1|L3) (n(2)), L3\(L1|L2) (n(3))
% Li - a set
% li - number of elements in Li
% Author: Nino Shervashidze - nino.shervashidze@tuebingen.mpg.de
% Copyright 2012 Nino Shervashidze

n=zeros(1,7);
i=1; j=1; k=1;

while i<=l1 && j <=l2 && k<=l3
  m=find_min(L1(i), L2(j), L3(k));
  n(m)=n(m)+1;
  switch m
   case 1
    i=i+1;
   case 2
    j=j+1;
   case 3
    k=k+1;
   case 4
    i=i+1; j=j+1;
   case 5
    i=i+1; k=k+1;
   case 6
    j=j+1; k=k+1;
   case 7
    i=i+1; j=j+1; k=k+1;
  end
end
%[i j k]
if i>l1 && j>l2 && k>l3
else 
  if i>l1 && j>l2
    n(3)=n(3)+l3-k+1; k=l3+1;
  else 
    if i>l1 && k>l3
      n(2)=n(2)+l2-j+1; j=l2+1;
    else 
      if j>l2 && k>l3
        n(1)=n(1)+l1-i+1; i=l1+1;
      else 
        if i>l1
          while j<=l2 && k<=l3
            if L2(j)<L3(k) n(2)=n(2)+1; j=j+1;	
            else 
              if L2(j)>L3(k) n(3)=n(3)+1; k=k+1;
              else n(6)=n(6)+1; j=j+1; k=k+1;
              end
            end
          end
        else 
          if j>l2
            while i<=l1 && k<=l3
              if L1(i)<L3(k) n(1)=n(1)+1; i=i+1;	
              else 
                if L1(i)>L3(k) n(3)=n(3)+1; k=k+1;
                else n(5)=n(5)+1; i=i+1; k=k+1;
                end
              end
            end
          else 
            if k>l3
              while i<=l1 && j<=l2
                if L1(i)<L2(j) n(1)=n(1)+1; i=i+1;	
                else
                  if L1(i)>L2(j) n(2)=n(2)+1; j=j+1;
                  else n(4)=n(4)+1; i=i+1; j=j+1;
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


if i>l1 && j>l2 && k>l3
else 
  if i>l1 && j>l2
    n(3)=n(3)+l3-k+1;
  else 
    if i>l1 && k>l3
      n(2)=n(2)+l2-j+1;
    else
      if j>l2 && k>l3
        n(1)=n(1)+l1-i+1;
      end
    end
  end
end

end

function [m]=find_min(a,b,c)
% auxiliary function: determine minimum (minima) from the set a,b,c 
% return 7 if a=b=c, 6 if b=c<a, 5 if a=c<b, 4 if a=b<c, 3 if c<a & c<b, 2 if b<a & b<c, 1 if a<b & a<c.

mini=a;
if b<mini mini=b; end
if c<mini mini=c; end

if mini==a
  if mini==b
    if mini==c
      m=7;
    else m=4;
    end
  else 
    if mini==c
      m=5;
    else m=1;
    end
  end
else 
  if mini==b
    if mini==c
      m=6;
    else m=2;
    end
  else m=3;
  end
end

end

