% v=wavread('H:\database\voice1\girl1\sad sthity1'); 
 function [l]=inpt(v)
 m=max(v)/50;
%  subplot(2,1,1);
%   plot(v); 
j=1;
c=0;
k=0;
for i=1:length(v)
    k=k+1;
    if v(i) > m
        break;
    end
end
for i=length(v):-1:k
    c=i;
    if v(i)>m
        break;
    end
end
for i=k:c
    l(j,1)=v(i);
        j=j+1;
    end        
%   subplot(2,1,2);
%   plot(l);   
  
