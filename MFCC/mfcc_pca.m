function [A,f]=mfcc1(s,fs)
b=[1 -0.98];
a=1;
s=filter(b,a,s);
% s=s./abs(max(s));
% s=s-mean(s); 
n=480;    % n=input('Enter the frame size');
ft=480/fs; fst=160/fs;
[r,cl]=size(s);
f=round((r-480)/160);
j=1;
k=n;
h=hamming(k);
for i=1:f
    d(:,i)=s(j:k);
    e(:,i)=d(:,i).*h;
    j=j+160;
    k=k+160;
    i=i+1;
end
  ft=fft(e);
  p=20;     % p=input('Enter the number of filter banks');
  m=melfb(p,n,fs);
  n2=1+floor(n/2);
  for i=1:f
    t=ft(:,i);
     z(:,i)=m*abs(t(1:n2)).^2;
     i=i+1;
  end
  v=dct(log(z));
  c=v';
  A= c(:,2:17);
  A = A(isfinite(A(:, 1)), :); % in order to check if there are any NaN values in the A matrix
  A= pca(A);
