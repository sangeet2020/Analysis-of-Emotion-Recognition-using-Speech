% function [A,f]=mfcc1(s,fs)
close all;clear all;clc;
[s,fs] = audioread('03a02Wb.wav');
% t=0;
% for r=0:9;
% for p=0:9;    
% file=sprintf('%s%d%d%d','D:\LPBF_Mar09\DB24_digit_II\0\',t,r,p);
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
  %disp(w);
  
%  for h=1:f
%  fid = fopen('mfcc_0_new.m','a');
%  count = fprintf(fid,'%f',w(h,:));
%  count=fprintf(fid,'\n');
%  end
%  end
% end
% u=w';
% 
%    r = vqlbg(u, 1);
%    disp(r);
% %     c=c(1:13);
% %     c=kmeans(v',1);
%  q=r';
%    w=(1:16);
% end
%    plot(linspace(0, (16000/2), 129), melfb(24, 256, 16000)'),
%   title('Mel-spaced filterbank'), xlabel('Frequency =(Hz)');
