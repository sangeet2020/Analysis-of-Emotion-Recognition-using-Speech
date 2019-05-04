close all;
clear all;
clc;
y=load('sad_savee_mfcc_delta_doubledelta.m');

tr0=y(:,:);
%i1=i1+200;
% tr1=X(8:14,:); tr2=X(15:21,:); tr3=X(21:28,:); 
% tr0=tr0';tr1=tr1';tr2=tr2';tr3=tr3';
 X=tr0';
d=size(X,1); % Dimension of data vector
N=size(X,2); % Number of Data points
L=16;%codebook size
[CB, p, DistHist] = vqsplit(X,d,N,L);
% P: (Original) Weight of each cluster i.e. the number of its vectors divided by total
%       number of vectors; (Modified) Population of each cluster
M=size(CB,2);
T=size(X,2);

% [I, dst]=VQIndex(X,T,CB,M) ;
% 
%  sigm1 = variance_calculation(X, d, I, M);

[mu,sigm,c]=gmm_estimate(X,M,100,CB);
%to store the mu
% q=zeros(13,8);
% r=zeros(13,8);
for i=1:16
fid = fopen('savee60_gmm_model_mfcc_delta_doubledelta.m','a');
 fprintf(fid,'%f ',mu(i,:)); fprintf(fid,'\n');
 fclose(fid);
end

%to store the sigm
for i=1:16
fid = fopen('savee60_gmm_model_mfcc_delta_doubledelta.m','a');
 fprintf(fid,'%f ',sigm(i,:)); fprintf(fid,'\n');
 fclose(fid);
end
%to store mixing coefficient c
fid = fopen('savee60_gmm_model_mfcc_delta_doubledelta.m','a');
 fprintf(fid,'%f ',c'); fprintf(fid,'\n');
 fclose(fid);
