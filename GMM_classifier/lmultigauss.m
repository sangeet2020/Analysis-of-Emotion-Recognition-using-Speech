 function [YM,Y]=lmultigauss(x,mus,sigm,c)
% [lYM,lY]=lmultigauss(X,mu,sigm,c)
% 
% computes multigaussian TOTAL log-likelihood
% 
% X   : (LxT) data (columnwise vectors)
% sigm: (LxM) variances vector  (diagonal of the covariance matrix)
% mu  : (LxM) means
% c   : (Mx1) the weights
%  tr0=load('f_lpc100test_all.m');
% x=tr0(21:30,:);
% x=x';
%  p=load('gmm_model.m');
%  mus=p(1:13,:)
%  sigm=p(14:26,:)
%  c=p(27,:)
%  c=c'
DEBUG=0;
DEBUG1=0;

[L,T]=size(x);
M=size(c,1);

if DEBUG [ size(x), size(mus), size(sigm), size(c)], end

% repeating, changing dimensions:
X1=permute(repmat(x',[1,1,M]),[1,3,2]);      % (T,L) -> (T,M,L) one per mixture

size_X1 = size(X1);

clear x;%Inserted by abk

Mu1=permute(repmat(mus,[1,1,T]),[3,2,1]);     % (L,M) -> (T,M,L)

clear mus;%Inserted by abk

if DEBUG size(X1), size(Mu1),pause; end

X1 = X1 - Mu1;

clear Mu1;%Inserted by abk

Sigm1=permute(repmat(sigm,[1,1,T]),[3,2,1]); % (L,M) -> (T,M,L)

%Mu1=permute(repmat(mus,[1,1,T]),[3,2,1]);     % (L,M) -> (T,M,L)

%clear mus;%Inserted by abk

if DEBUG size(Sigm1),pause; end

%Y=squeeze(exp( 0.5.*dot(X-Mu,(X-Mu)./Sigm))) % L dissapears: (L,T,M) -> (T,M)
%lY=-0.5.*dot_abk(X1-Mu1,(X1-Mu1)./Sigm1,3);   %                       -> (T,M)
lY=-0.5.*dot(X1,X1./Sigm1,3);   %                       -> (T,M)

clear Sigm1;%Inserted by abk

%size_lY = size(lY)

% c,const -> (T,M) and then multiply by old Y
lcoi=log(2.*pi).*(L./2)+0.5.*sum(log(sigm),1); % c,const -> (T,M)
lcoef=repmat(log(c')-lcoi,[T,1]);

%size_lcoef = size(lcoef)

clear c sigm;%Inserted by abk

if DEBUG1 lcoi,lcoef,lY,pause;end

clear lcoi;%Inserted by abk

YM=lcoef+lY;            % ( T,M ) one mixture per column

clear lcoef lY;%Inserted by abk

Y=lsum(YM,2)   ;            % add mixtures TOTAL LOG-LIKELIHOOD 

if DEBUG [ size(YM) NaN size(Y) ], end
  
