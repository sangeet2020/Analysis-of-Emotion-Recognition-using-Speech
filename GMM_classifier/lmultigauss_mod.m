function [YM,Y]=lmultigauss_mod(x,mus,sigm,c)
% [YM,Y]=lmultigauss(X,mu,sigm,c)
% 
% computes multigaussian log-likelihood
% 
% X   : (LxT) data (columnwise vectors)
% sigm: (LxM) variances vector  (diagonal of the covariance matrix)
% mu  : (LxM) means
% c   : (Mx1) the weights
%
%Outputs:
%YM            % ( T,M ) one mixture per column
%Y             % add mixtures ---> (T,1)
%
  
DEBUG=0;
DEBUG1=0;

[L,T]=size(x);
M=size(c,1);

if DEBUG sz_x=size(x), sz_mus=size(mus), sz_sigm=size(sigm), sz_c=size(c), pause; end

max_vol = 300000000;%32203990 20000000
half_max_vol = floor(max_vol/2);
vol = M*L;
vol1 = T*vol;
vol2 = vol1;
if vol1 >= max_vol
    YM = []; Y = [];
    T1 = T;
    while vol2>max_vol %half_max_vol 
        T1 = floor(T1/2);
        vol2 = T1*vol;
    end
    ts=1;
    while (ts+T1-1) < T
        te = ts + T1 -1;
        x1 = x(:,ts:te);
        % repeating, changing dimensions:
        X1=permute(repmat(x1',[1,1,M]),[1,3,2]);      % (T1,L) -> (T1,M,L) one per mixture
    
        Mu1=permute(repmat(mus,[1,1,T1]),[3,2,1]);     % (L,M) -> (T1,M,L)

        X1 = X1 - Mu1;

        clear Mu1;%Inserted by abk

        Sigm1=permute(repmat(sigm,[1,1,T1]),[3,2,1]); % (L,M) -> (T1,M,L)

        %Y=squeeze(exp( 0.5.*dot(X-Mu,(X-Mu)./Sigm))) % L dissapears: (L,T1,M) -> (T1,M)
        %lY=-0.5.*dot_abk(X1-Mu1,(X1-Mu1)./Sigm1,3);   %                       -> (T1,M)
        lY=-0.5.*dot_abk(X1,X1./Sigm1,3);   %                       -> (T1,M)

        clear Sigm1;%Inserted by abk

        % c,const -> (T,M) and then multiply by old Y
        lcoi=log(2.*pi).*(L./2)+0.5.*sum(log(sigm),1);%  (L,M) -> (1,M) 
        lcoef=repmat(log(c')-lcoi,[T1,1]);% c' (1, M) -> (T1,M)
    
        clear lcoi;

        YMM=lcoef+lY;            % ( T1,M ) one mixture per column
        
        clear lcoef;

        YY=lsum(YMM,2);                 % add mixtures ---> (T1,1)
    
        YM = [YM ; YMM];
    
        Y = [Y ; YY];
    
        ts = ts + T1;
    
    end%while (ts+T1-1) <= T
    %Rest of the length T of x
    T1 = T - ts + 1;
    x1 = x(:,ts:T);
    % repeating, changing dimensions:
    X1=permute(repmat(x1',[1,1,M]),[1,3,2]);      % (T1,L) -> (T1,M,L) one per mixture
    
    Mu1=permute(repmat(mus,[1,1,T1]),[3,2,1]);     % (L,M) -> (T1,M,L)

    X1 = X1 - Mu1;

    clear Mu1;%Inserted by abk

    Sigm1=permute(repmat(sigm,[1,1,T1]),[3,2,1]); % (L,M) -> (T1,M,L)

    %Y=squeeze(exp( 0.5.*dot(X-Mu,(X-Mu)./Sigm))) % L dissapears: (L,T1,M) -> (T1,M)
    %lY=-0.5.*dot_abk(X1-Mu1,(X1-Mu1)./Sigm1,3);   %                       -> (T1,M)
    lY=-0.5.*dot_abk(X1,X1./Sigm1,3);   %                       -> (T1,M)

    clear Sigm1;%Inserted by abk

    % c,const -> (T,M) and then multiply by old Y
    lcoi=log(2.*pi).*(L./2)+0.5.*sum(log(sigm),1);%  (L,M) -> (1,M) 
    lcoef=repmat(log(c')-lcoi,[T1,1]);% c' (1, M) -> (T1,M)
    
    clear lcoi;

    YMM=lcoef+lY;            % ( T1,M ) one mixture per column
    
    clear lcoef;

    YY=lsum(YMM,2);                 % add mixtures ---> (T1,1)
    
    YM = [YM ; YMM];
    
    Y = [Y ; YY];

else
    % repeating, changing dimensions:
    X1=permute(repmat(x',[1,1,M]),[1,3,2]);% (T,L) -> (T,M,L) one per mixture
    
    %size_X1_lmultigauss_mod = size(X1)

    clear x;%Inserted by abk
    
    Mu1=permute(repmat(mus,[1,1,T]),[3,2,1]);     % (L,M) -> (T,M,L)

    clear mus;%Inserted by abk

    if DEBUG sz_X1=size(X1), sz_Mu1=size(Mu1),pause; end
    %if DEBUG size(X1), size(Mu1); end

    X1 = X1 -Mu1;

    clear Mu1;%Inserted by abk

    Sigm1=permute(repmat(sigm,[1,1,T]),[3,2,1]); % (L,M) -> (T,M,L)

    if DEBUG sz_Sigm1=size(Sigm1),pause; end
    %if DEBUG sz_Sigm1=size(Sigm1), end

    %Y=squeeze(exp( 0.5.*dot(X-Mu,(X-Mu)./Sigm))) % L dissapears: (L,T,M) -> (T,M)
    %lY=-0.5.*dot_abk(X1-Mu1,(X1-Mu1)./Sigm1,3);   %                       -> (T,M)
    lY=-0.5.*dot(X1,X1./Sigm1,3);              %                       -> (T,M)

    clear Sigm1;%Inserted by abk

    if DEBUG size_lY = size(lY),pause; end
    %if DEBUG size_lY = size(lY); end

    % c,const -> (T,M) and then multiply by old Y
    lcoi=log(2.*pi).*(L./2)+0.5.*sum(log(sigm),1); % c,const -> (T,M)
    lcoef=repmat(log(c')-lcoi,[T,1]);

    if DEBUG sz_lcoi=size(lcoi),sz_c=size(c'), sz_lcoef = size(lcoef), sz_lY=size(lY), pause; end
    %if DEBUG size_lcoef = size(lcoef), end

    clear c sigm;%Inserted by abk

    if DEBUG1 lcoi, lcoef, lY, pause;end
    %if DEBUG1 lcoi,lcoef,lY;end

    clear lcoi;%Inserted by abk

    YM=lcoef+lY;            % ( T,M ) one mixture per column

    clear lcoef lY;%Inserted by abk

    Y=lsum(YM,2);                 % add mixtures ---> (T,1) 

    if DEBUG [ size(YM) NaN size(Y) ]; end

end
  