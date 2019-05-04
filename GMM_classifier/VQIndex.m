function [I, dst]=VQIndex(X,T,CB,M) 
% Distance function
% Returns the closest index of vectors in X to codewords in CB
% In other words:
% I is a vector. The length of I is equal to the number of columns in X.
% Each element of I is the index of closest codeword (column) of CB to
% coresponding column of X

%M=size(CB,2);
%T=size(X,2);
%N->T  L->M
LNThreshold=64*10000;

if M*T<LNThreshold
    D=zeros(M,T);
    for i=1:M
        D(i,:)=sum((repmat(CB(:,i),1,T)-X).^2,1);
    end
    [dst I]=min(D);
else
    I=zeros(1,T);
    dst=I;
    for i=1:T
        %{
        M = M
        sz_CB = size(CB)
        %}
        D=sum((repmat(X(:,i),1,M)-CB).^2,1);
        [dst(i) I(i)]=min(D);
    end
end

%{
%Original Program
function [I, dst]=VQIndex(X,N,CB,L) 
% Distance function
% Returns the closest index of vectors in X to codewords in CB
% In other words:
% I is a vector. The length of I is equal to the number of columns in X.
% Each element of I is the index of closest codeword (column) of CB to
% coresponding column of X

%L=size(CB,2);
%N=size(X,2);
LNThreshold=64*10000;

if L*N<LNThreshold
    D=zeros(L,N);
    for i=1:L
        D(i,:)=sum((repmat(CB(:,i),1,N)-X).^2,1);
    end
    [dst I]=min(D);
else
    I=zeros(1,N);
    dst=I;
    for i=1:N
        %{
        L = L
        sz_CB = size(CB)
        %}
        D=sum((repmat(X(:,i),1,L)-CB).^2,1);
        [dst(i) I(i)]=min(D);
    end
end
%}