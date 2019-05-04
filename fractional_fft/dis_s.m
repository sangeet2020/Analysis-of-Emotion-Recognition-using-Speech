function [Evec,Eval]=dis_s(N,ord)
%function [Evec,Eval]=dis_s(N)
%Returns sorted eigenvectors and eigenvalues of corresponding vectors

if nargin==1, ord=2;end;

%%Construct S Matrix
%S=diag(2*cos(2*pi/N*([0:N-1])))+diag(ones(1,N-1),1)+diag(ones(1,N-1),-1);
%S(1,N)=1;S(N,1)=1;
%%
S=creates(N,ord);

%%%%%%

%Construct P matrix

p=N;
r=floor(p/2);
P=zeros(p,p);

P(1,1)=1;
even=~rem(p,2);
for k=1:r-even,
	P(k+1,k+1)=1/sqrt(2);
	P(k+1,p-(k+1)+2)=1/sqrt(2);
end;

if (even), P(r+1,r+1)=1; end;

for k=r+1:p-1,
	P(k+1,k+1)=-1/sqrt(2);
	P(k+1,p-(k+1)+2)=1/sqrt(2);
end;

%%%%%%

CS=P*S*P';C2=CS(floor(1:N/2+1),floor(1:N/2+1));S2=CS(floor(N/2+2):N,floor(N/2+2):N);

[vc,ec]=eig(C2);[vs,es]=eig(S2);
qvc=[vc ;zeros(ceil(N/2-1),floor(N/2+1))];
SC2=P*qvc;	%Even Eigenvector of S

qvs=[zeros(floor(N/2+1),ceil(N/2-1));vs];
SS2=P*qvs;	%Odd Eigenvector of S

es=diag(es);ec=diag(ec);
[d,in]=sort(-ec);
SC2=SC2(:,in);
ec=ec(in);

[d,in]=sort(-es);
SS2=SS2(:,in);
es=es(in);

if rem(N,2)==0,
	S2C2=zeros(N,N+1);
	SS2(:,size(SS2,2)+1)=zeros(N,1);
	S2C2(:,[0:2:N]+1)=SC2;
	S2C2(:,[1:2:N]+1)=SS2;
	S2C2(:,N)=[];
else 
	S2C2=zeros(N,N);
	S2C2(:,[0:2:N]+1)=SC2;
	S2C2(:,[1:2:N-1]+1)=SS2;
end;

Evec=S2C2;
Eval=(-j).^[ 0:N-2 (N-1)+even];