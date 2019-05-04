% function F=dFRT(N,a,ord)
%function F=dFRT(N,a,ord) returns the NxN discrete fractional 
%Fourier transform matrix with fractional order 'a'. 
%The optional argument 'ord' is the order of approximation 
%of the S matrix (default 2). 

%Note: This Matlab file has some subfunctions for generating S_{2k}
%      matrices, eigenvector ordering etc. These functions are not
%      visible from the Matlab workspace.

[N fs]=audioread('00.wav');
a=0.98;
ord=2;
global Evec Eval ordp

if nargin==2, ord=2;end;

if (length(Evec)~=N | ordp~=ord),
	[Evec,Eval]=dis_s(N,ord);
	ordp=ord;
end;

even=~rem(N,2);
F=Evec*diag(exp(-j*pi/2*a*([0:N-2 N-1+even])))*Evec';

%%%%%%






%END