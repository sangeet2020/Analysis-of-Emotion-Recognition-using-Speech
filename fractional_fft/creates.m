function S=creates(N,ord)
%Creates S matrix of approximation order ord
%When ord=1, elementary S matrix is returned

ord=ord/2;
dum=[1 -2 1];s=0;
for k=1:ord,
	s=(-1)^(k-1)*prod(1:(k-1))^2/prod(1:2*k)*2*[0 dum(k+2:2*k+1) zeros(1,N-1-2*k) dum(1:k)]+s;
	dum=conv(dum,[1 -2 1]);
end;
S=cconvm(s)+diag(real(fft(s)));

%%%%%%
