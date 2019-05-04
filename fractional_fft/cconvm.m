function M=cconvm(v);
%Generates circular Convm matrix

v=v(:);N=length(v);
M=zeros(N,N);dum=v;
for k=1:N,
	M(:,k)=dum;
	dum=[dum(N); dum(1:N-1)];
end;

%%%%%%