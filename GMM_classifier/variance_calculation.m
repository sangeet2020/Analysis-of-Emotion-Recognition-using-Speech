%This function is written for to calculate the variances of the clusters 
%These variances will be used as the initialising variances for G.M.M
%classifier
function [initial_variance] = variance_calculation(X, L, I, M)
%L = size(X, 1);
%I is row vector giving index of codebook vector corresponding to each
%column vector of x; 1 col vector of x is a feature vector
%M is the size of the codebook and also the GMM

%initial_variance = [];
initial_variance = zeros(L, M);

p=ones(1,length(I)); % row vector of 1s of length(I)

for i=1:M
    p1 = i*p;
    p2 = (p1==I);
    p3 = logical(p2);
    %vectors_each_cluster = X(:,p3);
    %variance_each_cluster = var(vectors_each_cluster');
    variance_each_cluster = var(X(:,p3), 0, 2);
    %initial_variance(:, i) = [initial_variance variance_each_cluster'];
    initial_variance(:, i) = variance_each_cluster;  %(LXM)
end
    