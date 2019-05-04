function c = vqlbg(d, k) 
% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm 
% Inputs: 
% d contains training data vectors (one per column) (d is a matrix 
%19x139 in our case, the output of the mfcc function) 
% k is number of centroids required 
% 
% Outputs: 
% c contains the result VQ codebook (k columns, one for each centroids) 
% First, we initialize a few values which will be used in the function 
EPSILON = 0.01; 
THRESHOLD = 3; 
% We attribue the value Infinite to the distorsion because at the...
%beginning, we don't know it and it will decrease 
distorsion_prime = Inf; 
% We attribute the value 100 to the convergence (the convergence is the .. 
%expression to be evaluate with Epsilon in order to evaluate the 
%convergence) We put 100, but it's an arbitrary value. 
convergence = 100; 
% we compute the mean of the second dimension of the d-matrix. The result, 
%c, is a column-vector representating the mean of each row. 
c = mean (d,2); 
dimension_c = size(c); 
nb_columns_of_c = dimension_c(2); 
while (nb_columns_of_c < k) 
% We split each centroid 
y_plus = c*(1 + EPSILON); 
y_moins = c*(1 - EPSILON); 
c = [y_plus y_moins]; 
dimension_c = size(c); 
nb_columns_of_c = dimension_c(2); 
while (convergence > THRESHOLD) 
% Cluster vectors
z=disteu(d, c); 
[m, ind] = min(z, [], 2); 
% Find centroids 
for j = 1:nb_columns_of_c 
c(:, j) = mean(d(:, find(ind == j)), 2); 
end 
% Update distorsion 
distorsion = sum (m); 
% Update of the value 'convergence' 
convergence = (distorsion_prime - distorsion)/distorsion; 
% Update of the value distorsion_prime 
distorsion_prime = distorsion; 
end 
end 
