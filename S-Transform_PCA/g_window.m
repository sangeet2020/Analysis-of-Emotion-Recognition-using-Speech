%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gauss=g_window(length,freq,factor) 
 
% Function to compute the Gaussion window for  
% function Stransform. g_window is used by function 
% Stransform 
% 
%-----Inputs Needed-------------------------- 
% 
%   length-the length of the Gaussian window 
% 
%   freq-the frequency at which to evaluate 
%         the window. 
%   factor- the window-width factor 
% 
%-----Outputs Returned-------------------------- 
% 
%   gauss-The Gaussian window 
vector(1,:)=[0:length-1]; 
vector(2,:)=[-length:-1]; 
vector=vector.^2;     
vector=vector*(-factor*2*pi^2/freq^2); 
% Compute the Gaussion window 
gauss=sum(exp(vector)); 

