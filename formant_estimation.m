clc;clear all;
% INPUTS
%   x        (vector) of size Nx1 which contains signal
%   fs       (scalar) the sampling frequency
%   [ncoef]  (scalar) the number of coefficients. The default uses
%              ncoef = 2 + fs / 1000;
%             as a rule of thumb. 
% OUTPUTS
%   a        (vector) of size ncoefx1 which contains LPC coefficients
%   P        (scalar) variance (power) of the prediction error
%   e        (vector) of size Nx1 which contains residual error signals

% function [a P e] = spLpc(x, fs, ncoef)
[x fs] = audioread('00.wav');

%  if ~exist('ncoef', 'var') || isempty(ncoef)
ncoef = 2 + round(fs / 1000); % rule of thumb for human speech
%  end
 [a P] = lpc(x, ncoef);
%  if nargout > 2,
%     est_x = filter([0 -a(2:end)],1,x);    % Estimated signal
%     e = x - est_x;                        % Residual signal
%  end 
% end
% function [F] = spFormantsLpc(a, fs)
 r = roots(a);
 r = r(imag(r)>0.01);
 F = sort(atan2(imag(r),real(r))*fs/(2*pi));
% end