function lpccs = lpcc_coeff(signal,fs)

% clear all;clc;close all;
% [signal fs] = audioread('03a02Wb.wav');
lpccs = msf_lpcc(signal,fs);
formants = formant_estimation(signal,fs);
lpccs(:,17) = formants(1);
lpccs(:,18) = formants(2);
lpccs(:,19) = formants(3);
