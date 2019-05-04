function [red,green,blue] = spectro_func(y_initial,fs)

envelope = imdilate(abs(y_initial), true(1501, 1));
% Find the quiet parts.
quietParts = envelope < 0.06; % Or whatever value you want.
% Cut out quiet parts and plot.
yEdited = y_initial; % Initialize
yEdited(quietParts) = [];
%%
x = yEdited(:); 
xlen = length(x);                   % length of the signal
wlen = 256;      %1024 for savee and 256 for EmoDB                  % window length (recomended to be power of 2)
hop = wlen/4;                       % hop size (recomended to be power of 2)
nfft = 4096;                        % number of fft points (recomended to be power of 2)

% perform STFT
[S, f, t] = stft(x, wlen, hop, nfft, fs); % S = STFT matrix (only unique points, time across columns, freq across rows)
S_log = 20*log10(abs(S));



% Scaling intensity to [0, 1]
normA = S_log - min(S_log(:));
R_spec = normA./ max(normA(:));

%% Code for Spectrogram
h=imagesc(t,f,S_log);%[112.5 49.5 680 537]
colormap('jet');
shading interp;
axis xy;
axis off;

saveas(h,'spectro.png');

img1 = imread('spectro.png'); % Read image
img = imcrop(img1,[155.5 67.5 931 732]); %[116.5 50.5 677 533]
red = img(:,:,1); % Red channel
green = img(:,:,2); % Green channel
blue = img(:,:,3); % Blue channel

red=mean(mean(red));
green=mean(mean(green));
blue=mean(mean(blue));