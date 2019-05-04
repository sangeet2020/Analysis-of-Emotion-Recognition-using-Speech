% function [st,t,f] = st_temp(timeseries)
function st_pca = st_temp(timeseries)
% clc;clear variables; close all;
% 
% [timeseries, fs] = audioread('001.wav');

tic;

TRUE = 1;
FALSE = 0;

%%% DEFAULT PARAMETERS [change these for your particular application]
% verbose = TRUE;
% removeedge= FALSE;
analytic_signal = FALSE;
factor = 1;
%%% END of DEFAULT PARAMETERS

% Change to column vector
if size(timeseries,2) > size(timeseries,1)
    timeseries=timeseries';
end

minfreq = 100;
%maxfreq = fix(length(timeseries)/2); %fix same as floor function
maxfreq=8000;
samplingrate=16000;
freqsamplingrate=300;

% calculate the sampled time and frequency values from the two sampling rates
t = (0:length(timeseries)-1)*samplingrate;
spe_nelements =ceil((maxfreq - minfreq+1)/freqsamplingrate);
f = (minfreq + [0:spe_nelements-1]*freqsamplingrate)/(samplingrate*length(timeseries));
% disp(sprintf('The number of frequency voices is %d',spe_nelements));



% The actual S Transform function is here:
st = strans(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate,analytic_signal,factor);
st_matrix=abs(st);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_pca=pca(st_matrix');


%%%%%%%%%%%

%% end strans function

%%%%%%%%%%%%%%%


