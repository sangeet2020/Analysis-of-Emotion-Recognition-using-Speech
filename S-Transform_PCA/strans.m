function st = strans(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate,analytic_signal,factor);  

n=length(timeseries); % n=128
original = timeseries;

% If vector is real, do the analytic signal  
%An analytic signal is a complex-valued function that has no negative frequency components. 
%The real and imaginary parts of an analytic signal are real-valued functions related to each other by the Hilbert transform.

% if analytic_signal %analytic signal = false
% %    disp('Calculating analytic signal (using Hilbert transform)'); 
%    %  This is correct! 
%    ts_spe = fft(real(timeseries)); 
%    h = [1; 2*ones(fix((n-1)/2),1); ones(1-rem(n,2),1); zeros(fix((n-1)/2),1)]; 
%    ts_spe(:) = ts_spe.*h(:); 
%    timeseries = ifft(ts_spe); 
% end   
 
% Compute FFT's 
tic;vector_fft=fft(timeseries);tim_est=toc; 
vector_fft=[vector_fft,vector_fft]; 
tim_est = tim_est*ceil((maxfreq - minfreq+1)/freqsamplingrate)   ; 
% disp(sprintf('Estimated time is %f',tim_est));

% Preallocate the STOutput matrix 
st=zeros(ceil((maxfreq - minfreq+1)/freqsamplingrate),n); 
% Compute the mean 
% Compute S-transform value for 1 ... ceil(n/2+1)-1 frequency points 
% disp('Calculating S transform...'); 
if minfreq == 0 
   st(1,:) = mean(timeseries)*(1&[1:1:n]); 
else 
    st(1,:)=ifft(vector_fft(minfreq+1:minfreq+n).*g_window(n,minfreq,factor)); 
end 
 
%the actual calculation of the ST 
% Start loop to increment the frequency point 
i=0;
for banana=freqsamplingrate:freqsamplingrate:(maxfreq-minfreq) 
   st(banana/freqsamplingrate+1,:)=ifft(vector_fft(minfreq+banana+1:minfreq+banana+n).*g_window(n,minfreq+banana,factor)); 
   i=i+1;
end
% End loop to increment the frequency point 

% disp('Finished Calculation'); 
 
%%% end strans function 
