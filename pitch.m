% 
% % % Programm of cepstral analysis for finding the pitch of the pitch of the signal
% % 
% [j, fs]=audioread('00.wav');
% t=[0 : length(j)-1]/fs;
% y=fft(j);   % fs is sampling frequency.
% ms1=fs/1000;    % maximum speech Fx at 1000 Hz.
% ms20=fs/50;      % minimum speech Fx at 50 Hz.
% Y=fft(j.*hamming(length(j)));
% C=fft(log(abs(Y)+eps));
% q=(ms1 : ms20);
% [maxamp_at_pitch , fx]=max(abs(C(ms1 : ms20)));
% frequency_pitch=fs/(ms1+fs-1)





m=1;
[y,Fs]=audioread('00.wav'); % input: speech segment
max_value=max(abs(y));
y=y/max_value;
t=(1/Fs:1/Fs:(length(y)/Fs))*1000;
subplot(2,1,1);
plot(t,y);
%xtitle('A 30 millisecond segment of speech','time in milliseconds');
sum1=0;autocorrelation=0;
   for l=0:(length(y)-1)
    sum1=0;
    for u=1:(length(y)-l)
      s=y(u)*y(u+l);
      sum1=sum1+s;
    end
    autocor(l+1)=sum1;
  end
kk=(1/Fs:1/Fs:(length(autocor)/Fs))*1000;
subplot(2,1,2);
plot(kk,autocor);
%xtitle('Autocorrelation of the 30 millisecond segment of speech','time in milliseconds');
auto=autocor(21:160);
  max1=0;
  for uu=1:140
    if(auto(uu)>max1)
      max1=auto(uu);
      sample_no=uu;
    end 
  end
  pitch_period_To=(20+sample_no)*(1/Fs);
  pitch_freq_Fo(m)=1/pitch_period_To;
  m=m+1;

