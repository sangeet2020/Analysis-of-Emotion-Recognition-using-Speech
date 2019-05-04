% function E=feature_energy(window)
clc;clear all;
y=audioread('00.wav');
for i=1:length(y)
    e(i)=(y(i))^2;
end
subplot(2,1,1);
plot(y);
grid on;
subplot(2,1,2);
plot(e)
grid on;