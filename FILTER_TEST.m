clear all;
close all;
clc;


T=0:0.02:60;
fs=1/0.02;
x=3*sin(2*pi*0.3*T)+sin(2*pi*1*T);
plot(T,x);
k=0:(fs/length(T)):fs-(fs/length(T));
fft_x=abs(fft(x));
figure;
plot(k,fft_x);
figure;

[b,a] = butter(10,0.026,'high');
dataOut = filter(b,a,x);
plot(k,abs(dataOut));
figure;
fft1=abs(fft(dataOut))
plot(k,(fft1));

