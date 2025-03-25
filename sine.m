close all;
clear all; %#ok<CLALL>
clc;
t=0.1:0.1:100;
y=sin(t/2);
yy=awgn(y,1);
xx=plot(t,y);
BW=200e6;
c=3e8;
Tc=10^-3;
z=5*cos(pi*BW*(t)/Tc);
saveas(xx,'sine.jpg')
xy=abs(fft(yy));
xy1=zeros(ceil(length(xy)/2))
xy1=xy(1:ceil(length(xy)/2));
figure;
plot(t,xy1);

figure;
plot(t,z);
figure;
tt = 0:1/5e3:2;
yy = chirp(tt,0,1,10);
yy1=yy*5;
plot(tt,yy);