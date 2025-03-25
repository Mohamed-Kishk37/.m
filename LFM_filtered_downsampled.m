clear all;close all;clc;
%%
%% NOTES:

%(1)After inserting a new csv file, check which 
%%%%coulomn represents the echo and the reference signals.

%(2)plot the csv file and extract the start time of the first sawtooth signal,
%%%%then get the corresponding row number from the csv file, then update the
%%%%variable (start) with the new row number.

%(3)check the number of points in the csv file and update the variable
%%%%(number_of _points) with this value.

% T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\measurement 06-07-2017\tek0000CH2.csv',15,0,[15 0 10000 1]);
T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\target_3m.csv',15,0,[15 0 10000 2]);

% T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\MEASUREMENT20-06-2017\Sheet_3m\tek0000ALL.csv',14,0,[14 0 10014 2]);
% T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\mesurements 22-06-2017 2mixers\No target\tek0000ALL.csv',15,0,[15 0 10014 2]);

plot(T1(:,1),T1(:,2)/max(T1(:,2)),'r');
hold on;
% figure;
plot(T1(:,1),T1(:,3)/max(T1(:,3)));
title('PLOT')
figure;
%%
%%
%% converting csv file into a matrix
%%
q=T1(1,1);
qq=T1(2,1);
qqq=q-qq;
sampling_Rate=abs(1/qqq);
% SR=round(samplingRate);
LFM_Duration = 640e-9;
number_of_points=length(T1(:,2)); 
p=number_of_points;
start=1547;%63  58 1160
IF_BW=500e6;
Detection_Range = (3e8)*LFM_Duration/2; 
Range_Resolution=(3e8)/(2*sampling_Rate);
Period1=round(Detection_Range/Range_Resolution); 

Period=round(sampling_Rate*LFM_Duration);
freq_per_pin=IF_BW/Period;
stop=start+Period-1;

for i=1:Period:p-2*Period
T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\target_3m.csv',start+i,0,[start+i 0 stop+i 2]);
% T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\measurement 06-07-2017\tek0000CH2.csv',start+i,0,[start+i 0 stop+i 2]);

% T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\MEASUREMENT20-06-2017\Sheet_1m\tek0000ALL.csv',start+i,0,[start+i 0 stop+i 2]);

% T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\mesurements 22-06-2017 2mixers\No target\tek0000ALL.csv',start+i,0,[start+i 0 stop+i 2]);


if i==1
    j=1;
else
    j=round(i/Period)+1;
end    
 Matrx(j,:) = T(:,3); %#ok<*SAGROW>
 Reference_Matrix(j,:) = T(:,2);
%  zz=xcorr(Matrx(j,:),Reference_Matrix(j,:));
%  plot(zz,'k');
%  hold on;
 
 
 
 
 
plot(T(:,1),T(:,2),'r');
title(['scan no: ',num2str(j)])
figure;
% plot(T(:,1),T(:,3),'r');
% title(['scan no: ',num2str(j)])
% figure;
end
%%
%%
%%
% 
% z11=xcorr(Matrx(1,1:Period),Reference_Matrix(1,1:Period));
% z22=xcorr(Matrx(1,1:Period),Reference_Matrix(2,1:Period));
% z33=xcorr(Matrx(1,1:Period),Reference_Matrix(3,1:Period));
% z44=xcorr(Matrx(1,1:Period),Reference_Matrix(4,1:Period));
% z55=xcorr(Matrx(1,1:Period),Reference_Matrix(5,1:Period));
% % uncomment if j>5
% % z66=xcorr(Matrx(j,start:stop),Reference_Matrix(j,start+5*period:stop+5*period));
% plot(z11,'r');
% hold on;
% plot(z22,'g');
% hold on;
% plot(z33,'k');
% hold on;
% plot(z44,'b');
% hold on;
% plot(z55,'y');
% hold on; 
% figure;

%without fan
T11=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\without fan\tek0001ALL.csv',15,0,[15 0 10000 2]);
plot(T11(:,1),T11(:,2)/max(T11(:,2)));
title('without fan')
figure;

%so far
T22=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\so far\tek0003ALL.csv',15,0,[15 0 10000 2]);
plot(T22(:,1),T22(:,3)/max(T22(:,3)));
title('so far')
figure;

%no target
T33=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\no_target.csv',15,0,[15 0 10000 2]);
plot(T33(:,1),T33(:,3)/max(T33(:,3)));
title('no target')
figure;

%target 3m
T44=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\target_3m.csv',15,0,[15 0 10000 2]);
plot(T44(:,1),T44(:,3)/max(T44(:,3)));
title('target 3m')
figure;

%target 6m
T55=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\target_6m.csv',15,0,[15 0 10000 2]);
plot(T55(:,1),T55(:,3)/max(T55(:,3)));
title('target 6m')
figure;


z99=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T1(1547:9986,3)/max(T1(1547:9986,3)));
z11=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T11(1397:9986,3)/max(T11(1397:9986,3)));
z22=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T22(285:9986,3)/max(T22(285:9986,3)));
z33=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T33(1285:9986,3)/max(T33(1285:9986,3)));
z44=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T44(3129:9986,3)/max(T44(3129:9986,3)));
z55=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T55(1448:9986,3)/max(T55(1448:9986,3)));
% uncomment if j>5
% z66=xcorr(Matrx(j,start:stop),Reference_Matrix(j,start+5*period:stop+5*period));
plot(z11,'r');
hold on;
plot(z22,'m');
hold on;
plot(z33,'y');
hold on;
plot(z44,'g');
hold on;
plot(z55,'k');
hold on; 
plot(z99,'b');
hold on; 
figure;


z99=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T1(1547:3147,3)/max(T1(1547:3147,3)));
z11=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T11(1397:2997,3)/max(T11(1397:2997,3)));
z22=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T22(285:1885,3)/max(T22(285:1885,3)));
z33=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T33(1285:2885,3)/max(T33(1285:2885,3)));
z44=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T44(3129:4729,3)/max(T44(3129:4729,3)));
z55=xcorr(T1(1547:3147,2)/max(T1(1547:3147,2)),T55(1448:3048,3)/max(T55(1448:3048,3)));
% uncomment if j>5
% z66=xcorr(Matrx(j,start:stop),Reference_Matrix(j,start+5*period:stop+5*period));
plot(z11,'r');
hold on;
plot(z22,'m');
hold on;
plot(z33,'y');
hold on;
plot(z44,'g');
hold on;
plot(z55,'k');
hold on; 
plot(z99,'b');
hold on; 
figure;

%%
%%
%% Extracting Range
%%
[m,n]=size(Matrx);
% Tbin = abs(qqq);  % ns
Tbin=650e-9/Period;
T0 = 0; % ns
% c = 0.29979;  % m/ns
c=300000000;
range = c*(Tbin*(0:n-1) - T0)/2;  % Range Bins in meters
pcolor(range,1:j,transpose(Matrx).'), shading interp;
%colormap(hot)
title('CSV file content')
figure;
%% 
%%
%% applying MTI filter
%%
for i=1:j-1
MatrxMTI(i,:)=Matrx(i,:) - Matrx(i+1,:)
end
pcolor(range,1:j-1,transpose(MatrxMTI).'), shading interp;
%colormap(hot)
title('After MTI')
figure;
%%
%%
%% Applying SVD  
%%
[U S V]=svd(Matrx);
 S(1,1)=0;
z1=U*S*V;
pcolor(range,1:j,transpose(z1).'), shading interp;
%colormap(hot)
title('SVD')
%%
%%
% %% Pull out the raw scans (if saved)
% rawscansI = 1:Period;
% rawscansV =reshape(Matrx,[],Period)';%array of no of scans horizontal and range vertical
% [m,n]=size(Matrx);
% % Tbin = abs(qqq);  % ns
% Tbin=650e-9/Period;
% T0 = 0; % ns
% % c = 0.29979;  % m/ns
% c=300000000;
% range = c*(Tbin*(0:n-1) - T0)/2;  % Range Bins in meters
% figure;
% pcolor(range,1:j,rawscansV.'), shading interp;
% colormap(hot)
% xlabel('range (m)');
% ylabel('scan number)');
% title('CSV file content old method')
%%
%%
%% Variance
%%
% M=var(transpose(Matrx));
M=var((Matrx));
figure;
plot(range,M);
xlabel('range (m)');
ylabel('scan number)');
title('Variance')
%%
%% FFT1 
%%
% A=fft(transpose(Matrx));
% B=real(A);
% freq_range = (freq_per_pin*(0:n-1) - 0); %
% freq_Hz = (sampling_Rate*[0:n-1]/n);     % same as freq
% % n_2=ciel(n/2);
% figure;
% pcolor(freq_Hz,1:m,transpose(B)), shading interp;
% colormap(hot)
% title('FFT1')
% figure;
% plot(freq_range,B(:,2))
% title('FFT1')
%%
%%
% %% FFT with MTI
% A=fft(transpose(MatrxMTI));
% B=real(A);
% freq_range = (freq_per_pin*(0:n-1) - 0); %
% freq_Hz = (sampling_Rate*[0:n-1]/n);     % same as freq
% % n_2=ciel(n/2);
% figure;
% pcolor(freq_Hz,1:m-1,transpose(B)), shading interp;
% colormap(hot)
% title('FFT with MTI')
% figure;
% plot(freq_Hz,B(:,2))
% title('FFT with MTI')
%%
%%
%% FFT
%%
fs=sampling_Rate;
% fs=IF_BW
A=fft(transpose(Matrx));
B=abs(A);
fre=0:fs/n:fs-(fs/n);
figure;
plot(fre,B(:,2)), shading interp;
% pcolor(fre,1:(m),transpose(B)), shading interp;
% colormap(hot)
title('FFT')
%%
%%
%% FFT11
%%
fs=sampling_Rate;
% fs=IF_BW
C=fft(transpose(A));
DD=abs(C);
fre=0:fs/n:fs-(fs/n);
figure;
plot(fre,DD(2,:)), shading interp;
%colormap(hot)
title('FFT11')
%%
%%
%% FFT2
%%
fs=sampling_Rate;
% fs=IF_BW
C=fft((Matrx));
D=abs(C);
fre=0:fs/n:fs-(fs/n);
figure;
plot(fre,D(2,:)), shading interp;
%colormap(hot)
title('FFT2')
%%
%%
%% FFT3
%%
fs=sampling_Rate;
% fs=IF_BW
E=fft(T1(:,2));
F=abs(E);
% fre=0:fs/(2*p):(fs/2)-(fs/(2*p)); %% half band
fre=0:fs/(p):(fs)-(fs/(p));
fre1=0:fs/p:fs-2*(fs/p);
figure;
plot(fre,(F(1:p)))
title('FFT3')
%%
%%







% 
% 
% 
% T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\no_target.csv',15,0,[15 0 10000 2]);
% 
% % T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\MEASUREMENT20-06-2017\Sheet_3m\tek0000ALL.csv',14,0,[14 0 10014 2]);
% % T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\mesurements 22-06-2017 2mixers\No target\tek0000ALL.csv',15,0,[15 0 10014 2]);
% 
% plot(T1(:,1),T1(:,2)/max(T1(:,2)),'r');
% hold on;
% % figure;
% plot(T1(:,1),T1(:,3)/max(T1(:,3)));
% title('PLOT')
% figure;
% %%
% %%
% %% converting csv file into a matrix
% %%
% q=T1(1,1);
% qq=T1(2,1);
% qqq=q-qq;
% sampling_Rate=abs(1/qqq);
% % SR=round(samplingRate);
% LFM_Duration = 640e-9;
% number_of_points=length(T1(:,2)); 
% p=number_of_points;
% start=1160;%63  58
% IF_BW=500e6;
% Detection_Range = (3e8)*LFM_Duration/2; 
% Range_Resolution=(3e8)/(2*sampling_Rate);
% Period1=round(Detection_Range/Range_Resolution); 
% 
% Period=round(sampling_Rate*LFM_Duration);
% freq_per_pin=IF_BW/Period;
% stop=start+Period-1;
% 
% for i=1:Period:p-2*Period
% T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\Measurements 09-07-2017\9-7-2017\no_target.csv',start+i,0,[start+i 0 stop+i 2]);
% % T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\measurement 06-07-2017\tek0000CH2.csv',start+i,0,[start+i 0 stop+i 2]);
% 
% % T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\MEASUREMENT20-06-2017\Sheet_1m\tek0000ALL.csv',start+i,0,[start+i 0 stop+i 2]);
% 
% % T=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\mesurements 22-06-2017 2mixers\No target\tek0000ALL.csv',start+i,0,[start+i 0 stop+i 2]);
% 
% 
% if i==1
%     j=1;
% else
%     j=round(i/Period)+1;
% end    
%  Matrx(j,:) = T(:,3); %#ok<*SAGROW>
%  Reference_Matrix(j,:) = T(:,2);
%  zz=xcorr(Matrx(j,:),Reference_Matrix(j,:));
%  plot(zz,'k');
%  hold on;
%  
% % plot(T(:,1),T(:,2),'r');
% % title(['scan no: ',num2str(j)])
% % figure;
% % plot(T(:,1),T(:,3),'r');
% % title(['scan no: ',num2str(j)])
% % figure;
% end














% %% using BPF applied for every range bin
% nn1=8; % filter order
% fs1=200E6;%sampling frequency
% s11=(4E9)/fs1;%normalized pass frequency
% s21=(6E9)/fs1;%normalized stop frequency
% [b11,a11] = butter(nn1,[s11 s21],'bandpass');%first bandpass butter worth filter from 3.1-5.1 GHz
% y=zeros(m,n);
% 
% for i=1:n
%   y(:,i) = filter(b11,a11,j(:,i));
% end
% %%%%%%%%%%%%%%%%%%%%
% q3=var(transpose(y));
% %%%%%%%%%%%%%%%%%%%%%
% %using MTI filter applied for all scans in every range bins
% filter1=zeros(m,n-1000);
% a=[1, -1];
% b=1;
% for i = 1:m
%   filter1(i,:) = filter(a,b,y(i,:));
% end
% %%%%%%%%%%%%%%%%%%%%
% q4=var(transpose(filter1));
% %%%%%%%%%%%%%%%%%%%%
% %using BPF applied for every scan in all range bins (human breathing rate
% %0.1~1 Hz.
% nn=3;
% fs=50;
% s1=(0.2)/fs;%normalized pass frequency
% s2=(2)/fs;%normalized pass frequency0
% [b1,a1] = butter(nn,[s1 s2],'bandpass');%escond bandpass butter worth filter from 0.3-0.8 hz
% y5=zeros(m,n-1000);
% %% with MTI
% for i=1:m
%   y5(i,:) = filter(b1,a1,filter1(i,:));
% end
% %%%%%%%%%%%%%%%%%%%%%
% q5=var(transpose(y5));
% %%%%%%%%%%%%%%%%%%%%%
% figure;
% pcolor(range,1:n-1000,y5.'), shading interp;
% colormap(hot)
% xlabel('range (m)');
% ylabel('scan number');title('after BPF with 2nd MTI');









% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % figure;
% % plot(F(:,1),F(:,2));    %/max(F(:,2))
% % xlabel('Frequency (Hz)');
% % ylabel('Amplitude (Volt)');grid;
% 
% 
% 
% % 
% % 
% % % hold on;
% % % plot(t(:,1),t(:,3)/max(t(:,3)),'r');
% % % xlabel('time (seconds)');
% % % ylabel('amplitude (Volts)');grid;
% % 
% % figure;
% % plot(t(:,2)/max(t(:,2)));
% % xlabel('time (seconds)');
% % ylabel('amplitude (Volts)');grid;
% % 
% % % hold on;
% % % plot(t(:,3)/max(t(:,3)),'r');
% % % xlabel('time (seconds)');
% % % ylabel('amplitude (Volts)');grid;
% % 
% % figure;
% % plot(t(121:323,1),t(121:323,2));
% % xlabel('aa');
% % 
% % 
% fs=200E6;
% x=T(:,2);
% x_ds=downsample(x,16);
% fs=12.5E6;
% y=fft(x_ds);
% k=length(x_ds);
% fre=0:fs/k:fs-(fs/k);
% fre1=0:fs/k:fs-2*(fs/k);
% %y(1)=0;
% figure;
% plot(fre,abs(y))
%  
% z=diff(x_ds);
% % for n=2:k
% %     z=x(n)-x(n-1);
% % end
% z1=fft(z);
% 
% figure;
% plot(fre1,abs(z1))
% % nn1=8; % filter order
% % s11=(0.004E6)/fs;%normalized pass frequency
% % s21=(4E6)/fs;%normalized stop frequency
% % [b11,a11] = butter(nn1,[s11 s21],'bandpass');%first bandpass butter worth filter from 3.1-5.1 GHz
% %   phi = filter(b11,a11,z);
% %   z2=fft(phi);
% % 
% %   figure;
% % plot(fre1,abs(z2))
% % xlabel('z2');
% %   phi1 = filter(b11,a11,x_ds);
% %   z3=fft(phi1);
% % 
% %   figure;
% % plot(fre,abs(z3))
