clear all;close all;clc; 
%==========================
%==========================================================================
% DESCRIPTION : GENERATE MULTIPLE (Ref)-SIGNAL USED IN SIMULATION OF
%               PULSE-COMPRESSION TECHNIQUE (USING CAPTURED INVERTED LFM 
%               SAMPLES FROM WFT EVERY 100 ns)
% NOTE        : CAPTURED DATA IS ACTIVE LOW (INVERTED) AND SHOULD BE
%               FLIPPED BEFORE ANALYSIS
%==========================================================================

%% READ CAPTURED I AND Q SAMPLES FROM FILES
%------------------------------------------
prompt1 = 'Enter LR (I&Q) component File (with extension): ';
File_Name_I = input(prompt1,'s');
fidI = fopen(File_Name_I);
PC_CFAR = textscan(fidI,'%c');
fclose(fidI);
% CI=CI_Q(5:8);
% CQ=CI_Q(1:4);

%% GENERATE I'S OFFSET-BINARY, TWO'S COMPLEMENT, AND SIGNED SAMPLES
%------------------------------------------------------------------

% % Convert Cell to Array
PC_CFAR_Array = cell2mat(PC_CFAR);
File_size= length(PC_CFAR_Array);
File_size_Q=File_size/2;
PC_CFAR_Array_dec=zeros(File_size_Q,2);
for i=1:File_size_Q
    PC_CFAR_Array_dec(i,:)=[hex2dec(PC_CFAR_Array(1+((i-1)*2))) hex2dec(PC_CFAR_Array(2+((i-1)*2)))];
end
PC_CFAR_Array_Magnitude=zeros(File_size_Q,1);
for i=1:File_size_Q
    PC_CFAR_Array_Magnitude(i)=[PC_CFAR_Array_dec(i,2) + PC_CFAR_Array_dec(i,1)*16];
end

figure
plot(PC_CFAR_Array_Magnitude)

z=length(PC_CFAR_Array_Magnitude);
z1=round(z/18);
ArraySize=zeros(z1,18);
% ArraySize
for i=1:z1
   for j=1:18 
    Array256(i,j)=PC_CFAR_Array_Magnitude(((i-1)*18)+j);
   end 
end
A=Array256(2500:length( Array256)-2500,2:18);
B=transpose(A);B=diff(B);
figure;pcolor(B.'), shading interp;
colormap(hot)
xlabel('range (m)');
ylabel('scan number)');
A=transpose(B);
fft_A=transpose(B);
for i=1:16
    fft_A(:,i)=abs(fft(A(:,i)));
end
fs=1000;
k=0:fs/length(A(:,1)):fs-(fs/length(A(:,1)));

figure;pcolor(k,1:16,fft_A.'), shading interp;
colormap(hot)
xlabel('range (m)');
ylabel('scan number)');

% filter1=zeros(m,n-1000);
% a=[1, -1];
% b=1;
% filter1=zeros(17,length(B));
% for i = 1:17
%   filter1(i,:) = filter(a,b,B(i,:));
% end
% figure;pcolor(filter1.'), shading interp;
% colormap(hot)
% xlabel('range (m)');
% ylabel('scan number)');
% nn=3;
% fs=1000;
% s1=(0.5)/fs;%normalized pass frequency
% s2=(9)/fs;%normalized pass frequency
% [b1,a1] = butter(nn,[s1 s2],'bandpass');%escond bandpass butter worth filter from 0.3-0.8 hz
% figure;freqz(b1,a1);
% 
% y5=zeros(17,length(B));
% for i=1:17
%   y5(i,:) = filter(b1,a1,A(:,i));
% end
% figure;plot(A(:,7));
% hold on;plot(y5(7,:),'r')
% 
% 
% figure;
% pcolor(y5.'), shading interp;
% colormap(hot)
% xlabel('range (m)');
% ylabel('scan number)');
% m2 = var(transpose(y5));
% m3=skewness(transpose(y5));
% [p1 p2]=max(m2);
% Target_range=p2*0.75%target's range bin
% range=[1:17]*0.75;
% figure;
% plot(range,m2);
% figure;
% plot(range,m3);
% y5=transpose(y5);
% FFT_y5=abs(fft(y5));
% k=0:fs/length(A(:,1)):fs-(fs/length(A(:,1)));
% figure;
% plot(k,FFT_y5(:,p2));
% xlim([0 2]);[f1 f2]=max(FFT_y5(:,p2));
% Target_breathing_rate=k(f2)
% figure;plot(k,FFT_A(:,7));xlim([0 2]);
% figure;plot(k,FFT_A(:,8));xlim([0 2]);
% figure;plot(k,FFT_A);figure;plot(k,A);
% y5=transpose(y5);
% [U S V]=svd(y5);
% S(1,1)=0;
% B1=U*S*V;
% figure;
% pcolor(B1.'), shading interp;
% colormap(hot)
% xlabel('range (m)');
% ylabel('scan number)');
% FFT_B1=abs(fft(B1));
% figure;
% plot(k,FFT_B1(:,p2));
% xlim([0 2]);[f1 f2]=max(FFT_B1(:,p2));
% Target_breathing_rate=k(f2)
% x=abs(mean(m3));
% if x<0.1
%     target='false alarm...there is no target'
% end
% if x>0.1
%     target='Human being...target is here'
% end
