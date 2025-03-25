%%window function for 
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

Bytes_per_PRT_with_Header=18;
Bytes_per_PRT_no_Header=17;

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
z1=round(z/Bytes_per_PRT_with_Header);
ArraySize=zeros(z1,Bytes_per_PRT_with_Header);
% ArraySize
for i=1:z1
   for j=1:Bytes_per_PRT_with_Header 
    Array256(i,j)=PC_CFAR_Array_Magnitude(((i-1)*Bytes_per_PRT_with_Header)+j);
   end 
end
A=Array256(:,2:Bytes_per_PRT_with_Header);
B=transpose(A);
figure;pcolor(B.'), shading interp;
colormap(hot)
xlabel('range(m) (txt plot)');
ylabel('scan number)');
% filter1=zeros(m,n-1000);

xzz=length(A(:,1));
% xzz=max(A(1,:));
%%
%%Peak Search
for i=1: xzz

[aa1 aa2]=max(A(i,:));

maxArray(i,1) = aa1 ; 
end
%%
%Normalization
maxarray=max(maxArray);
for i=1: xzz
AA(i,:)=A(i,:)/maxArray(i,1);

end

%%

%%Squaring (Power density)
for i=1: xzz
AAA(i,:)=(AA(i,:)).^2;

  
end

for ii=1: xzz
    for iii=1:Bytes_per_PRT_no_Header
    cvc= AAA(ii,iii);
if cvc <.5 
    AAA(ii,iii)=0;
end 
    end
end

BB=transpose(AAA);
figure;pcolor(BB.'), shading interp;
colormap(hot)
xlabel('range(m) (txt plot Power Density)');
ylabel('scan number)');

%%


a=[1, -1];
b=1;
filter1=zeros(Bytes_per_PRT_no_Header,length(B));
for i = 1:Bytes_per_PRT_no_Header
  filter1(i,:) = filter(a,b,B(i,:));
end
nn=3;
fs=1000;%/16;
s1=(0.5)/fs;%normalized pass frequency
s2=(6)/fs;%normalized pass frequency
[b1,a1] = butter(nn,[s1 s2],'bandpass');%escond bandpass butter worth filter from 0.3-0.8 hz
y5=zeros(Bytes_per_PRT_no_Header,length(B));
for i=1:Bytes_per_PRT_no_Header
  y5(i,:) = filter(b1,a1,filter1(i,:));
end
figure;pcolor(y5.'), shading interp;
colormap(hot)
xlabel('range(m)(BPF)');
ylabel('scan number)');
[U S V]=svd(y5);
S(1,1)=0;
B1=U*S*V;
figure;pcolor(B1.'), shading interp;
xlabel('range(m)(SVD)' );
ylabel('scan number)');
m2 = var(transpose(B1));
[p1 p2]=max(m2);Target_range=p2*0.75%target's range bin
 fs=1000;%/16;
 y5=transpose(y5);
FFT_y5=abs(fft(y5));
k=0:fs/length(A(:,1)):fs-(fs/length(A(:,1)));
% figure;plot(k,FFT_y5(:,p2));xlim([0 16]);
figure;plot(k,FFT_y5(:,1));xlim([0 16]);
figure;plot(k,FFT_y5(:,2));xlim([0 16]);
figure;plot(k,FFT_y5(:,3));xlim([0 16]);
figure;plot(k,FFT_y5(:,4));xlim([0 16]);
figure;plot(k,FFT_y5(:,5));xlim([0 16]);
figure;plot(k,FFT_y5(:,6));xlim([0 16]);
figure;plot(k,FFT_y5(:,7));xlim([0 16]);
figure;plot(k,FFT_y5(:,8));xlim([0 16]);
figure;plot(k,FFT_y5(:,9));xlim([0 16]);
figure;plot(k,FFT_y5(:,10));xlim([0 16]);
figure;plot(k,FFT_y5(:,11));xlim([0 16]);
xlabel('fft of the target' );
% figure;plot(k,FFT_A(:,7));xlim([0 2]);
% figure;plot(k,FFT_A(:,8));xlim([0 2]);
% figure;plot(k,FFT_A);figure;plot(k,A);