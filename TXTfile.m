% hh% 
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
A=Array256(:,2:18);
fs=60;
FFT_A=abs(fft(A));
k=0:fs/length(A(:,1)):fs-(fs/length(A(:,1)));
figure;plot(k,FFT_A(:,6));xlim([0 2]);
figure;plot(k,FFT_A(:,7));xlim([0 2]);
figure;plot(k,FFT_A(:,8));xlim([0 2]);
figure;plot(k,FFT_A);figure;plot(k,A);