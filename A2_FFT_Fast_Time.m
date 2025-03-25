function [fre, M] = A2_FFT_Fast_Time(Matrx1, sampling_Rate, n,IF_BW_1)
fs=sampling_Rate;
M=fft((transpose(Matrx1)));
B1=abs(M);
% fre=0:fs/n:fs/(IF_BW_1*(fs/n))-(fs/n);

fre=0:fs/n:fs-(fs/n);
M=B1;
FastTimeFFTview=B1;
end

