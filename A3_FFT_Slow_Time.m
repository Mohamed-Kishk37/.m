function [fre, M] = A3_FFT_Slow_Time(Matrx, sampling_Rate, n, PRF_1, j, Skip)

fs1=PRF_1/Skip;
M=fft(transpose(Matrx));
DD1=abs(M);
%  fre=-fs1/2:fs1/j:(fs1/2)-fs1/j;
fre=0:fs1/j:(fs1)-fs1/j;
M=DD1;
end

