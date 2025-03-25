function [M, range ,n] = A5_SVD(Matrx1, LFM_Duration, Period_1)

[U S V]=svd(Matrx1);

 S(1,1)=0;
%   S(2,2)=0;
M=U*S*V;
[m,n]=size(M);
Tbin=LFM_Duration/Period_1;
T0 = 0; % ns
c=300000000;
range = c*(Tbin*(0:n-1) - T0)/2;  % Range Bins in meters
M=transpose(M);
end

