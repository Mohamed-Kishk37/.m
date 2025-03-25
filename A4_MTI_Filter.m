function [ j, M, range ,n] = A4_MTI_Filter(Matrx1, LFM_Duration, Period_1, j, k,sampling_Rate,IF_BW_1)
    if     k==1
         for i=1:j-1
             M(:,i)=Matrx1(:,i) - Matrx1(:,i+1) ;
         end    
    elseif k==2
         for i=1:j-2
             M(:,i)=Matrx1(:,i) - 2* Matrx1(:,i+1) + Matrx1(:,i+2);
         end    
    elseif k==3
         for i=1:j-3
             M(:,i)=Matrx1(:,i) - 3* Matrx1(:,i+1) + 3* Matrx1(:,i+2) - Matrx1(:,i+3);
         end 
    end            
[dd,d]=size(M);
% Tbin = abs(qqq);  % ns
Tbin=LFM_Duration/Period_1;
T0 = 0; % ns
c=300000000;
range = (c*(Tbin*(0:dd-1) - T0)/2)/(IF_BW_1/sampling_Rate);  % Range Bins in meters
n=d;
j=j-k;
end

