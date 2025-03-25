
function [n range] =A1_Scan_Matrix_Plot(Matrx1,LFM_Duration,Period_1,j,IF_BW_1,sampling_Rate)
%% Extracting Range

 [m,n]=size(Matrx1);
assignin('base', 'Matrx1', Matrx1);
% Tbin = abs(qqq);  % ns
Tbin=LFM_Duration/Period_1;
T0 = 0;
c=300000000;
% n=n/(IF_BW_1/sampling_Rate);
range = (c*(Tbin*(0:n-1) - T0)/2)/(IF_BW_1/sampling_Rate);  % Range Bins in meters

% y1=fft(T11(:,2));                           
% fre=0:fs/C11:fs-(fs/C11);                   assignin('base', 'fre', fre); %-(fs/L)
% fre1=0:fs/C11:fs;
% % Tbin=LFM_Duration/Period_1;               assignin('base', 'Tbin', Tbin);
% Fm=IF_BW_1/LFM_Duration; 
% % Fm=sampling_Rate/LFM_Duration;
% T0 = 0;
% c=300000000;
% nn=length(T11);    
% assignin('base', 'nn', nn);
% range = c*(fre1.')/(2*Fm);



assignin('base', 'range', range);

end