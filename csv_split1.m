T11=csvread('180131_011203_Ch1.csv',0,3,[0 3 1000 4]); %
q=T11(1,1);
qq=T11(2,1);
qqq=q-qq;
fs=abs(1/qqq);
lfm_duration= 640e-9;
samples_per_period=fs*lfm_duration;
 stop=round((((1000000/samples_per_period)-1)*samples_per_period)+5)
% stop=55996421
i=1;

 start=6
% start=54996486

% stop=999941 %%  =(((1000000/64)-1)*64)+5  ....  This number is multiples of 64(min no. of points) to gerantee to start the next csv file with the start of the new period, change it when changing (LFM BW) or (LFM Duration)    
counter=1
for i=1:1:125
file = sprintf('csv%d.csv',i);
T11=csvread('180131_011203_Ch1.csv',start,3,[start 3 stop 4]); %
dlmwrite(file,T11,'coffset',3,'precision',9,'roffset',0);
type csv1.csv;
x=start;
start= stop+1
stop=2*stop-x+1
counter =counter+1
end
