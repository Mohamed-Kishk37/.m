% clear all; close all; clc;

% T11=csvread('csv1.csv',20,3,[20 3 219 4]);    %One PRT
T11=csvread('D:\1XlinxRADAR19102016\Masters\DrMabrouk\#Measurements 2 17-01-2019 with matlab code\RECORDS TWR TEK DPO 15012019\TX HORN RX HORN\Bre 1.5m Lab 1\New folder\csv1.csv',1,3,[1 3 249799 4]);   %249620  999799 is the number of points per file 
% plot(T11(:,1),T11(:,2)/max(T11(:,2)),'r');title('CSV file content')
xlabel('Time');ylabel('Normalized Amplitude');
q=T11(1,1);
qq=T11(2,1);
qqq=q-qq;

sampling_Rate=abs(1/qqq);                
fs=sampling_Rate;
LFM_Duration=0.001;         %SET LFM duration(PRT)
IF_BW_1= 200e6;             %SET LFM Bandwidth
number_of_points=249800;    %SET Number of points per file
start_1=0;                  %SET start point from file


offset=1; %Start from this file
step=1;   %Jump (step) files  
max=40;    %Maximum number of files to read  
 
Skip1=1;   %Skipped number of PRT in the same file, to speedup processing 


Detection_Range = (3e8)*LFM_Duration/2;
PRF_1=1/LFM_Duration;
Range_Resolution=(3e8)/(2*IF_BW_1);
% Period_1=length(T11(:,2));
Period_1=round(1*sampling_Rate*LFM_Duration);
p=number_of_points;
stop_1=start_1+Period_1-1;
freq_per_pin=sampling_Rate/Period_1;
Periods_per_frame=number_of_points/Period_1;
% freq_per_pin=IF_BW_1/Period_1;
% Range_Resolution=(3e8)/(2*sampling_Rate);

%%
C11=length(T11(:,2))-1;
y1=abs(fft(T11(:,2)));
y2=y1(3748:length(y1));% 2.25m cable compensation 
C22=length(y2);

fre=0:fs/C11:fs-(fs/C11);                  
fre1=0:fs/C11:fs;
Fm=IF_BW_1/LFM_Duration; 
T0 = 0;
c=300000000;    
range = c*(fre1.')/(2*Fm)-2.25;                  
%range = c*(Tbin*(0:nn-1) - T0)/2;           
range1 =Range_Resolution*(1:Period_1-1);
Rlength=length(range);
% range2=range1;
range2=range(3748 : Rlength);    % 2.25m cable compensation                      

figure;
plot(range2,transpose(abs(y2)));title('1');
xlabel('range');ylabel('(Amplitude)');
%%
% figure;
% L1=ceil(length(T11)/10);   
% 
% y11=y1(1:L1,1);                             
% range111=range(-1:L1,1);                     
% plot(range111,transpose(abs(y11)));title('2');
% xlabel('No.of Range Cells');ylabel('(No.Of periods)');
%%
figure;
L11=ceil(length(T11)/100);
y111=y1(1:L11,1);                            
range111=range(1:L11,1);    
plot(range111,transpose(abs(y111)));
title('Range error correction');title('3');
xlabel('No.of Range Cells');ylabel('(No.Of periods)');

%%
% y1=((T11(1:L,2))/(max(abs(T11(1:L,2)))));
figure;
% yy1=y1(1:L/2,1);                            assignin('base', 'yy1', yy1);
plot(fre1,transpose(abs(y1)));title('4');
xlabel('No.of Range Cells');ylabel('(No.Of periods)');
%%
FullPath1 = ('D:\1XlinxRADAR19102016\Masters\DrMabrouk\#Measurements 2 17-01-2019 with matlab code\RECORDS TWR TEK DPO 15012019\TX HORN RX HORN\Bre 1.5m Lab 1\New folder\csv1.csv');
T111=csvread(FullPath1,0,3,[0 3 1000 4]); %

%% Converting Row Data to fast and slow time Matrix
%%
%% File config

offset1=offset; %Start from this file
step1=step;   %Jump (step) files  
max1=max;    %Maximum number of files to read  
 
Skip=Skip1;   %Skipped number of PRT in the same file, to speedup processing 

count=1;
x=0;
N = max1/step1;
GetPathName1 = ('D:\1XlinxRADAR19102016\Masters\DrMabrouk\#Measurements 2 17-01-2019 with matlab code\RECORDS TWR TEK DPO 15012019\TX HORN RX HORN\Bre 1.5m Lab 1\New folder\');
StringConcat=strcat(GetPathName1,'csv%d.csv');
NewStringConcat=strrep(StringConcat,'\','/'); % string replace    %Fixes the error due to path name (The path name contains (\) character which gives an error)
% Matrx1= zeros();   
% Period_1=ceil(Period_1/(IF_BW_1/sampling_Rate));
% Matrx1=zeros(980:200);
step_1=round((Period_1)*Skip);
zzz=round(p/(step_1)); 
j=0;
for k=offset1:step1:max1
file = sprintf(NewStringConcat,k);
      x=j;
    for i=1:step_1:(p)%-1*step_1
%         tic;
        T1=csvread(file,start_1+i,3,[start_1+i 3 stop_1+i 4]);
%          toc;     
        if i==1
            j=1+x;
        else

            j=round(i/(step_1))+1+x;  assignin('base', 'j', j);
%             plot(T1(:,1),T1(:,2)/max(T1(:,2)),'r');title('CSV file content');
%             drawnow;            
            V=['--File# ',num2str(k),'/',num2str(max1),'  --PRT # ',num2str(j),'/',num2str(zzz*max1)];
            disp(V);
        end
%         tic;
        Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>   +((k-1)*zzz/max1)
%         toc;

    end
%     plot(T1(:,1),T1(:,2)/max(T1(:,2)),'r');title('CSV file content1');
% drawnow;
end

% plot(T1(:,1),T1(:,2)/max(T1(:,2)),'r');title('CSV file content');

% 
%% Extracting Range
[n, range] =A1_Scan_Matrix_Plot(Matrx1,LFM_Duration,Period_1,j,IF_BW_1,sampling_Rate);
figure;
pcolor(range,0:j-1,transpose(Matrx1).'), shading interp;
xlabel('No.of Range Cells');ylabel('(No.Of periods)');
%for (Maximize) plot functions

fs11=PRF_1/Skip;
fre11=0:fs11/j:(fs11)-fs11/j;
mm=Matrx1(:,20);
mmm=abs(fft(mm));
mmm=mmm/max(mmm);
figure;
plot(fre11,mmm);title('Slow Time FFT');
xlim([0 2]);
xlabel('freq(Hz)');ylabel('(Amplitude)');
%filter coe order=50 pass=5  stop=10 Hz
FilterCoe=[0.000290487516503184,0.000618250743032147,0.00113207802906472,0.00174005242535650,0.00229736517390446,0.00257499176282377,0.00228189797842897,0.00111552186378814,-0.00116323382905009,-0.00464049005484343,-0.00916392947568759,-0.0142835514326711,-0.0192387274856467,-0.0230083451995839,-0.0244276936746564,-0.0223613048714801,-0.0159096625449513,-0.00461159696044747,0.0113950303722794,0.0312941411453700,0.0536174947910755,0.0763839077459492,0.0973375655853440,0.114250451791702,0.125240629528860,0.129051199803209,0.125240629528860,0.114250451791702,0.0973375655853440,0.0763839077459492,0.0536174947910755,0.0312941411453700,0.0113950303722794,-0.00461159696044747,-0.0159096625449513,-0.0223613048714801,-0.0244276936746564,-0.0230083451995839,-0.0192387274856467,-0.0142835514326711,-0.00916392947568759,-0.00464049005484343,-0.00116323382905009,0.00111552186378814,0.00228189797842897,0.00257499176282377,0.00229736517390446,0.00174005242535650,0.00113207802906472,0.000618250743032147,0.000290487516503184];
mmmfilter=conv(mmm,FilterCoe);
mmmfilter1=mmmfilter(51:5050);

Matrx2=Matrx1;
j2=j;
range2=range;

%% FFT1
[fre, M] = A2_FFT_Fast_Time(Matrx1, sampling_Rate, n,IF_BW_1)
figure;
pcolor(fre,0:j-1,(M).'),shading interp;title('FFT(Fast time)');
%for (Maximize) plot functions
B1=M;
j1=j;
fre1=fre;

%% FFT slow time
[fre, M] = A3_FFT_Slow_Time(M, sampling_Rate, n, PRF_1, j, Skip)
figure;
pcolor(range,fre,transpose(M()).'), shading interp; 
% colormap(hot)
title('FFT(Slow time)');
%for (Maximize) plot functions
DD3=M;
j3=j;
range3=range;
fre3=fre;
figure;
plot(T1(:,1),T1(:,2),'r');
 
%% Applying SVD  

% [M, range ,n] = A5_SVD(Matrx1, LFM_Duration, Period_1);
% 
% axes(handles.axes6);
% pcolor(range,0:j-1,(M).'), shading interp;
% %colormap(hot)
% title('SVD')

%% FFT1 after SVD
% [fre, M] = A2_FFT_Fast_Time(Matrx1, sampling_Rate, n)
% 
% axes(handles.axes4);
% pcolor(fre,0:j-1,(M).'),shading interp;title('FFT(Fast time)');
% %for (Maximize) plot functions
% B1=M;
% j1=j;
% fre1=fre;

 %% Breating Filter
% 
% % 
% % %using BPF applied for every scan in all range bins (human breathing rate
% % %0.1~1 Hz.
% % nn=3;
% % fss=250;
% % s1=(0.1)/fss;%normalized pass frequency
% % s2=(2)/fss;%normalized pass frequency0
% % [b1,a1] = butter(nn,[s1 s2],'bandpass');%escond bandpass butter worth filter from 0.3-0.8 hz
% % %y5=zeros(10001);
% % with MTI
% %  z=T1(:,2);
% % for i=1:10001
% % %   y6 = filter(b1,a1,y5(i,:));
% % 
% %   y5 = filter(b1,a1,(z(:,i)));
% % end
% % figure;
% % plot(frquency_Hz,(y5));
% % title('with BPF')
% 
% 
% %% MTI
% %%
% if get(handles.MTIcheck,'value')== 1
%     k = getappdata(handles.figure1 , 'k');
% [ j, M, range ,n] = A4_MTI_Filter(M, LFM_Duration, Period_1, j ,k,sampling_Rate,IF_BW_1);                 
% assignin('base', 'M2', M);
% assignin('base', 'range2', range);
% assignin('base', 'j2', j);
% assignin('base', 'n2', n);
% axes(handles.axes6);
% pcolor(range,0:(j-1),(M).'), shading interp;title('A');
% xlabel('Fast time(Range(m))');ylabel('Slow time(No.Of Scans)');
% 
%     
% % assignin('base', 'B2', B2);
% end
% 
% 
% 
% 
% %% Applying SVD 1 
% %%
% if get(handles.SVDcheck1,'value')== 1
% 
% [M, range ,n] = A5_SVD(Matrx1, LFM_Duration, Period_1);
% assignin('base', 'M3', M);
% assignin('base', 'range3', range);
% assignin('base', 'j3', j);
% assignin('base', 'n3', n);
%     if get(handles.MTIcheck,'value')== 1
%         k = getappdata(handles.figure1 , 'k');
%         M=M(1:21 ,1:150-k); % Resize matrix according to MTI k pulse canceller
%     end
% axes(handles.axes6);
% pcolor(range,0:j-1,(M).'), shading interp;
% colormap(hot);
% title('SVD');
% end
% 
% %% FFT slow time
% %%
% [fre, M] = A3_FFT_Slow_Time(M, sampling_Rate, n, PRF_1, j, Skip)
% 
% axes(handles.axes5);
% pcolor(range,fre,transpose(M).'), shading interp; 
% % colormap(hot)
% title('FFT(Slow time)');
% 
% %for (Maximize) plot functions
% DD3=M;
% j3=j;
% range3=range;
% fre3=fre;
% setappdata(handles.figure1, 'DD3', DD3);
% setappdata(handles.figure1, 'j3', j3);
% setappdata(handles.figure1, 'range3', range);
% setappdata(handles.figure1, 'fre3', fre3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % setappdata(handles.figure1, 'DD1', DD1)
% % setappdata(handles.figure1, 'j', j)
% % setappdata(handles.figure1, 'range', range)
% % setappdata(handles.figure1, 'fre', fre)
% 
% %%
% % --- Executes during object creation, after setting all properties.
% function axes1_CreateFcn(hObject, eventdata, handles)
% function axes2_CreateFcn(hObject, eventdata, handles)
% function axes3_CreateFcn(hObject, eventdata, handles)
% function axes4_CreateFcn(hObject, eventdata, handles)
% function axes5_CreateFcn(hObject, eventdata, handles)
% function Maximize1_Callback(hObject, eventdata, handles)
% T11 = getappdata(handles.figure1 , 'T11');
% figure;
% plot(T11(:,1),T11(:,2)/max(T11(:,2)),'r');title('CSV file content')
% function Maximize2_Callback(hObject, eventdata, handles)
% T1 = getappdata(handles.figure1 , 'T1');
% figure;
% plot(T1(:,1),T1(:,2),'r');
% function Maximize3_Callback(hObject, eventdata, handles)
% Matrx2 = getappdata(handles.figure1 , 'Matrx2');
% j2 = getappdata(handles.figure1 , 'j2');
% range2 = getappdata(handles.figure1 , 'range2');
% figure;
% pcolor(range2,0:j2-1,transpose(Matrx2).'), shading interp;title('Matrix');colormap(hot);
% function Maximize4_Callback(hObject, eventdata, handles)
% fre1 = getappdata(handles.figure1 , 'fre1');
% B1 = getappdata(handles.figure1 , 'B1');
% j1 = getappdata(handles.figure1 , 'j1');
% range = getappdata(handles.figure1 , 'range');
% figure;
% pcolor(fre1,0:j1-1,(B1).'), shading interp;colormap(hot);
% function Maximize5_Callback(hObject, eventdata, handles)
% fre3 = getappdata(handles.figure1 , 'fre3');
% DD3 = getappdata(handles.figure1 , 'DD3');
% j3 = getappdata(handles.figure1 , 'j3');
% range3 = getappdata(handles.figure1 , 'range3');
% figure;
% pcolor(range3,fre3,transpose(DD3).'), shading interp; colormap(hot);
% function Maximize6_Callback(hObject, ~, handles)
% function PeriodsSkipped_Callback(hObject, eventdata, handles)
% function PeriodsSkipped_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function PeriodsPerFrame_Callback(hObject, eventdata, handles)
% function PeriodsPerFrame_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function iCounter_CreateFcn(hObject, eventdata, handles)
% function cancel1_Callback(hObject, eventdata, handles)
% function figure1_CloseRequestFcn(hObject, eventdata, handles)
% % Hint: delete(hObject) closes the figure
% delete(hObject);
% function offset_Callback(hObject, eventdata, handles)
% function offset_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function step_Callback(hObject, eventdata, handles)
% function step_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function Max_Callback(hObject, eventdata, handles)
% function Max_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function axes6_CreateFcn(hObject, eventdata, handles)
% % Hint: place code in OpeningFcn to populate axes6
% 
% % Hints: get(hObject,'String') returns contents of PeriodsPerFrame1 as text
% %        str2double(get(hObject,'String')) returns contents of PeriodsPerFrame1 as a double
% 
% function PathNameBox1_Callback(hObject, eventdata, handles)
% % Hints: get(hObject,'String') returns contents of PathNameBox1 as text
% %        str2double(get(hObject,'String')) returns contents of PathNameBox1 as a double
% function PathNameBox1_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function PathNameBtn1_Callback(hObject, eventdata, handles)
% [filename1 pathname1] = uigetfile({'D:/1XlinxRADAR19102016/Masters/Dr.Mabrouk/#measurements 28-12-2017/*.csv'},'file selector')
% fullpathname1=strcat(pathname1, filename1);
% set(handles.PathNameBox1,'String',pathname1);
% set(handles.PathNameBox2,'String',fullpathname1);
% setappdata(handles.figure1, 'pathname1', pathname1);
% function PathNameBox2_Callback(hObject, eventdata, handles)
% function PathNameBox2_CreateFcn(hObject, eventdata, handles)
% % Hint: edit controls usually have a white background on Windows.
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function PathNameBtn2_Callback(hObject, eventdata, handles)
% [filename2 pathname2] = uigetfile({'D:/1XlinxRADAR19102016/Masters/Dr.Mabrouk/#measurements 28-12-2017/*.csv'},'file selector');
% fullpathname2=strcat(pathname2, filename2);
% set(handles.PathNameBox2,'String',fullpathname2);
% setappdata(handles.figure1, 'pathname2', pathname2);
% 
% function Auto1_Callback(hObject, eventdata, handles)
% function Manual1_Callback(hObject, eventdata, handles)
% 
% 
% %% some important functions
% %%  use it for MTI, NO MTI
% % if get(handles.radioButton1, 'value')
% %   % Do whatever you want to do if radio button 1 is selected.
% % elseif get(handles.radioButton2, 'value')
% %   % Do whatever you want to do if radio button 2 is selected.
% % end
% %% Use it to store data in handles structure   
% %  handles = guidata(hObject);
% % handles.counter = 1;
% % guidata(hObject, handles);
% %%%  EX
% % function plusPushButtonCallback(buttonHandle, eventData, handles)
% % % At THIS point, "handles" is just a copy of guidata for the GUI
% % handles.counter = handles.counter + 1;
% % % At THIS point, "handles" is different to the guidata for the GUI
% % guidata(buttonHandle, handles);
% % % At THIS point, "handles" has been put back in as new guidata for the GUI   
% %% hobject eventdata handles usage
% % The only time I use hObject is if I have two objects sharing a callback,
% % e.g. pushbutton1 and pushbutton2 both call foo_Callback. In this scenario,
% % I would use hObject to determine which one called or if we needed to modify the calling object.
% % handles is a structure holding ALL of the handles for the GUI.
% % hObject, is JUST the handle for the GUI object used.
% % eentdata contains specific event data - like what cell got selected on a table/what key got pressed.
% 
% % --- Executes on button press in MTIcheck.
% function MTIcheck_Callback(hObject, eventdata, handles)
% % Hint: get(hObject,'Value') returns toggle state of MTIcheck
% 
% % --- Executes on selection change in PulseCancellerPopup.
% function PulseCancellerPopup_Callback(hObject, eventdata, handles)
% % Hints: contents = cellstr(get(hObject,'String')) returns PulseCancellerPopup contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from PulseCancellerPopup
%     contents = cellstr(get(hObject,'String'));
%     popChoice=contents(get(hObject,'value'));
%     if     (strcmp(popChoice,'2 pulse canceller'))
%          k=1;
%     elseif (strcmp(popChoice,'3 pulse canceller'))     
%          k=2;
%     elseif (strcmp(popChoice,'4 pulse canceller'))     
%          k=3;
%     else
%          k=1;
%     end   
%   setappdata(handles.figure1, 'k', k);
%   
% % --- Executes during object creation, after setting all properties.
% function PulseCancellerPopup_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% cellstr(set(hObject,'String'));
% 
% % --- Executes on button press in SVDcheck.
% function SVDcheck_Callback(hObject, eventdata, handles)
% % Hint: get(hObject,'Value') returns toggle state of SVDcheck
% 
% % --- Executes on button press in SVDcheck1.
% function SVDcheck1_Callback(hObject, eventdata, handles)
% % Hint: get(hObject,'Value') returns toggle state of SVDcheck1
% 
% 
% function FFTlength1_Callback(hObject, eventdata, handles)
% % Hints: get(hObject,'String') returns contents of FFTlength1 as text
% %        str2double(get(hObject,'String')) returns contents of FFTlength1 as a double
% 
% % --- Executes during object creation, after setting all properties.
% function FFTlength1_CreateFcn(hObject, eventdata, handles)
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% 
% function FFTlength2_Callback(hObject, eventdata, handles)
% function FFTlength2_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to FFTlength2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function FFTlength11_Callback(hObject, eventdata, handles)
% % hObject    handle to FFTlength11 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of FFTlength11 as text
% %        str2double(get(hObject,'String')) returns contents of FFTlength11 as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function FFTlength11_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to FFTlength11 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% function FFTlength22_Callback(hObject, eventdata, handles)
% function FFTlength22_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to FFTlength22 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
