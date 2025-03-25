function varargout = TWR(varargin)
% TWR MATLAB code for TWR.fig
%      TWR, by itself, creates a new TWR or raises the existing
%      singleton*.
%
%      H = TWR returns the handle to a new TWR or the handle to
%      the existing singleton*.
%
%      TWR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWR.M with the given input arguments.
%
%      TWR('Property','Value',...) creates a new TWR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TWR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TWR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TWR

% Last Modified by GUIDE v2.5 05-Jan-2018 19:27:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TWR_OpeningFcn, ...
                   'gui_OutputFcn',  @TWR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TWR is made visible.
function TWR_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for TWR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TWR wait for user response (see UIRESUME)
% uiwait(handles.figure1);

loadState(handles);

% --- Outputs from this function are returned to the command line.
function varargout = TWR_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SaveDefaults.
function SaveDefaults_Callback(hObject, eventdata, handles)

saveState(handles);



% --- Executes on button press in RestoreDefualts.
function RestoreDefualts_Callback(hObject, eventdata, handles)

loadState(handles);



function saveState(handles)

state.LFMdurationPRT=get(handles.LFMdurationPRT ,'string');
state.LFMbw=get(handles.LFMbw ,'string');
state.Start=get(handles.Start ,'string');
state.NoOfPoints=get(handles.NoOfPoints ,'string');
state.PeriodsSkipped=get(handles.PeriodsSkipped ,'string');
state.offset=get(handles.offset ,'string');
state.step=get(handles.step ,'string');
state.Max=get(handles.Max ,'string');

save state.mat state


function loadState(handles)

 load('state.mat','state')
 
 set(handles.LFMdurationPRT ,'string',state.LFMdurationPRT);
 set(handles.LFMbw ,'string',state.LFMbw);
 set(handles.Start ,'string',state.Start);
 set(handles.NoOfPoints ,'string',state.NoOfPoints);
 set(handles.PeriodsSkipped ,'string',state.PeriodsSkipped);
  set(handles.offset ,'string',state.offset);
 set(handles.step ,'string',state.step);
 set(handles.Max ,'string',state.Max);
 
function MaxRange_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of MaxRange as text
%        str2double(get(hObject,'String')) returns contents of MaxRange as a double


% --- Executes during object creation, after setting all properties.
function MaxRange_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NoOfPoints_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of NoOfPoints as text
%        str2double(get(hObject,'String')) returns contents of NoOfPoints as a double


% --- Executes during object creation, after setting all properties.
function NoOfPoints_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SamplingRate_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of SamplingRate as text
%        str2double(get(hObject,'String')) returns contents of SamplingRate as a double


% --- Executes during object creation, after setting all properties.
function SamplingRate_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RangeResolution_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of RangeResolution as text
%        str2double(get(hObject,'String')) returns contents of RangeResolution as a double


% --- Executes during object creation, after setting all properties.
function RangeResolution_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeriodNoOfBins_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of PeriodNoOfBins as text
%        str2double(get(hObject,'String')) returns contents of PeriodNoOfBins as a double


% --- Executes during object creation, after setting all properties.
function PeriodNoOfBins_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqPerBin_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of FreqPerBin as text
%        str2double(get(hObject,'String')) returns contents of FreqPerBin as a double


% --- Executes during object creation, after setting all properties.
function FreqPerBin_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PRF_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of PRF as text
%        str2double(get(hObject,'String')) returns contents of PRF as a double


% --- Executes during object creation, after setting all properties.
function PRF_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Stop_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Stop as text
%        str2double(get(hObject,'String')) returns contents of Stop as a double


% --- Executes during object creation, after setting all properties.
function Stop_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LFMdurationPRT_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of LFMdurationPRT as text
%        str2double(get(hObject,'String')) returns contents of LFMdurationPRT as a double



% --- Executes during object creation, after setting all properties.
function LFMdurationPRT_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LFMbw_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of LFMbw as text
%        str2double(get(hObject,'String')) returns contents of LFMbw as a double


% --- Executes during object creation, after setting all properties.
function LFMbw_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Start_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Start as text
%        str2double(get(hObject,'String')) returns contents of Start as a double


% --- Executes during object creation, after setting all properties.
function Start_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in ReadFile.
function ReadFile_Callback(hObject, eventdata, handles)
% close all;
% clear all;
% T11=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\#Measurements 02-01-2018\600M 100MSPS 20M Point\180103_011317_Ch1.csv',15,3,[15 3 10000 4]); %
% T11=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\#measurements_25_12_2017\LFM_measurements\sheet_new_filter\6m.csv',15,0,[15 0 1900 1]); %
% T11=csvread('E:\#MASTERS\fixed target 5m\180103_205727_Ch1.csv',15,3,[15 3 10000 4]); %
 %  T11=csvread('E:\#MASTERS\movement\180103_214840_Ch1.csv',0,3,[0 3 1000 4]); %

set(handles.iCounter,'String','Reading...','enable','on','BackgroundColor','red')

% set(handles.Refresh1,'enable','off')  
%  set(handles.LFMbw,'enable','off') 
%  set(handles.LFMdurationPRT,'enable','off') 
%  set(handles.Start,'enable','off') 
%  set(handles.NoOfPoints,'enable','off') 
%  set(handles.PeriodsSkipped,'enable','off') 
 
%   r step = 1:steps
 T11=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\#measurements 28-12-2017\CSV SPLIT\csv1.csv',0,3,[0 3 24000 4]); %

  



axes(handles.axes1);
plot(T11(:,1),T11(:,2)/max(T11(:,2)),'r');title('CSV file content')
q=T11(1,1);
qq=T11(2,1);
qqq=q-qq;
% global sampling_Rate;
sampling_Rate=abs(1/qqq);
set(handles.SamplingRate,'String',sampling_Rate)

% global fs;
fs=sampling_Rate;

% global LFM_Duration;
LFM_Duration=str2num(get(handles.LFMdurationPRT,'String'));
Detection_Range = (3e8)*LFM_Duration/2;
set(handles.MaxRange,'String',Detection_Range)


PRF_1=1/LFM_Duration;
set(handles.PRF,'String',PRF_1)


% %%%%%%%%%global Period_1;
Period_1=round(2*sampling_Rate*LFM_Duration);
set(handles.PeriodNoOfBins,'String',Period_1)

% number_of_points=length(T11(:,2));
% set(handles.NoOfPoints,'String',number_of_points)

number_of_points=str2num(get(handles.NoOfPoints,'String'));
 
p=number_of_points;

start_1=str2num(get(handles.Start,'String'));
stop_1=start_1+Period_1-1;
set(handles.Stop,'String',stop_1)

IF_BW_1=str2num(get(handles.LFMbw,'String'));
freq_per_pin=IF_BW_1/Period_1;
set(handles.FreqPerBin,'String',freq_per_pin)

Range_Resolution=(3e8)/(2*IF_BW_1);
set(handles.RangeResolution,'String',Range_Resolution)

% IF_BW_1=str2num(get(handles.LFMbw,'String'));
% freq_per_pin=IF_BW_1/Period_1;
% set(handles.FreqPerBin,'String',freq_per_pin)



Periods_per_frame=number_of_points/Period_1;
set(handles.PeriodsPerFrame,'String',Periods_per_frame)

Skip=str2num(get(handles.PeriodsSkipped,'String'));

assignin('base', 'T11', T11);
assignin('base', 'fs', fs);


fs=sampling_Rate;
E1=fft(T11(1:24001,2));
F1=abs(E1);
F1=F1/max(F1);
% % fre=0:fs/(2*p):(fs/2)-(fs/(2*p)); %% half band
fre=0:fs/24001:(fs)-(fs/24001);
% fre1=0:fs/p:fs-2*(fs/p);
figure;
plot(fre,(F1(1:24001)),'r');
title('FFT4')

%%%%%%%%%%%%%% Converting Row Data to fast and slow time Matrix
% global Matrx1;
% global j;





offset1=str2num(get(handles.offset,'String'));
step1=str2num(get(handles.step,'String'));
max1=str2num(get(handles.Max,'String'));

count=0
x=0;
% i=1;
% j=0;


N = max1/step1;
% Creation of a timer waitbar object
TWB = Timerwaitbar(N);
% Loop
% for n = 1:N   
%     % Simulation of a task 
% %     pause(0.1);      
%     % Timer waitbar update
%     TWB.update();    
%     % Loop break if manual cancelation
%     if TWB.isinterrupted()
%         break
%     end    
% end

% Object deletion


for k=offset1:step1:max1
    
    file = sprintf('D:/1XlinxRADAR19102016/Masters/Dr.Mabrouk/#measurements 28-12-2017/CSV SPLIT1/csv%d.csv',k);
%       TWB.update();    
    set(handles.Refresh1,'UserData',1)
    for i=1:(Period_1)*Skip:(p)-2*(Period_1)*Skip

        T1=csvread(file,start_1+i,3,[start_1+i 3 stop_1+i 4]);
        if i==1
            j=1+x;
        else
            j=round(i/((Period_1)*Skip))+1+x;
        end
        Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>
        
        %          if get(handles.cancel1,'UserData')==1
        %              break;
        %          end
    end
   
    TWB.update();
    % Loop break if manual cancelation
    if TWB.isinterrupted()
        break
    end
    guidata(hObject, handles);
    set(handles.Refresh1,'UserData',1)
    x=j;
    count=count+1
%     if get(handles.cancel1,'UserData')==1
%          break;
%     end   
set(handles.ifiles,'String',num2str(count)) 

end

TWB.delete();



% for i=1:(Period_1)*Skip:(p)-2*Period_1
% T1=csvread('csv2.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     if i==1
%         j=1+x;
%     else
%         j=round(i/((Period_1)*Skip))+1+x;
%     end    
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>   
% end
% x=j;
% for i=1:(Period_1)*Skip:(p)-2*Period_1
% T1=csvread('csv3.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     if i==1
%         j=1+x;
%     else
%         j=round(i/((Period_1)*Skip))+1+x;
%     end    
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>   
% end
% x=j;
% for i=1:(Period_1)*Skip:(p)-2*Period_1
% T1=csvread('csv4.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     if i==1
%         j=1+x;
%     else
%         j=round(i/((Period_1)*Skip))+1+x;
%     end    
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>   
% end
% x=j;


%%%%%%%%%%%%%%%%%%%
% Skip=str2num(get(handles.PeriodsSkipped,'String'));
% 
% for i=1:(Period_1)*Skip:(p/4)-2*Period_1
% T1=csvread('csv.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     if i==1
%         j=1
%     else
%         j=round(i/((Period_1)*Skip))+1
%     end    
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>   
% end
% x=i+(Period_1)*Skip
% 
% for i=x:(Period_1)*Skip:(p/2)-2*Period_1  
%     T1=csvread('E:\#MASTERS\movement\180103_214840_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     j=round(i/(Period_1)/Skip)+1  
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>
% end
% x=i+(Period_1)*Skip
% for i=x:(Period_1)*Skip:(p-(p/4))-2*Period_1  
%     T1=csvread('E:\#MASTERS\movement\180103_214840_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     j=round(i/(Period_1)/Skip)+1  
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>
% end
% x=i+(Period_1)*Skip
% for i=x:(Period_1)*Skip:(p)-2*Period_1  
%     T1=csvread('E:\#MASTERS\movement\180103_214840_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%     j=round(i/(Period_1)/Skip)+1  
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>
% end
% x=i+(Period_1)*Skip





%%%%%%%%%%%%%%%%%
% Skip=str2num(get(handles.PeriodsSkipped,'String'));
% 
% for i=1:(Period_1)*Skip:p-2*Period_1
%    
% %     T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\#Measurements 02-01-2018\600M 100MSPS 20M Point\180103_011317_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
% %     T1=csvread('D:\1XlinxRADAR19102016\Masters\Dr.Mabrouk\#measurements_25_12_2017\LFM_measurements\sheet_new_filter\6m.csv',start_1+i,0,[start_1+i 0 stop_1+i 1]);
% %     T1=csvread('E:\#MASTERS\fixed target 5m\180103_205727_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
%   
% T1=csvread('E:\#MASTERS\movement\180103_214840_Ch1.csv',start_1+i,3,[start_1+i 3 stop_1+i 4]);
% 
%     if i==1
%         j=1;
%     else
%         j=round(i/(Period_1)/Skip)+1
%     end    
%     Matrx1(j,:) = T1(:,2); %#ok<*SAGROW>
%     
%     
% end




 set(handles.iCounter,'String','Done','enable','on','BackgroundColor','green')
 set(handles.Refresh1,'enable','on')
 set(handles.LFMbw,'enable','on') 
 set(handles.LFMdurationPRT,'enable','on') 
 set(handles.Start,'enable','on') 
 set(handles.NoOfPoints,'enable','on') 
 set(handles.PeriodsSkipped,'enable','on')  
  
axes(handles.axes2);
plot(T1(:,1),T1(:,2),'r');
title(['scan no: ',num2str(j),'  (One complete period)'])






% 
setappdata(handles.figure1, 'Matrx1', Matrx1)  % Stores contents of Matrx1 figure1 object using the name identifier, 'Matrx1' for sharing this variable with other functions
setappdata(handles.figure1, 'LFM_Duration', LFM_Duration)
setappdata(handles.figure1, 'Period_1', Period_1)
setappdata(handles.figure1, 'j', j)
setappdata(handles.figure1, 'fs', fs)
setappdata(handles.figure1, 'sampling_Rate', sampling_Rate)
setappdata(handles.figure1, 'T11', T11)
setappdata(handles.figure1, 'T1', T1)




% --- Executes on button press in Refresh1.
function Refresh1_Callback(hObject, eventdata, handles)
% close all  ;

Matrx1 = getappdata(handles.figure1 , 'Matrx1'); %Retrieve the data from figure1 object and display it.
LFM_Duration = getappdata(handles.figure1 , 'LFM_Duration');
Period_1 = getappdata(handles.figure1 , 'Period_1');
j = getappdata(handles.figure1 , 'j');
fs = getappdata(handles.figure1 , 'fs');
sampling_Rate = getappdata(handles.figure1 , 'sampling_Rate');



%%%%%%%%%%%%%% Extracting Range
[m,n]=size(Matrx1);
assignin('base', 'Matrx1', Matrx1);

% Tbin = abs(qqq);  % ns
Tbin=LFM_Duration/Period_1;
T0 = 0; % ns
% c = 0.29979;  % m/ns
c=300000000;
range = c*(Tbin*(0:n-1) - T0)/2;  % Range Bins in meters
assignin('base', 'range', range);

axes(handles.axes3);
pcolor(range,1:j,transpose(Matrx1).'), shading interp;title('Matrix')
xlabel('Fast time(Range(m))')
ylabel('Slow time(No.Of Scans)')


%%%%%%%%MTI FILTER

% for i=1:j-1
% Matrx2(i,:)=Matrx1(i,:) - Matrx1(i+1,:);
% end
% [dd,d]=size(Matrx1);
% % Tbin = abs(qqq);  % ns
% Tbin=LFM_Duration/Period_1;
% T0 = 0; % ns
% % c = 0.29979;  % m/ns
% c=300000000;
% range = c*(Tbin*(0:d-1) - T0)/2;  % Range Bins in meters
% 
% axes(handles.axes6);
% pcolor(range,1:j-1,transpose(Matrx2).'), shading interp;title('Matrix')
% xlabel('Fast time(Range(m))')
% ylabel('Slow time(No.Of Scans)')

%%%%%%%%%%%%%%%%%%%%% FFT1
fs=sampling_Rate;
A1=fft((transpose(Matrx1)));
B1=abs(A1);
B1=B1/max(B1);
fre=0:fs/n:fs-(fs/n);
axes(handles.axes4);
pcolor(range,1:j,(B1).'), shading interp;
title('FFT(Fast time)'); 

setappdata(handles.figure1, 'B1', B1)
%%%%%%%%%%%%%%%%%% FFT2


assignin('base', 'j', j); % Add variable to Matlab Workspace

x=j
for i=1:(j-1)
B2(:,i)=B1(:,i) - B1(:,i+1);
assignin('base', 'B2', B2);
end
[dd,d]=size(B2);
assignin('base', 'd', d);
assignin('base', 'dd', dd);
% Tbin = abs(qqq);  % ns
Tbin=LFM_Duration/Period_1;
T0 = 0; % ns
% c = 0.29979;  % m/ns
c=300000000;
range = c*(Tbin*(0:dd-1) - T0)/2;  % Range Bins in meters
assignin('base', 'range', range);
axes(handles.axes6);
pcolor(range,1:(j-1),(B2).'), shading interp;title('A')
xlabel('Fast time(Range(m))')
ylabel('Slow time(No.Of Scans)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=sampling_Rate;
C1=fft(transpose(B2));
DD1=abs(C1);
DD1=DD1/max(DD1);

fre=0:fs/n:fs-(fs/n); 
axes(handles.axes5);
pcolor(range,1:j-1,transpose(DD1).'), shading interp;hold on; 
% colormap(hot)
title('FFT(Slow time)')

setappdata(handles.figure1, 'DD1', DD1)
%%%%%%%%%%%%%%%%%%%%%%%%%%


setappdata(handles.figure1, 'range', range)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes2

% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes3

% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes4

% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes5

% --- Executes on button press in Maximize1.
function Maximize1_Callback(hObject, eventdata, handles)
T11 = getappdata(handles.figure1 , 'T11');
figure;
plot(T11(:,1),T11(:,2)/max(T11(:,2)),'r');title('CSV file content')

% --- Executes on button press in Maximize2.
function Maximize2_Callback(hObject, eventdata, handles)
T1 = getappdata(handles.figure1 , 'T1');
figure;
plot(T1(:,1),T1(:,2),'r');



% --- Executes on button press in Maximize3.
function Maximize3_Callback(hObject, eventdata, handles)
Matrx1 = getappdata(handles.figure1 , 'Matrx1');
j = getappdata(handles.figure1 , 'j');
range = getappdata(handles.figure1 , 'range');
figure;
pcolor(range,1:j,transpose(Matrx1).'), shading interp;title('Matrix')



% --- Executes on button press in Maximize4.
function Maximize4_Callback(hObject, eventdata, handles)
B1 = getappdata(handles.figure1 , 'B1');
j = getappdata(handles.figure1 , 'j');
range = getappdata(handles.figure1 , 'range');
figure;
pcolor(range,1:j,(B1).'), shading interp;


% --- Executes on button press in Maximize5.
function Maximize5_Callback(hObject, eventdata, handles)
DD1 = getappdata(handles.figure1 , 'DD1');
j = getappdata(handles.figure1 , 'j');
range = getappdata(handles.figure1 , 'range');
figure;
pcolor(range,1:j-1,transpose(DD1).'), shading interp; 


% --- Executes on button press in Maximize6.
function Maximize6_Callback(hObject, eventdata, handles)


function PeriodsSkipped_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of PeriodsPerFrame1 as text
%        str2double(get(hObject,'String')) returns contents of PeriodsPerFrame1 as a double


% --- Executes during object creation, after setting all properties.
function PeriodsSkipped_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PeriodsPerFrame_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of PeriodsSkipped as text
%        str2double(get(hObject,'String')) returns contents of PeriodsSkipped as a double


% --- Executes during object creation, after setting all properties.
function PeriodsPerFrame_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function iCounter_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in cancel1.
function cancel1_Callback(hObject, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);



function offset_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double


% --- Executes during object creation, after setting all properties.
function step_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Max as text
%        str2double(get(hObject,'String')) returns contents of Max as a double


% --- Executes during object creation, after setting all properties.
function Max_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate axes6
