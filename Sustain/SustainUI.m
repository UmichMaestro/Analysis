function varargout = SustainUI(varargin)
% SUSTAINUI MATLAB code for SustainUI.fig
%      SUSTAINUI, by itself, creates a new SUSTAINUI or raises the existing
%      singleton*.
%
%      H = SUSTAINUI returns the handle to a new SUSTAINUI or the handle to
%      the existing singleton*.
%
%      SUSTAINUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUSTAINUI.M with the given input arguments.
%
%      SUSTAINUI('Property','Value',...) creates a new SUSTAINUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SustainUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SustainUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SustainUI

% Last Modified by GUIDE v2.5 03-Nov-2017 22:05:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SustainUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SustainUI_OutputFcn, ...
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


% --- Executes just before SustainUI is made visible.
function SustainUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SustainUI (see VARARGIN)

% Choose default command line output for SustainUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SustainUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SustainUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fName, pName] = uigetfile('*.msm*');
handles.path = [pName fName];
[f,A,start,finish] = OpenBinary(handles.path);
handles.frq = f;
handles.A = A;
guidata(hObject, handles)
set(handles.start, 'String', start);
set(handles.finish, 'String', finish);
set(handles.loops, 'String', 3);
plot(handles.axes1, A');
set(handles.axes1,'XTick',linspace(0,500,51));
hideMean(handles);
if start ~= 0 && finish ~= 0 
    disp('show lines')
    show_mean_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% function WriteSustainPoints( path, start, finish )
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));

newPath = [handles.path '2']; % msm2 => segmented
copyfile(handles.path, newPath); %make copy
WriteSustainPoints(newPath, start, finish);


% --- Executes on button press in playOriginal.
function playOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to playOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Synthesis2( frq, A, fs, window )
fs = 44100;
sig = Synthesis2(handles.frq, handles.A, fs);
sig = sig/32;
player = audioplayer(sig, fs);
%guidata(hObject, handles)
playblocking(player);

% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%function [ out ] = Sustain( frq, A, fs, start, finish, loop )
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
duration = str2double(get(handles.loops,'String'));
random = true; % get(handles.randomCheckbox, 'Value');
fs = 44100;
sig = Sustain(handles.frq, handles.A, fs, start, finish, duration, random);
sig = sig/32;
player = audioplayer(sig, fs);
%guidata(hObject, handles)
playblocking(player);




% --- Executes on button press in segment.
function segment_Callback(hObject, eventdata, handles)
% hObject    handle to segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% delete(handles.rect)
[x,] = ginput(2);
start = round(x(1));
finish = round(x(2));
set(handles.start, 'String', start);
set(handles.finish, 'String', finish);
drawMean(start, finish, handles);

% --- Executes on button press in show_mean.
function show_mean_Callback(hObject, eventdata, handles)
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
drawMean(start, finish, handles)


% --- Executes on button press in hide_mean.
function hide_mean_Callback(hObject, eventdata, handles)
hideMean(handles)


function drawMean(start, finish, handles)
hideMean(handles)
if start == 0 && finish == 0
    return;
end

residuals = zeros(size(handles.A, 1), finish-start+1);

for i = 1:5
    data = handles.A(i,start:finish); % data of i-th harmonics in sustain portion
    
    % draw mean
    meanValue = mean(data);
    rectangle('Position', [start meanValue finish-start 0],'LineWidth',2.0, 'LineStyle',':');
    
    % draw line
    c = polyfit([start:finish], data, 1);
    x = linspace(start, finish, 2);
    y = c(1)*x + c(2);
    line(x,y,'Color','red','LineWidth',2.0);
    
    residuals(i,:) = detrend(data, 'linear');
end

function hideMean(handles)
overlays = findobj(handles.axes1,'LineWidth',2.0);
if isempty(overlays)==false
   delete(overlays) 
end



% --- Executes on button press in runLPC.
function runLPC_Callback(hObject, eventdata, handles)
% hObject    handle to runLPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
order = str2double(get(handles.order,'String'));

% coefficients experiment
% coeff = zeros(16, 16);
% gs = zeros(16, 1);
% for i=1:16
%     model = GenerateModel(handles.A, start, finish, i);
%     gs(i) = model(1, 2);
%     coeff(i,1:i) = model(1, 3:end);
% end

model = GenerateModel(handles.A, start, finish, order);
handles.model = model;
guidata(hObject, handles);

figure(2)
data = handles.A(1,:);
plot(data')
grid on

sustain = data(start:finish);
duration = finish-start+1;
% trend line
c = polyfit([start:finish], sustain, 1);
x = linspace(start, finish, 2);
y = c(1)*x + c(2);
line(x,y,'Color','red');
% residual
r = detrend(sustain, 'linear');
line(start:finish, r, 'Color','black');

figure(3)
% linear prediction
g = sqrt(model(1,2));
a = [1, model(1,3:end)];
noise = g * randn(1, duration);
est = filter(a, 1, noise);
plot(noise, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1);
line(1:duration, r, 'Color','blue', 'LineWidth', 2); % r is detrended data
line(1:duration, est, 'Color', 'red', 'LineWidth', 2);




% --- Executes on button press in play_with_model.
function play_with_model_Callback(hObject, eventdata, handles)
% hObject    handle to play_with_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

runLPC_Callback(hObject, eventdata, handles)

start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
duration = str2double(get(handles.loops,'String'));
ar = get(handles.randomCheckbox, 'Value');
model = handles.model;
sustain = duration*100;
fs = 44100;
A = EnvelopWithModel(handles.A, start, finish, sustain, model, ar);

figure(7)
plot(A')

fs = 44100;
sig = Synthesis2(handles.frq, A, fs);
sig = sig/32;
player = audioplayer(sig, fs);
%guidata(hObject, handles)
playblocking(player);






%%%%%%%%%%%%%%%%%%%%%%%% text fields %%%%%%%%%%%%%%%%%%%


function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start as text
%        str2double(get(hObject,'String')) returns contents of start as a double


% --- Executes during object creation, after setting all properties.
function start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loops_Callback(hObject, eventdata, handles)
% hObject    handle to loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loops as text
%        str2double(get(hObject,'String')) returns contents of loops as a double


% --- Executes during object creation, after setting all properties.
function loops_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function finish_Callback(hObject, eventdata, handles)
% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of finish as text
%        str2double(get(hObject,'String')) returns contents of finish as a double


% --- Executes during object creation, after setting all properties.
function finish_CreateFcn(hObject, eventdata, handles)
% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a double


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in randomCheckbox.
function randomCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to randomCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of randomCheckbox


% --- Executes on button press in exportButton.
function exportButton_Callback(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[folder,name] = fileparts(handles.path);
outPath = sprintf('%s/%s',folder,name);

fs = 44100;
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
duration = str2double(get(handles.loops,'String'));
sustain = duration*100;

disp('no sustain. pure msm.');
sig = Synthesis2(handles.frq, handles.A, fs);
sig = sig/max(sig);
audiowrite(sprintf('%s-no_sustain.wav', outPath), sig, fs);

disp('loop sustain');
sig = Sustain(handles.frq, handles.A, fs, start, finish, duration, false);
sig = sig/max(sig);
audiowrite(sprintf('%s-loop.wav', outPath), sig, fs);

disp('mean sustain');
model = GenerateModel(handles.A, start, finish, 1);
A = EnvelopWithModel(handles.A, start, finish, sustain, model, false);
sig = Synthesis2(handles.frq, A, fs);
sig = sig/max(sig);
audiowrite(sprintf('%s-mean.wav', outPath), sig, fs);

disp('AR with 1 coeff');
A = EnvelopWithModel(handles.A, start, finish, sustain, model, true);
sig = Synthesis2(handles.frq, A, fs);
sig = sig/max(sig);
audiowrite(sprintf('%s-AR1.wav', outPath), sig, fs);

for i = [4,8,12,16]
    disp(sprintf('AR with %d coeff', i));
    model = GenerateModel(handles.A, start, finish, i);
    A = EnvelopWithModel(handles.A, start, finish, sustain, model, true);
    sig = Synthesis2(handles.frq, A, fs);
    sig = sig/max(sig);
    audiowrite(sprintf('%s-AR%d.wav', outPath, i), sig, fs);
end
disp('DONE');
