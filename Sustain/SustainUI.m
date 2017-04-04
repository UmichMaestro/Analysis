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

% Last Modified by GUIDE v2.5 03-Apr-2017 20:15:43

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


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%function [ out ] = Sustain( frq, A, fs, start, finish, loop )
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
loops = str2double(get(handles.loops,'String'));
fs = 44100;
sig = Sustain(handles.frq, handles.A, fs, start, finish, loops);
sig = sig/32;
player = audioplayer(sig, fs);
%guidata(hObject, handles)
playblocking(player);



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


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fName, pName] = uigetfile('*.msm');
handles.path = [pName fName];
[f,A,start,finish] = OpenBinary(handles.path);
handles.frq = f;
handles.A = A;
guidata(hObject, handles)
set(handles.start, 'String', start);
set(handles.finish, 'String', finish);
set(handles.loops, 'String', 2);
plot(handles.axes1, A');
set(handles.axes1,'XTick',linspace(0,500,51));




% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% function WriteSustainPoints( path, start, finish )
start = str2double(get(handles.start,'String'));
finish = str2double(get(handles.finish,'String'));
WriteSustainPoints(handles.path, start, finish);



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
