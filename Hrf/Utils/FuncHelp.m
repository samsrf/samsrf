function varargout = FuncHelp(varargin)
% FUNCHELP MATLAB code for FuncHelp.fig
%      FUNCHELP, by itself, creates a new FUNCHELP or raises the existing
%      singleton*.
%
%      H = FUNCHELP returns the handle to a new FUNCHELP or the handle to
%      the existing singleton*.
%
%      FUNCHELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FUNCHELP.M with the given input arguments.
%
%      FUNCHELP('Property','Value',...) creates a new FUNCHELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FuncHelp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FuncHelp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FuncHelp

% Last Modified by GUIDE v2.5 09-Aug-2023 09:45:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FuncHelp_OpeningFcn, ...
                   'gui_OutputFcn',  @FuncHelp_OutputFcn, ...
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


% --- Executes just before FuncHelp is made visible.
function FuncHelp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FuncHelp (see VARARGIN)

% Choose default command line output for FuncHelp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FuncHelp wait for user response (see UIRESUME)
% uiwait(handles.figure1);

Pn = [fileparts(which('FuncHelp.m')) '/../'];
Ds = dir([Pn '*.']);
Ds = {Ds.name}';

% Loop thru functions
L = {}; f = 1;
for i = 1:length(Ds)
    % Only include actual folders
    if Ds{i}(1) ~= '.'
        % Add folder name
        L{f} = Ds{i};
        f = f + 1;
        % Find functions
        Fs = dir([Pn Ds{i} '/*.m']);
        Fs = {Fs.name}';
        % Loop thru functions
        for j = 1:length(Fs)
            L{f} = ['   ' Fs{j}(1:end-2)];
            f = f + 1;
        end
    end
end
set(handles.listbox1, 'Value', 1);
set(handles.listbox1, 'String', L);
set(handles.text2, 'String', '');


% --- Outputs from this function are returned to the command line.
function varargout = FuncHelp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

contents = cellstr(get(hObject,'String')); 
f = get(hObject,'Value');
s = contents{f};
if s(1) == ' '
    % Display help text
    set(handles.text2, 'String', help(s(4:end)));
else
    % Folder name has no text
    set(handles.text2, 'String', '');
end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
