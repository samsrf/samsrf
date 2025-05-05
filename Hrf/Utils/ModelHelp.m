function varargout = ModelHelp(varargin)
% MODELHELP MATLAB code for ModelHelp.fig
%      MODELHELP, by itself, creates a new MODELHELP or raises the existing
%      singleton*.
%
%      H = MODELHELP returns the handle to a new MODELHELP or the handle to
%      the existing singleton*.
%
%      MODELHELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELHELP.M with the given input arguments.
%
%      MODELHELP('Property','Value',...) creates a new MODELHELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModelHelp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModelHelp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Last Modified by GUIDE v2.5 15-Apr-2022 03:38:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModelHelp_OpeningFcn, ...
                   'gui_OutputFcn',  @ModelHelp_OutputFcn, ...
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


%% --- Executes just before ModelHelp is made visible.
function ModelHelp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModelHelp (see VARARGIN)

% Choose default command line output for ModelHelp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ModelHelp wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Retrieve analysis function
contents = cellstr(get(handles.popupmenu1,'String')); 
analfunc = contents{get(handles.popupmenu1,'Value')}; 
analfunc = analfunc(strfind(analfunc, '(')+1:strfind(analfunc, ')')-1);
set(handles.listbox1, 'String', ModelHelpList(analfunc));
set(handles.text2, 'String', ModelHelpText(analfunc, 'Name'));

%% --- Outputs from this function are returned to the command line.
function varargout = ModelHelp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Retrieve analysis function
contents = cellstr(get(hObject,'String')); 
analfunc = contents{get(hObject,'Value')}; 
analfunc = analfunc(strfind(analfunc, '(')+1:strfind(analfunc, ')')-1);
set(handles.listbox1, 'Value', 1);
set(handles.listbox1, 'String', ModelHelpList(analfunc));
set(handles.text2, 'String', ModelHelpText(analfunc, 'Name'));


%% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Retrieve analysis function
contents = cellstr(get(handles.popupmenu1,'String')); 
analfunc = contents{get(handles.popupmenu1,'Value')}; 
analfunc = analfunc(strfind(analfunc, '(')+1:strfind(analfunc, ')')-1);
% Retrieve parameter name
contents = cellstr(get(hObject,'String'));
param = contents{get(hObject,'Value')};
% Set text about parameter
set(handles.text2, 'String', ModelHelpText(analfunc, param));


%% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
