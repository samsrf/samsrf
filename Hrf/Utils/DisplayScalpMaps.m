function varargout = DisplayScalpMaps(varargin)
% DISPLAYSCALPMAPS MATLAB code for DisplayScalpMaps.fig
%      DISPLAYSCALPMAPS, by itself, creates a new DISPLAYSCALPMAPS or raises the existing
%      singleton*.
%
%      H = DISPLAYSCALPMAPS returns the handle to a new DISPLAYSCALPMAPS or the handle to
%      the existing singleton*.
%
%      DISPLAYSCALPMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAYSCALPMAPS.M with the given input arguments.
%
%      DISPLAYSCALPMAPS('Property','Value',...) creates a new DISPLAYSCALPMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DisplayScalpMaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DisplayScalpMaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DisplayScalpMaps

% Last Modified by GUIDE v2.5 03-Oct-2023 14:39:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplayScalpMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplayScalpMaps_OutputFcn, ...
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


% --- Executes just before DisplayScalpMaps is made visible.
function DisplayScalpMaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DisplayScalpMaps (see VARARGIN)

% Choose default command line output for DisplayScalpMaps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DisplayScalpMaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global Srf Model dcm_obj ph handles

samsrf_clrscr;

% Load map
[f,p] = uigetfile('eeg_*.mat');
load([p filesep f]);
if ~strcmpi(Srf.Values{1}, 'R^2') && ~strcmpi(Srf.Values{1}, 'nR^2')
    set(handles.edit1, 'String', '-Inf', 'Enable', 'off');
end

% Adjust plots
axes(handles.axes1);
samsrf_sensors(Srf, Srf.Values{get(handles.listbox1,'Value')}, str2double(get(handles.edit1,'String')), nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts));

% Adjust input fields
set(handles.listbox1, 'String', Srf.Values); % List of values
set(handles.slider1, 'Min', 0, 'Max', 1, 'Value', 0, 'SliderStep', [1 1]*double(nanmin(diff(unique(Srf.TimePts)))/range(Srf.TimePts)));

% global handles

% Vertex selection hack 
dcm_obj = datacursormode(gcf); % Data cursor mode object
set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function
dcm_obj.removeAllDataCursors();

% Close request function
crfcn = @closereq;
set(gcf, 'CloseRequestFcn', crfcn);

% Figure for pRF profiles
ph = figure('Units', 'normalized', 'Position', [.0469 .4769 .2917 .3889]);


% --- Outputs from this function are returned to the command line.
function varargout = DisplayScalpMaps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Map selection menu
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

global Srf

% Plot scalp map
axes(handles.axes1);
samsrf_sensors(Srf, Srf.Values{get(handles.listbox1,'Value')}, str2double(get(handles.edit1,'String')), nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts));
% Set clipping range
cr = eval(['[' get(handles.edit2,'String') ']']);
if ischar(cr)
    cr = [];
    set(handles.edit2, 'String', cr);
end
if ~isempty(cr)
    if length(cr) == 1
        cr = [-1 1]*abs(cr);
    elseif length(cr) > 2
        cr = cr(1:2);
    end
    caxis(cr);
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



%% Time point slider
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global Srf

% Update slider
t = nanmin(Srf.TimePts) + get(hObject,'Value')*range(Srf.TimePts);
u = unique(Srf.TimePts);
d = abs(u - t);
t = u(d == min(d));
set(hObject, 'Value', (t-nanmin(Srf.TimePts))/range(Srf.TimePts));
set(handles.text3, 'String', ['Time point:   ' num2str(t)]);

% Plot scalp map
axes(handles.axes1);
samsrf_sensors(Srf, Srf.Values{get(handles.listbox1,'Value')}, str2double(get(handles.edit1,'String')), nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts));
% Set clipping range
cr = eval(['[' get(handles.edit2,'String') ']']);
if ischar(cr)
    cr = [];
    set(handles.edit2, 'String', cr);
end
if ~isempty(cr)
    if length(cr) == 1
        cr = [-1 1]*abs(cr);
    elseif length(cr) > 2
        cr = cr(1:2);
    end
    caxis(cr);
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% R^2 threshold
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

global Srf

% Plot scalp map
axes(handles.axes1);
hold off
samsrf_sensors(Srf, Srf.Values{get(handles.listbox1,'Value')}, str2double(get(handles.edit1,'String')), nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts));
% Set clipping range
cr = eval(['[' get(handles.edit2,'String') ']']);
if ischar(cr)
    cr = [];
    set(handles.edit2, 'String', cr);
end
if ~isempty(cr)
    if length(cr) == 1
        cr = [-1 1]*abs(cr);
    elseif length(cr) > 2
        cr = cr(1:2);
    end
    caxis(cr);
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Clipping range
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

global Srf

% Plot scalp map
axes(handles.axes1);
samsrf_sensors(Srf, Srf.Values{get(handles.listbox1,'Value')}, str2double(get(handles.edit1,'String')), nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts));
% Set clipping range
try
    cr = eval(['[' get(handles.edit2,'String') ']']);
    if ~isempty(cr)
        if length(cr) == 1
            cr = [-1 1]*abs(cr);
            set(handles.edit2, 'String', num2str(cr));
        elseif length(cr) > 2
            cr = cr(1:2);
        end
        caxis(cr);
    end
catch
    cr = [];
    set(handles.edit2, 'String', cr);
end    


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Vertex selection function
function txt = gvf(empt, event_obj)

global Srf Model dcm_obj handles ph

% Selected position in mesh space coordinates
pos = get(event_obj, 'Position'); 
% Selected "vertex" number
v = find(Srf.Vertices(:,1) == pos(1) & Srf.Vertices(:,2) == pos(2), 1); 
% No data string
txt = ''; 
if ~isempty(v)
    s = Srf.SensorNums(v); % Sensor number
    t = nanmin(Srf.TimePts) + get(handles.slider1,'Value')*range(Srf.TimePts); % Time point
    x = find(Srf.SensorNums==s & Srf.TimePts == t); % Corresponding column index
    m = get(handles.listbox1,'Value'); % Current value map
    b = find(strcmpi(Srf.Values, 'Beta')); % Which row is beta?
    % Plot figure
    figure(ph);
    if isfield(Srf, 'X_coords')
        % Plot reverse correlation pRF profile
        samsrf_showprf(Srf, x);
        hold on
        scatter(Srf.Data(2,x), Srf.Data(3,x), 120, 'k', 'filled', 'markeredgecolor', 'w');
        scatter(Srf.Data(2,x), Srf.Data(3,x), 60, '*w');
        hold off
        caxis([-1 1]*nanmean(abs(Srf.Data(b, Srf.TimePts == t))));
    elseif isfield(Model, 'Param1')
        % Plot forward model pRF profile
        samsrf_showprf(Srf, x, Model);
    else
        % Plot response values
        plot(Srf.Data(:,x), 'ko-', 'MarkerFaceColor', 'w');
        grid on
        set(gca, 'xtick', 1:length(Srf.Values), 'xticklabels', Srf.Values);
        r = Srf.Data(:, Srf.TimePts == t);
        if sum(Srf.Data(:) < 0)
            % Voltage data
            yl = [-1 1];
        else
            % Power data
            yl = [0 1];
        end
        ylim(yl*nanmax(abs(r(:))));
        ylabel('Response');
    end
    % Information
    title([Srf.Sensors{s} ': ' Srf.Values{m} ' = ' num2str(Srf.Data(m,x))]);
end
dcm_obj.removeAllDataCursors();


%% Closing function
function closereq(src, evnt)
% Clear global variables when figure is closed

clear global Srf Model dcm_obj handles ph 
delete(src);
