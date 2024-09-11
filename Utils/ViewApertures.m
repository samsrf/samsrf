function varargout = ViewApertures(varargin)
% VIEWAPERTURES MATLAB code for ViewApertures.fig
%      VIEWAPERTURES, by itself, creates a new VIEWAPERTURES or raises the existing
%      singleton*.
%
%      H = VIEWAPERTURES returns the handle to a new VIEWAPERTURES or the handle to
%      the existing singleton*.
%
%      VIEWAPERTURES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWAPERTURES.M with the given input arguments.
%
%      VIEWAPERTURES('Property','Value',...) creates a new VIEWAPERTURES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewApertures_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewApertures_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 03/08/2024 - Added support for multi-condition apertures (DSS)
% 05/09/2024 - Now allows direct input of filename (DSS)
%

% Last Modified by GUIDE v2.5 23-May-2018 15:22:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewApertures_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewApertures_OutputFcn, ...
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


% --- Executes just before ViewApertures is made visible.
function ViewApertures_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewApertures (see VARARGIN)

% Choose default command line output for ViewApertures
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ViewApertures wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global ApFrm ApXY ApCond

% Load apertures
ApCond = [];
if isempty(varargin) 
    [ApName, ApPath] = uigetfile('aps_*.mat', 'Select apertures');
else
    ApName = [varargin{1} '.mat'];
    ApPath = '';
end
if ApName ~= 0
    ApName = ApName(1:end-4);
    load([ApPath ApName]);
end

% Plot first frame
subplot(handles.axes1);
if ~isempty(ApXY)
    % Vectorised apertures
    disp('Loading vectorised apertures...');
    if size(ApFrm,1) == size(ApXY,1)+1
        disp(' Multi-condition apertures!');
        ApCond = ApFrm(1,:);
        ApFrm = ApFrm(2:end,:);
    end
    polarpatch(ApXY(:,1), ApXY(:,2), ApFrm(:,1));
    if isempty(ApCond)
        txy = max(ApXY) .* [.65 .85];
    else
        txy = max(ApXY) .* [.45 .85];
    end
else
    % Movie apertures
    ApXY = []; 
    disp('Loading movie apertures...');
    [x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
    y = flipud(y);
    curfrm = ApFrm(:,:,1);
    if var(curfrm(:)) == 0
        imshow(curfrm); 
    else
        contourf(x, y, curfrm, 'edgecolor', 'none');
    end
    if isempty(ApCond)
        txy = [80 10];
    else
        txy = [60 10];
    end
end
axis square
axis off
if nanmin(ApFrm(:)) < 0
    colormap berlin
    set(gca, 'Clim', [-1 1]);
else
    colormap gray
    set(gca, 'Clim', [0 1]);
end
if ~isempty(ApXY)
    set(handles.slider1, 'Value', 1/size(ApFrm,2), 'SliderStep', [1/(1+size(ApFrm,2)) 1/(1+size(ApFrm,2))]);
else
    set(handles.slider1, 'Value', 1/size(ApFrm,3), 'SliderStep', [1/(1+size(ApFrm,3)) 1/(1+size(ApFrm,3))]);
end
vc = '1';
if ~isempty(ApCond)
    vc = [vc ' (' num2str(ApCond(1)) ')'];
end
text(txy(1), txy(2), vc, 'color', [.5 .5 .5], 'fontsize', 10);

% Close request function
crfcn = @closereq;
set(hObject, 'CloseRequestFcn', crfcn);


% --- Outputs from this function are returned to the command line.
function varargout = ViewApertures_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global ApFrm ApXY ApCond

v = get(hObject,'Value');
if ~isempty(ApXY)
    f = ceil(v*size(ApFrm,2));
else
    f = ceil(v*size(ApFrm,3));
end
if f == 0 f = 1; end
subplot(handles.axes1);
if ~isempty(ApXY)
    % Vectorised apertures
    polarpatch(ApXY(:,1), ApXY(:,2), ApFrm(:,f));
    if isempty(ApCond)
        txy = max(ApXY) .* [.65 .85];
    else
        txy = max(ApXY) .* [.45 .85];
    end
else
    % Movie apertures
    [x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
    y = flipud(y);
    curfrm = ApFrm(:,:,f);
    if var(curfrm(:)) == 0
        imshow(curfrm); 
    else
        contourf(x, y, curfrm, 'edgecolor', 'none');
    end
    if isempty(ApCond)
        txy = [80 10];
    else
        txy = [60 10];
    end
end
axis square
axis off
if nanmin(ApFrm(:)) < 0
    colormap berlin
    set(gca, 'Clim', [-1 1]);
else
    colormap gray
    set(gca, 'Clim', [0 1]);
end
vc = num2str(f);
if ~isempty(ApCond)
    vc = [vc ' (' num2str(ApCond(f)) ')'];
end
text(txy(1), txy(2), vc, 'color', [.5 .5 .5], 'fontsize', 10);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

global ApFrm ApXY ApCond

try
    while get(hObject, 'Value')
        v = get(handles.slider1,'Value'); % Get slider value
        if ~isempty(ApXY)
            % Vectorised apertures
            f = ceil(v*size(ApFrm,2)); % Convert to frame number
        else
            % Movie apertures
            f = ceil(v*size(ApFrm,3)); % Convert to frame number
        end
        f = f + 1; % Advance by one frame
        if ~isempty(ApXY)
            % Vectorised apertures
            if f > size(ApFrm,2)
                f = 1;
            end
            v = f / size(ApFrm,2); % Convert back to slider value
        else
            % Movie apertures 
            if f > size(ApFrm,3)
                f = 1;
            end
            v = f / size(ApFrm,3); % Convert back to slider value
        end
        set(handles.slider1, 'Value', v);

        % Plot next frame
        subplot(handles.axes1);
        if ~isempty(ApXY)
            % Vectorised apertures
            polarpatch(ApXY(:,1), ApXY(:,2), ApFrm(:,f));
            if isempty(ApCond)
                txy = max(ApXY) .* [.65 .85];
            else
                txy = max(ApXY) .* [.45 .85];
            end
        else
            % Movie apertures
            [x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
            y = flipud(y);
            curfrm = ApFrm(:,:,f);
            if var(curfrm(:)) == 0
                imshow(curfrm); 
            else
                contourf(x, y, curfrm, 'edgecolor', 'none');
            end
            if isempty(ApCond)
                txy = [80 10];
            else
                txy = [60 10];
            end
        end
        axis square
        axis off
        if nanmin(ApFrm(:)) < 0
            colormap berlin
            set(gca, 'Clim', [-1 1]);
        else
            colormap gray
            set(gca, 'Clim', [0 1]);
        end
        vc = num2str(f);
        if ~isempty(ApCond)
            vc = [vc ' (' num2str(ApCond(f)) ')'];
        end
        text(txy(1), txy(2), vc, 'color', [.5 .5 .5], 'fontsize', 10);
        pause(0.1);
    end
catch
end


% Clear global ApFrm when figure is closed
function closereq(src, evnt)
clear global 
delete(src);
close all