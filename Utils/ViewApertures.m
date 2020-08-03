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
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewApertures

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

global ApFrm

% Load apertures
[ApName, ApPath] = uigetfile('aps_*.mat', 'Select apertures');
if ApName ~= 0
    ApName = ApName(1:end-4);
    load([ApPath ApName]);
end

% Plot first frame
subplot(handles.axes1);
[x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
y = flipud(y);
curfrm = ApFrm(:,:,1);
if var(curfrm(:)) == 0
    imshow(curfrm); 
else
    contourf(x, y, curfrm, 'edgecolor', 'none');
end
axis square
axis off
if nanmin(ApFrm(:)) < 0
    colormap hotcold
    set(gca, 'Clim', [-1 1]);
else
    colormap gray
    set(gca, 'Clim', [0 1]);
end
set(handles.slider1, 'Value', 1/size(ApFrm,3), 'SliderStep', [1/(1+size(ApFrm,3)) 1/(1+size(ApFrm,3))]);
text(80, 10, '1', 'color', [.5 .5 .5], 'fontsize', 15);

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

global ApFrm

v = get(hObject,'Value');
f = ceil(v*size(ApFrm,3));
if f == 0 f = 1; end
subplot(handles.axes1);
[x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
y = flipud(y);
curfrm = ApFrm(:,:,f);
if var(curfrm(:)) == 0
    imshow(curfrm); 
else
    contourf(x, y, curfrm, 'edgecolor', 'none');
end
axis square
axis off
if nanmin(ApFrm(:)) < 0
    colormap hotcold
    set(gca, 'Clim', [-1 1]);
else
    colormap gray
    set(gca, 'Clim', [0 1]);
end
text(80, 10, num2str(f), 'color', [.5 .5 .5], 'fontsize', 15);


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

global ApFrm

try
    while get(hObject, 'Value')
        v = get(handles.slider1,'Value'); % Get slider value
        f = ceil(v*size(ApFrm,3)); % Convert to frame number
        f = f + 1; % Advance by one frame
        if f > size(ApFrm,3)
            f = 1;
        end
        v = f / size(ApFrm,3); % Convert back to slider value
        set(handles.slider1, 'Value', v);

        % Plot next frame
        subplot(handles.axes1);
        [x,y] = meshgrid(1:size(ApFrm,2), 1:size(ApFrm,1));
        y = flipud(y);
        curfrm = ApFrm(:,:,f);
        if var(curfrm(:)) == 0
            imshow(curfrm); 
        else
            contourf(x, y, curfrm, 'edgecolor', 'none');
        end
        axis square
        axis off
        if nanmin(ApFrm(:)) < 0
            colormap hotcold
            set(gca, 'Clim', [-1 1]);
        else
            colormap gray
            set(gca, 'Clim', [0 1]);
        end
        text(80, 10, num2str(f), 'color', [.5 .5 .5], 'fontsize', 15);
        pause(0.1);
    end
catch
end


% Clear global ApFrm when figure is closed
function closereq(src, evnt)
clear global 
delete(src);
close all