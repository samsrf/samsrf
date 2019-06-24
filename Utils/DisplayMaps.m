function varargout = DisplayMaps(varargin)
% 
% ******* SamSrf pRF map display tool ******* 
%
% Run this GUI to display maps and make figures from them.
%
%  If you wish to change default ROI to limit the mesh, 
%  change them below where indicated by '%%% %%% %%%'.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplayMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplayMaps_OutputFcn, ...
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


%% --- Executes just before DisplayMaps is made visible.
function DisplayMaps_OpeningFcn(hObject, eventdata, handles, varargin)
global Srf RoiName Pval

%% Default parameters (Feel free to edit)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
RoiName = 'occ'; % ROI name to limit the mesh (without hemisphere prefix)
Pval = 0.00000001; % Starting p-value with which to threshold maps
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

%% Welcome message
[vn vd] = samsrf_version;
clc; 
disp('****************************************************************************');
disp('     Welcome to the Seriously Annoying MatLab Surfer Map Display Tool!');
disp('     by D. S. Schwarzkopf from the University of Auckland, New Zealand');
new_line;
disp(['                 Version ' num2str(vn) ', Released on ' vd]);
disp('      (see SamSrf/ReadMe.txt for what is new in this version)');
disp('****************************************************************************');
new_line;

% Choose default command line output for DisplayMaps
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%% Load first maps
LoadMapData(handles, RoiName);
% Draw map
RedrawMaps(handles, true);

%% Close request function
crfcn = @SamSrfCloseReq;
set(hObject, 'CloseRequestFcn', crfcn);


%% Output function
function varargout = DisplayMaps_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%% Load map file
function LoadMapData(handles, RoiName)
global Srf Pval R2Thrsh

%% Load surface data file
[fn,pn] = uigetfile('*.mat', 'Load Srf');
load([pn filesep fn]);
Srf = samsrf_expand_srf(Srf);
vx = samsrf_loadlabel([pn filesep Srf.Hemisphere '_' RoiName]);
% Remove non-ROI vertices if desired
if ~isnan(vx)
    disp(['Only displaying ROI ' Srf.Hemisphere '_' RoiName '!']);
    roivx = true(size(Srf.Vertices,1),1);
    roivx(vx) = false;
    Srf.Vertices(roivx,:) = NaN;
end
% Turn off toggle button if no smoothed maps
if isfield(Srf, 'Raw_Data')
    Srf.Smoothed_Data = Srf.Data;
    set(handles.togglebutton1, 'Enable', 'on');
else
    set(handles.togglebutton1, 'Enable', 'off');
end

%% Update map list
% If polar & eccentricity need to be calculated
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
    Values = Srf.Values;
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
else
    Values = Srf.Values;
end
% Set to first item in the list
set(handles.popupmenu2, 'String', Values, 'Value', 1); 

%% Update R^2 threshold
if isfield(Srf, 'Y')  
    % If time course is saved
    Tval = tinv(1 - Pval/2, size(Srf.Y,1)-2); % t-statistic for this p-value
    R = sign(Tval) .* (abs(Tval) ./ sqrt(size(Srf.Y,1)-2 + Tval.^2)); % Convert t-testistic into correlation coefficient
    R2Thrsh = R.^2; % Variance explained 
    set(handles.edit1, 'String', num2str(R2Thrsh)); % Set threshold
else
    % Just default to zero
    set(handles.edit1, 'String', '0'); % Set threshold
end

%% Redraw maps
function RedrawMaps(handles, IsNewFig)
global Srf Paths FigHdl PatchHdl 

% Reopen, if figure was closed 
if ~ishandle(FigHdl)
    IsNewFig = true;
end
% Which surface mesh?
contents = cellstr(get(handles.popupmenu1,'String'));
Mesh = contents{get(handles.popupmenu1,'Value')};
% Which map? 
contents = cellstr(get(handles.popupmenu2,'String'));
MapNum = get(handles.popupmenu2,'Value');
MapType = contents{MapNum};
% What R^2 threshold?
if strcmpi(Srf.Values{1}, 'R^2')
    R2Thrsh = str2double(get(handles.edit1, 'String'));
    if strcmpi(get(handles.togglebutton1, 'Enable'), 'on')
        if get(handles.togglebutton1,'Value') 
            Srf.Data = Srf.Raw_Data;
            Srf.Data(1,Srf.Data(1,:) <= R2Thrsh) = NaN;
        else
            Srf.Data = Srf.Smoothed_Data;
            Srf.Data(1,Srf.Data(1,:) <= R2Thrsh) = NaN;
        end
    end
else
    % No R^2 in this Srf
    set(handles.edit1, 'String', '0', 'Enable', 'off');
    R2Thrsh = 0;
end
% Cut-offs for maps
CutOffStr = get(handles.edit2, 'String');
ThrshVec = eval(['[' CutOffStr ']']);
ThrshVec = [R2Thrsh ThrshVec];
% Eccentricity range
ThrshVec = [ThrshVec eval(['[' get(handles.edit4, 'String') ']'])];
% Camera angle
if IsNewFig
    FigHdl = figure;
    if ~exist('def_views', 'var')
        % Focus on early visual cortex
        if Srf.Hemisphere(1) == 'l'
            CamView = [36 -20 1.8];
        else
            CamView = [-50 -6 1.6];
        end
    else
        % Use default camera angle
        if Srf.Hemisphere(1) == 'l'
            CamView = def_views(:,1)';
        else
            CamView = def_views(:,2)';
        end
    end
else
    figure(FigHdl);
    CamView = [get(gca, 'View') 1];
end
% Draw the map
if get(handles.popupmenu2,'Value') <= size(Srf.Data,1) || strcmpi(MapType, 'Polar') || strcmpi(MapType, 'Eccentricity')
    figure(FigHdl);
    if IsNewFig
        % Draw new mesh
        PatchHdl = samsrf_surf(Srf, Mesh, ThrshVec, Paths, CamView, MapNum);
        set(gcf, 'Position', [306  501  560  420]);
        set(gca,'CameraViewAngleMode', 'Manual');
        daspect([1 1 1]);
    else
        % Just update colours
        samsrf_surf(Srf, Mesh, ThrshVec, Paths, CamView, MapNum, PatchHdl);
    end
    FigHdl = gcf;
else
    figure(FigHdl);
    if IsNewFig
        % Draw new mesh
        PatchHdl = samsrf_surf(Srf, Mesh, Inf, Paths, CamView, Srf.Values{1});
        set(gcf, 'Position', [306  501  560  420]);
        set(gca,'CameraViewAngleMode', 'Manual');
        daspect([1 1 1]);
    else
        % Just update colours
        samsrf_surf(Srf, Mesh, Inf, Paths, CamView, Srf.Values{1}, PatchHdl);
    end
    FigHdl = gcf;
    set(FigHdl, 'Name', 'Map undefined for unsmoothed data!');
end


%% Clear global variables when figure is closed
function SamSrfCloseReq(src, evnt)
new_line;
disp('SamSrf says:'); 
disp(' "I hope you''ll come back to be annoyed by me again soon...');
disp('  Good bye for now!"');
new_line;

clear global 
delete(src);


%% Draw ROIs or paths
function pushbutton1_Callback(hObject, eventdata, handles)
global Paths

[pf, pn] = uigetfile('*.label; *path*; del_*.mat', 'MultiSelect', 'on');
if isscalar(pf) && pf == 0
    pf = {};
    Paths = {};
else
    if ischar(pf)
        pf = {pf};
    end

    Paths = {};
    for i = 1:length(pf)
        Paths{i} = [pn filesep pf{i}];
    end
    Paths = Paths';
end
RedrawMaps(handles, false);


%% Smoothed or raw maps?
function togglebutton1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton1
global Srf

if get(hObject,'Value') 
    % Use raw data
    Srf.Data = Srf.Raw_Data;
    set(hObject, 'String', 'Smoothed Maps');
else
    % Use smoothed data
    Srf.Data = Srf.Smoothed_Data;
    set(hObject, 'String', 'Unsmoothed Maps');
end
% Update map list
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0')
    Values = Srf.Values;
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
else
    Values = Srf.Values;
end
set(handles.popupmenu2, 'String', Values);
RedrawMaps(handles, false);


%% Which surface mesh?
function popupmenu1_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global FigHdl

if ishandle(FigHdl)
    close(FigHdl);
end
RedrawMaps(handles, false);

%% 
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% R^2 threshold 
function edit1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
RedrawMaps(handles, false);

%%
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Cut-off values (min/max for most maps)
function edit2_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% If undefined set to R^2 threshold
if isempty(get(hObject, 'String'))
    R2Thrsh = get(handles.edit1, 'String'); 
    set(hObject, 'String', R2Thrsh);  
end
RedrawMaps(handles, false);

%%
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Which map to display?
function popupmenu2_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

% Which map? 
contents = cellstr(get(handles.popupmenu2,'String'));
MapType = contents{get(handles.popupmenu2,'Value')};
% Default cut-offs
if strcmpi(MapType, 'R^2')
    set(handles.edit2, 'String', '0 1');
elseif strcmpi(MapType, 'Polar')
    set(handles.edit2, 'String', '0 0');
elseif strcmpi(MapType, 'Eccentricity') 
    if exist('SamSrf_defaults.mat', 'file')
        load('SamSrf_defaults.mat');
        set(handles.edit2, 'String', ['0 ' num2str(def_eccen)]);
    else
        set(handles.edit2, 'String', '0 10');
    end
elseif strcmpi(MapType, 'Sigma') || strcmpi(MapType, 'Fwhm') || strcmpi(MapType, 'nSigma') ...
       || strcmpi(MapType, 'Centre') || strcmpi(MapType, 'Surround') || strcmpi(MapType, 'Sigma1') || strcmpi(MapType, 'Sigma2')   
    if exist('SamSrf_defaults.mat', 'file')
        load('SamSrf_defaults.mat');
        set(handles.edit2, 'String', ['0 ' num2str(def_eccen/2)]);
    else
        set(handles.edit2, 'String', '0 5');
    end
elseif strcmpi(MapType, 'Spread') 
    set(handles.edit2, 'String', '0 50');
elseif strcmpi(MapType, 'Beta') || strcmpi(MapType, 'Suppression')
    set(handles.edit2, 'String', '0 0.5');
elseif strcmpi(MapType, 'Cmf')
    set(handles.edit2, 'String', '0 10');
elseif strcmpi(MapType, 'Field Sign')
    set(handles.edit2, 'String', '0 3');
elseif strcmpi(MapType, 'x0') || strcmpi(MapType, 'y0')
    set(handles.edit2, 'String', '0 10');
else
    set(handles.edit2, 'String', '0 0.5');
end
RedrawMaps(handles, false);

%%
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Open new figure
function pushbutton4_Callback(hObject, eventdata, handles)

RedrawMaps(handles, true);


%% Clone cameras in all figures
function pushbutton5_Callback(hObject, eventdata, handles)
global FigHdl

figure(FigHdl);
samsrf_clonecam;


%% Eccentricity range
function edit4_Callback(hObject, eventdata, handles)

RedrawMaps(handles, false);

%%
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Load new map
function pushbutton6_Callback(hObject, eventdata, handles)
global RoiName FigHdl

if ishandle(FigHdl)
    close(FigHdl);
end
set(handles.edit2, 'String', '');
set(handles.edit4, 'String', '');
set(handles.edit2, 'String', '0 1');
LoadMapData(handles, RoiName);
RedrawMaps(handles, true);
