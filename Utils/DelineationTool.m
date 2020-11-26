function varargout = DelineationTool(varargin)
% 
% ******* SamSrf Delineation tool ******* 
%
% Run this GUI to delineate regions of interest in retinotopic maps.
%  (Sorry, this tool currently does not cater for other sensory modalities).
%
%  If you wish to change default parameters change them in this script.
%  Look for the section indicated by '%%% %%% %%%'.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DelineationTool_OpeningFcn, ...
                   'gui_OutputFcn',  @DelineationTool_OutputFcn, ...
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


% --- Executes just before DelineationTool is made visible.
function DelineationTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DelineationTool (see VARARGIN)

% Choose default command line output for DelineationTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DelineationTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% Global variables
global SrfName Roi Srf R2Thrsh EccThrsh Curv Vertices Points Paths RoiList RoiSeeds RoiColours PolRgb EccRgb FsRgb CfRgb hp he hf hc pp pe pf pc IsFsMap ActPrct

load('SamSrf_defaults.mat');
if ~exist('def_disproi')
    def_disproi = NaN; 
end

%% Default parameters (Change at your own leisure/peril!)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
Mesh = 'sphere'; % Which cortex model to use
RoiName = def_disproi; % ROI name without hemisphere
Pval = 0.0001; % Starting p-value with which to threshold maps
ActPrct = [5 95]; % Percentiles to threshold activation maps
RoiList = {'V1' 'V2v' 'V3v' 'V4' 'V2d' 'V3d' 'V3A' 'V3B' 'LO1' 'LO2' 'VO1' 'VO2' 'TO1' 'TO2' 'V6' 'IPS0' 'IPS1' 'IPS2'}'; % ROI list
RoiColours = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [0 1 0]; [0 0 1]; hsv(6); jet(6)]; % Colours of ROIs
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
 
%% Initialise variousn variables
Points = []; % Initialise way point vector
Paths = {}; % Initialise path vector
RoiSeeds = NaN(length(RoiList),1); % Initialise list of ROI seed vertices
IsFsMap = true; % Always start with field sign map
% Eccentricity range
EccStr = get(handles.edit2,'String');
EccVec = eval(['[' EccStr ']']);
EccThrsh = EccVec(1:2); 

%% Load map data & mesh folder & region of interest
% Load data
SrfName = uigetfile('*h_*.mat', 'Select pRF map');
SrfName = SrfName(1:end-4);
load(SrfName);
[Srf, vx] = samsrf_expand_srf(Srf);

% What ROI to display
if RoiName(1) == '<' || RoiName(1) == '>'
    % If ROI defined by coordinates
    if length(RoiName) == 1
        error('You must define inflated mesh coordinate in def_disproi!');
    end
    switch RoiName(2)
        case 'X'
            wv = Srf.Inflated(:,1);
        case 'Y'
            wv = Srf.Inflated(:,2);
        case 'Z'
            wv = Srf.Inflated(:,3);
        otherwise
            error('Invalid inflated mesh coordinate specified in def_disproi!');
    end
    if length(RoiName) < 3
        error('You must define inflated mesh cut-off coordinate in def_disproi!');
    end
    wc = str2double(RoiName(3:end));
    if RoiName(1) == '<'
        vx = find(wv < wc);
    elseif RoiName(1) == '>'
        vx = find(wv > wc);
    end
    disp(['Only displaying inflated mesh vertices with ' RoiName]);
elseif isnan(RoiName)
    % Use ROI from Srf
    disp('Using ROI from Srf if it exists');
else
    % Load region of interest
    Roi = [SrfName(1:2) '_' RoiName];
    vx = samsrf_loadlabel(Roi);
    disp(['Only displaying ROI ' Roi]);
end
if ~isnan(vx)
    roivx = true(size(Srf.Vertices,1),1);
    roivx(vx) = false;
    Srf.Vertices(roivx,:) = NaN;
end
if ~isfield(Srf, 'Values')
    disp('This file does not contain a pRF map!');
    return
end
Srf.Values{end+1} = 'Custom';
if sum(strcmpi(Srf.Values, 'Beta')) == 1
    Srf.Data(end+1,:) = Srf.Data(strcmpi(Srf.Values, 'Beta'),:);
end
% Check if field sign exists
if sum(strcmpi(Srf.Values, 'Field Sign')) == 0
    IsFsMap = false;
    FsMapName = 'Beta map';
else
    FsMapName = 'Field Sign map';
end

% Determine starting threshold
if isfield(Srf, 'Y')  
    % If time course is saved
    Tval = tinv(1 - Pval/2, size(Srf.Y,1)-2); % t-statistic for this p-value
    R = sign(Tval) .* (abs(Tval) ./ sqrt(size(Srf.Y,1)-2 + Tval.^2)); % Convert t-testistic into correlation coefficient
    R2 = R.^2; % Variance explained 
else
    % Just default to zero
    R2 = 0;
end
set(handles.edit1, 'String', num2str(R2)); % Set threshold box
R2Thrsh = str2double(get(handles.edit1, 'String')); % R^2 threshold

%% Load curvature
if isfield(Srf, 'Curvature') 
    Curv = Srf.Curvature;
else
    disp('No curvature data!');
    Curv = 1;
end
% Greyscale for curvature
CurvGrey = gray(11); % Grey scale colour map
CurvGrey = CurvGrey(1:10,:); % Remove black & white
% Transform curvature
Curv = -Curv + 0.5;
Curv(Curv <= 0) = 0.000001;
Curv(Curv > 1) = 1;
Curv = ceil(Curv * size(CurvGrey,1));

%% Load mesh
Mesh = lower(Mesh);
Mesh(1) = upper(Mesh(1));
if isfield(Srf, Mesh)
    Vertices = getfield(Srf, Mesh);
else
    disp(['Mesh ' Mesh ' not found!']);
    Vertices = Srf.Vertices;
end
Vertices(isnan(Srf.Vertices(:,1)),:) = NaN; 
Faces = Srf.Faces;
if length(Curv) == 1
    Curv = repmat(Curv, 1, size(Vertices,1));
end

%% R^2 threshold
set(handles.edit1,'String',num2str(R2Thrsh));
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
    r = Srf.Data(1,:) <= 0.01 | isnan(Srf.Data(1,:));
else 
    r = isnan(Srf.Data(1,:));
end

%% Display maps
CalculateMapRgb;
[hc pc] = DisplayMesh('Cortex & ROIs', Srf, Vertices, Faces, CfRgb);
[hf pf] = DisplayMesh(FsMapName, Srf, Vertices, Faces, FsRgb);
[he pe] = DisplayMesh('Eccentricity map', Srf, Vertices, Faces, EccRgb);
[hp pp] = DisplayMesh('Polar map', Srf, Vertices, Faces, PolRgb);
set(hc, 'Units', 'normalized');
fpos = get(hc, 'Position');
set(hc, 'Units', 'normalized', 'Position', [fpos(3)*0.25 fpos(2:4)]);
set(hf, 'Units', 'normalized', 'Position', [fpos(3)*0.75 fpos(2:4)]);
set(he, 'Units', 'normalized', 'Position', [fpos(3)*1.25 fpos(2:4)]);
set(hp, 'Units', 'normalized', 'Position', [fpos(3)*1.75 fpos(2:4)]);

%% Welcome message
[vn vd] = samsrf_version;
clc; 
disp('****************************************************************************');
disp('   Welcome to the Seriously Annoying MatLab Surfer Map Delineation Tool!');
disp('    by D. S. Schwarzkopf from the University of Auckland, New Zealand');
new_line;
disp(['                 Version ' num2str(vn) ', Released on ' vd]);
disp('      (see SamSrf/ReadMe.md for what is new in this version)');
disp('****************************************************************************');
new_line;
disp(['Displaying map: ' SrfName]);
new_line;

%% Close request function
crfcn = @closereq;
set(hObject, 'CloseRequestFcn', crfcn);


%% Display a map on the mesh
function [h p] = DisplayMesh(Name, Srf, Vertices, Faces, Rgb)
% Open figure with only permitted tools
h = figure('Name', Name, 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure');
tb = findall(h);
for i = 3:length(tb) 
    set(tb(i), 'Visible', 'off'); 
end
set(findall(tb,'ToolTipString','Data Cursor'), 'Visible', 'on');
set(findall(tb,'ToolTipString','Rotate 3D'), 'Visible', 'on');
% Display the mesh
p = patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Rgb, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis off;
ax = gca;
ax.Clipping = 'off';
% Focus on early visual cortex
if Srf.Hemisphere(1) == 'l'
    set(gca, 'view', [47 -32]);
else
    set(gca, 'view', [-30 -40]);
end
zoom(2);
set(gca, 'projection', 'perspective');

% Ensure camera position is constant
set(gca,'CameraViewAngleMode', 'Manual');
daspect([1 1 1]);

% Coupled rotations
rc = rotate3d;
rc.ActionPostCallback = @postrotcallback;
rc.Enable = 'on';

% Vertex selection hack
dcm_obj_cf = datacursormode(h);
set(dcm_obj_cf, 'UpdateFcn', @getvertexfun);
% THE FOLLOWING LINE NEEDS TO BE UNCOMMENTED IF YOU ALWAYS GET A TOOL TIP WHEN YOU CLICK ON ANY VERTEX!!!
% set(dcm_obj_cf, 'DisplayStyle', 'Window');

% Close request function
crfcn = @closereq;
set(h, 'CloseRequestFcn', crfcn);


%% Clears all data cursors everywhere
function ClearDataCursors
global hp hf he hc
dcm_obj = datacursormode(hp);
dcm_obj.removeAllDataCursors();
dcm_obj = datacursormode(hf);
dcm_obj.removeAllDataCursors();
dcm_obj = datacursormode(he);
dcm_obj.removeAllDataCursors();
dcm_obj = datacursormode(hc);
dcm_obj.removeAllDataCursors();


%% Calculate RGB colours for the maps in the delineation tool 
function CalculateMapRgb(LoadSavedDelins) 
global SrfName Srf Paths R2Thrsh EccThrsh Curv CfRgb PolRgb EccRgb FsRgb RoiList RoiSeeds RoiColours pc IsFsMap ActPrct

if nargin == 0
    LoadSavedDelins = true;
end

% Remove rubbish
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
    r = Srf.Data(1,:) <= 0.01 | isnan(Srf.Data(1,:));
else 
    r = isnan(Srf.Data(1,:));
end

% Greyscale for curvature
figure; CurvGrey = gray(11); close; % Grey scale colour map
CurvGrey = CurvGrey(1:10,:); % Remove black & white

% Colours of paths
WhitePath = [1 1 1];
BlackPath = [0 0 0];
YellowPath = [1 1 0]; 

% Vertices on paths
Vs = [];
if LoadSavedDelins
    % Load any previously saved paths
    if exist(['del_' SrfName '.mat'], 'file')
        Vs = samsrf_loadpaths(['del_' SrfName '.mat']);
        disp('Loaded previously saved paths.');
    end
else
    % Retrieve path vertices
    for i = 1:length(Paths)
        Vs = [Vs; Paths{i}];
    end
end

% Eccentricity map
Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
Rho(Srf.Data(1,:) <= R2Thrsh) = NaN;
Rho(Rho < EccThrsh(1)) = NaN;
Rho = Rho / EccThrsh(2);
Rho(Rho > 1) = NaN; 
Pha = round(Rho * 360);
Pha = mod(Pha, 360);
Pha(Pha==0) = 360;
Pha(r|isnan(Pha)) = 360 + Curv(r|isnan(Pha));
figure; Cmap = [colormap(fscol); CurvGrey]; close
EccRgb = Cmap(Pha,:);
EccRgb(Vs,:) = repmat(BlackPath, size(Vs,1), 1); % Draw paths

% Polar map
Pha = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
Pha(Srf.Data(1,:) <= R2Thrsh | isnan(Rho)) = NaN;
Cmap = colormap(fscol); 
Cmap = [Cmap; CurvGrey];
if Srf.Hemisphere(1) == 'l'
    Pha = mod(ceil(Pha + 270), 360) + 1;
else
    Pha = -Pha;
    Pha = mod(ceil(Pha + 90), 360) + 1;
end
Pha(Pha == 0) = 360;
Pha(r|isnan(Pha)) = 360 + Curv(r|isnan(Pha));
PolRgb = Cmap(Pha,:);
PolRgb(Vs,:) = repmat(WhitePath, size(Vs,1), 1); % Draw paths

if sum(strcmpi(Srf.Values, 'Field Sign')) == 0 
    IsFsMap = false;
end
if IsFsMap
    % Field sign map
    Fs = sign(Srf.Data(strcmpi(Srf.Values, 'Field Sign'), :));
    Fs(Srf.Data(1,:) <= R2Thrsh | isnan(Rho)) = 0;
    FsRgb = CurvGrey(Curv,:);
    FsRgb(Fs'==-1,:) = (FsRgb(Fs'==-1,:) + repmat([0 0 1],sum(Fs==-1),1)) / 2;
    FsRgb(Fs'==+1,:) = (FsRgb(Fs'==+1,:) + repmat([1 0 0],sum(Fs==+1),1)) / 2;
    FsRgb(Vs,:) = repmat(YellowPath, size(Vs,1), 1); % Draw paths
else
    % Activation data
    Fs = Srf.Data(end, :);
    Fs(Srf.Data(1,:) <= R2Thrsh | isnan(Rho)) = NaN;
    Thrsh = prctile(abs(Fs), ActPrct);
    Fs(Fs>0 & Fs<Thrsh(1)) = Thrsh(1);
    Fs(Fs>0 & Fs>Thrsh(2)) = Thrsh(2);
    Fs(Fs<0 & Fs>-Thrsh(1)) = -Thrsh(1);
    Fs(Fs<0 & Fs<-Thrsh(2)) = -Thrsh(2);
    Fs(Fs>0) = Fs(Fs>0) - Thrsh(1);
    Fs(Fs<0) = Fs(Fs<0) + Thrsh(1);
    mThr = Thrsh - Thrsh(1);
    Pha = round(Fs / abs(mThr(2)) * 100) + 100;
    Cmap = [flipud(colormap(winter(100))); colormap(hot(100)); CurvGrey];
    Pha(Pha==0) = 1;
    Pha(r|isnan(Pha)|Pha==100) = 200 + Curv(r|isnan(Pha)|Pha==100);
    FsRgb = Cmap(Pha,:);
    FsRgb(Vs,:) = repmat(BlackPath, size(Vs,1), 1); % Draw paths
end

% Cortex & ROIs
Cmap = [CurvGrey; RoiColours];
CfRgb = Cmap(Curv,:);
% Load ROI labels
if exist(['del_' SrfName '.mat'], 'file')
    for i = 1:length(RoiList)
        if ~isnan(RoiSeeds(i))
            % Load label
            Rvs = samsrf_loadlabel(['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_' RoiList{i}]);
            if isnan(Rvs)
                % ROI label no longer exists
                RoiSeeds(i) = NaN; % Remove from labelled list
            else
                % Paint the map
                Rgb = (repmat(Cmap(i+size(CurvGrey,1),:),length(Rvs),1) + CfRgb(Rvs,:)) / 2; % What colour?
                CfRgb(Rvs,:) = Rgb; % Colour the ROI
                % Update map
                set(pc, 'FaceVertexCData', CfRgb); 
            end
        end
    end
end
CfRgb(Vs,:) = repmat(WhitePath, size(Vs,1), 1); % Draw paths

% Display paths on maps
RedrawAllPaths;


%% Draw a path between two selected vertices
function DrawPaths
global Srf Vertices Points Paths PolRgb EccRgb FsRgb CfRgb pp pe pf pc

% Colours of paths
WhitePath = [1 1 1];
BlackPath = [0 0 0];
YellowPath = [1 1 0]; 

% Create connected path
Paths{length(Paths) + 1} = []; % New path index
if length(Points) > 1
    for p = 2:length(Points)
        pv = Points(p-1); % Previous vertex
        cv = Points(p); % Current vertex
        p2c = Vertices(cv,:) - Vertices(pv,:); % Vector from previous to current vertex
        vsp = []; % Vertices on new path
        for s = 0:.01:1
            cs = Vertices(pv,:) + p2c*s; % One small step along vector
            d = sqrt((Vertices(:,1)-cs(1)).^2 + (Vertices(:,2)-cs(2)).^2 + (Vertices(:,3)-cs(3)).^2); % Euclidian distance of all vertices from current step
            nv = find(d == min(d), 1); % New vertex on path
            nvn = samsrf_neighbours(nv, Srf.Faces);
            vsp = [vsp; nv; nvn(1:2:end)];
        end
        Paths{end} = [Paths{end}; vsp];
    end
end
Points = [];
ClearDataCursors;

% Colour vertices on path
PolRgb(Paths{end},:) = repmat(WhitePath, length(Paths{end}), 1);
EccRgb(Paths{end},:) = repmat(BlackPath, length(Paths{end}), 1);
FsRgb(Paths{end},:) = repmat(YellowPath, length(Paths{end}), 1);
CfRgb(Paths{end},:) = repmat(WhitePath, length(Paths{end}), 1);

% Redraw paths
RedrawAllPaths;


%% Display paths on all maps
function RedrawAllPaths
global PolRgb EccRgb FsRgb CfRgb pp pe pf pc
set(pp, 'FaceVertexCData', PolRgb); 
set(pe, 'FaceVertexCData', EccRgb); 
set(pf, 'FaceVertexCData', FsRgb); 
set(pc, 'FaceVertexCData', CfRgb); 


%% Fill in the ROIs assigned to the list bounded by the paths & save 
function OutRoi = FloodFillRoi(RoiNum, RoiSeed)
global SrfName RoiList RoiColours Srf Vertices Paths CfRgb hc pc 

% Only fill if a vertex has been assigned
WhitePath = [1 1 1];
Cmap = RoiColours;
% Labelled vertices
Fvs = zeros(1,size(Vertices,1)); 
% Switch to cortex & ROI figure
figure(hc);
OutRoi = RoiNum;

% Flood filling
fc = 0;
Bvx = RoiSeed; % Border vertex indeces (initialised as selected vertex)
Fvs(RoiSeed) = 1; % Filled vertices (initialised as selected vertex)

% While there are border vertices
while ~isempty(Bvx)        
    % Label all vertices on surrounding path
    Cvx = find(Fvs == 1); % Vertices on the border
    Bvx = []; % New border vertex indeces
    nb = samsrf_neighbours(Cvx, Srf.Faces); % Immediate neighbour vertices
    bn = (CfRgb(nb,1) == 1 & CfRgb(nb,2) == 1 & CfRgb(nb,3) == 1) | Fvs(nb)' == 2; % Neighbour vertices on paths or previously labelled ones
    if sum(bn) > 0
        Fvs(Cvx) = 2; % Neighbours on paths so label current vertex as already labelled
        Fvs(nb(bn)) = 2; % Label paths or neighbours to paths
    end
    nb(bn) = []; % Remove those on paths or neighbouring paths
    if sum(Fvs(nb) > 0) < length(nb)
        Bvx = nb;
    end
    Fvs(Fvs == 1) = 2; % Label previously labelled vertices
    Fvs(Bvx) = 1; % Label new vertices

    % Count filling steps
    fc = fc + 1;
    if fc == 100
        disp('Warning: Flood filling is taking unusually long...');
    end    
    if fc > 200
        OutRoi = -1;
        break
    end
end

% Only if filling went okay
if OutRoi > 0
    % Paint the map
    Rgb = (repmat(Cmap(RoiNum,:),sum(Fvs>0),1) + CfRgb(Fvs>0,:)) / 2; % What colour?
    CfRgb(Fvs>0,:) = Rgb; % Colour the ROI
    % Redraw paths
    for i = 1:length(Paths)
        CfRgb(Paths{i},:) = repmat(WhitePath, length(Paths{i}), 1);
    end
    % Update map
    set(pc, 'FaceVertexCData', CfRgb); 

    % Save new ROI
    if ~exist(['ROIs_' SrfName(4:end)], 'dir') 
        mkdir(['ROIs_' SrfName(4:end)]);
    end
    samsrf_srf2label(Srf, ['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_' RoiList{RoiNum}], 1, find(Fvs));
end


%% Clone all cameras after a rotation occurred
function postrotcallback(obj,evd)
global Points SelectedPath hp he hf hc
Points = [];
SelectedPath = 0;
ClearDataCursors;
CalculateMapRgb(false);
% Get current view
c = get(gca, 'CameraPosition');
% Loop through figures
for h = [hp he hf hc]
    figure(h);
    set(gca, 'CameraPosition', c);
end


%% What happens when a vertex is selected
function txt = getvertexfun(empt, event_obj)
global Vertices Points CfRgb pc
pos = get(event_obj, 'Position');
v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1);
txt = ''; % No data info 
if isempty(Points)
    % Add first point
    Points = v;
    CfRgb(v,:) = [1 0 0];
    set(pc, 'FaceVertexCData', CfRgb); 
    disp(['Selected new vertex #' num2str(v)]);
else
    % Check that new point is not too far away
    d = sqrt((Vertices(v,1)-Vertices(Points(end),1)).^2 + (Vertices(v,2)-Vertices(Points(end),2)).^2 + (Vertices(v,3)-Vertices(Points(end),3)).^2);
    if d < 40
        Points = [Points; v];
        CfRgb(v,:) = [1 0 0];
        set(pc, 'FaceVertexCData', CfRgb); 
    else
        disp(['Too far from previous vertex! ' num2str(d) ' mm']);
    end
end


%% Clear global variables when figure is closed
function closereq(src, evnt)
clear global 
delete(src);
close all


%% --- Outputs from this function are returned to the command line.
function varargout = DelineationTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Add new path
function pushbutton1_Callback(hObject, eventdata, handles)
% Update message
set(handles.text2, 'String', 'Busy...');
pause(0.001);
% Connect & draw the paths
DrawPaths; 
% Update message
set(handles.text2, 'String', 'Added new path.');


%% Clear last path
function pushbutton2_Callback(hObject, eventdata, handles)
global Points Paths PolRgb EccRgb FsRgb CfRgb pp pe pf pc

% Update message
set(handles.text2, 'String', 'Busy...');
pause(0.001);

% Clear previous path & all waypoints
Points = [];
if ~isempty(Paths)
    if length(Paths) == 1
        Paths = {}; % Only one path currently exists so delete
    else
        Paths = Paths(1:end-1); % Delete most recent path
    end
end
CalculateMapRgb(false);
 
% Redraw paths
RedrawAllPaths;

% Update message
set(handles.text2, 'String', 'Cleared last path.');


%% Clear everything
function pushbutton3_Callback(hObject, eventdata, handles)
global Points Paths RoiSeeds PolRgb EccRgb FsRgb CfRgb pp pe pf pc

% Update message
set(handles.text2, 'String', 'Busy...');
pause(0.001);

% Clear paths & waypoints
Points = [];
Paths = {};
ClearDataCursors;

% Clear ROI seed vertices
RoiSeeds(:) = NaN;

% Recalculate all the maps colours & don't load delineations!
CalculateMapRgb(false); 

% Display paths on map
set(pp, 'FaceVertexCData', PolRgb); 
set(pe, 'FaceVertexCData', EccRgb); 
set(pf, 'FaceVertexCData', FsRgb); 
set(pc, 'FaceVertexCData', CfRgb); 

% Update message
set(handles.text2, 'String', 'Cleared everything!');


%% (Un-)Label ROI
function pushbutton4_Callback(hObject, eventdata, handles)
global SrfName Srf Points RoiList RoiSeeds

% Add already assigned vertices to ROI list
RoiSelList = RoiList;
for i = 1:length(RoiList)
    if ~isnan(RoiSeeds(i))
        RoiSelList{i} = [RoiList{i} ' ~ #' num2str(RoiSeeds(i))];
    end
end

% Select new ROI
RoiNum = listdlg('ListString', RoiSelList, 'SelectionMode', 'single', 'Name', 'Which ROI?');
if ~isempty(RoiNum) 
    if isempty(Points)
        % Remove assigned vertex
        if ~isnan(RoiSeeds(RoiNum))
            RoiSeeds(RoiNum) = NaN;
            CalculateMapRgb(false); 
            % Update message
            set(handles.text2, 'String', ['Unlabelled' RoiList{RoiNum}]);
            % Delete unlabelled ROI!
            delete(['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_' RoiList{RoiNum} '.label']);
        end
    else
        % Which vertex is seed?
        RoiSeeds(RoiNum) = Points(end);
        % Update message
        set(handles.text2, 'String', 'Filling ROI...');
        pause(0.001);
        % Fill & save ROI
        Roi = FloodFillRoi(RoiNum, RoiSeeds(RoiNum));
        % Update message
        if Roi == -1
            set(handles.text2, 'String', 'Something''s amiss!');
            disp('Something is amiss: Flood filling took too way long!');
            RoiSeeds(RoiNum) = NaN;
        else
            set(handles.text2, 'String', ['Labelled ' RoiList{RoiNum}]);
        end
    end
end

% Clear way points etc
Points = [];
ClearDataCursors;


%% Save delineations
function pushbutton7_Callback(hObject, eventdata, handles)
global SrfName Paths RoiList RoiSeeds Srf
% Saves paths & ROI seeds to disc
Vs = [];
for i = 1:length(Paths)
    Vs = [Vs; Paths{i}];
end
save(['del_' SrfName], 'SrfName', 'Vs', 'Paths', 'RoiList', 'RoiSeeds');
% Update message
set(handles.text2, 'String', ['Saved del_' SrfName]);

% Combine ventral & dorsal quadrants of V2 and V3 if they exist 
for a = 2:3
    IsOkay = samsrf_combine_labels({['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_V' num2str(a) 'v'] ...
                                    ['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_V' num2str(a) 'd']}, ...
                                    ['ROIs_' SrfName(4:end) filesep Srf.Hemisphere '_V' num2str(a)]);
    if IsOkay    
        % Update message
        set(handles.text2, 'String', 'Combined Quadrants');
    else
        disp(['Couldn''t combine V' num2str(a) ' quadrants!']);
    end
end


%% Create smart path
function pushbutton9_Callback(hObject, eventdata, handles)
global Srf Vertices Points Paths R2Thrsh EccThrsh PolRgb EccRgb FsRgb CfRgb

% Colours of paths
WhitePath = [1 1 1];
BlackPath = [0 0 0];
YellowPath = [1 1 0]; 

if ~isempty(Points)
    % Loop through seed points
    for p = 1:length(Points)
        % Find vertices around seed
        Vs = samsrf_georoi(Points(p), 10, Vertices, Srf.Faces);

        % Field sign map
        Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
        Rho(Srf.Data(1,:) <= R2Thrsh) = NaN;
        Rho(Rho < EccThrsh(1)) = NaN;
        Rho = Rho / EccThrsh(2);
        Rho(Rho > 1) = NaN; 
        Fs = sign(Srf.Data(strcmpi(Srf.Values, 'Field Sign'), :));
        Fs(Srf.Data(1,:) <= R2Thrsh | isnan(Rho)) = 0;

        % Search through local neighbourhood for reversals
        NewPath = [];
        for v = 1:length(Vs)
            % Neighbours of current vertex
            Nv = samsrf_neighbours(Vs(v), Srf.Faces);
            % If any neighbours have another field sign
            if find(Fs(Nv) ~= Fs(Vs(v)))
                NewPath = [NewPath; Vs(v)]; % Add to path
            end
        end

        % Add new path
        Paths{length(Paths) + 1} = NewPath;

        % Colour vertices on path
        PolRgb(Paths{end},:) = repmat(WhitePath, length(Paths{end}), 1);
        EccRgb(Paths{end},:) = repmat(BlackPath, length(Paths{end}), 1);
        FsRgb(Paths{end},:) = repmat(YellowPath, length(Paths{end}), 1);
        CfRgb(Paths{end},:) = repmat(WhitePath, length(Paths{end}), 1);

        % Redraw paths
        RedrawAllPaths;
    end
    
    % Clean up
    Points = [];
    ClearDataCursors;
end

%% Swap between activation & field sign map
function pushbutton10_Callback(hObject, eventdata, handles)
global Srf hf IsFsMap ActPrct

if sum(strcmpi(Srf.Values, 'Field Sign')) == 0 
    IsFsMap = false;
end
if IsFsMap || sum(strcmpi(Srf.Values, 'Field Sign')) == 0 
    % Load activation map
    [fn pn] = uigetfile([Srf.Hemisphere '_*.mat']);
    if fn ~= 0
        ActMap = load([pn filesep fn]);
        ActMap.Srf = samsrf_expand_srf(ActMap.Srf);
        if isfield(ActMap.Srf, 'Values')
            Values = ActMap.Srf.Values;
        else
            Values = {};
            for i = 1:size(ActMap.Srf.Data,1)
                Values{i} = ['Volume #' num2str(i)];
            end
        end
        if size(ActMap.Srf.Data,1) == 1
            dt = 1;
        else
            dt = listdlg('ListString', Values, 'SelectionMode', 'single');
        end
        if ~isempty(dt)
            Srf.Data(end,:) = ActMap.Srf.Data(dt,:);
            set(hf, 'name', [Values{dt} ' map (' num2str(ActPrct(1)) '-' num2str(ActPrct(2)) ' %ile)']);
            IsFsMap = false;
            if sum(strcmpi(Srf.Values, 'Field Sign')) == 1 
                set(handles.pushbutton10, 'String', 'Load Field Sign');
            end
        end
    end
else
    % Load field sign map
    set(hf, 'name', 'Field sign map');
    IsFsMap = true;
    set(handles.pushbutton10, 'String', 'Load Activation');
end

% Redraw all
CalculateMapRgb(false);
RedrawAllPaths;


%% Selects the path at the current vertex
function pushbutton12_Callback(hObject, eventdata, handles)
global Points Paths PolRgb EccRgb FsRgb CfRgb SelectedPath

% Colours of paths
WhitePath = [1 1 1];
BlackPath = [0 0 0];
YellowPath = [1 1 0]; 

if ~isempty(Points) && ~isempty(Paths) 
    % Find selected path
    SelectedPath = 0;
    for i = 1:length(Paths)
        if find(Paths{i} == Points(end))
            SelectedPath = i;
        end
    end
    
    Points = [];
    ClearDataCursors;

    % Colour vertices on path
    PolRgb(Paths{SelectedPath},:) = repmat(YellowPath, length(Paths{SelectedPath}), 1);
    EccRgb(Paths{SelectedPath},:) = repmat(WhitePath, length(Paths{SelectedPath}), 1);
    FsRgb(Paths{SelectedPath},:) = repmat(BlackPath, length(Paths{SelectedPath}), 1);
    CfRgb(Paths{SelectedPath},:) = repmat(YellowPath, length(Paths{SelectedPath}), 1);

    % Redraw paths
    RedrawAllPaths;

    % Update message
    set(handles.text2, 'String', 'Selected this path.');
else
    % Update message
    set(handles.text2, 'String', 'No path to select here!');
end

%% Delete the selected path
function pushbutton13_Callback(hObject, eventdata, handles)
global Points Paths SelectedPath

Points = [];
ClearDataCursors;
if SelectedPath > 0
    if ~isempty(Paths)
        if length(Paths) == 1
            Paths = {}; % Only one path currently exists so delete
        else
            px = 1:length(Paths); % Path indeces
            px(SelectedPath) = []; % Remove selected path index
            Paths = Paths(px); % Selected path has been removed
        end
    end
    CalculateMapRgb(false);

    % Redraw paths
    RedrawAllPaths;

    % Update message
    set(handles.text2, 'String', 'Cleared this path.');
else
    % Update message
    set(handles.text2, 'String', 'No path selected!');
end
SelectedPath = 0;


%% R^2 threshold 
function edit1_Callback(hObject, eventdata, handles)
global R2Thrsh

R2Thrsh = str2double(get(hObject,'String'));
CalculateMapRgb(false);
% Redraw paths
RedrawAllPaths;

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Eccentricity range
function edit2_Callback(hObject, eventdata, handles)
global EccThrsh

EccStr = get(hObject,'String');
EccVec = eval(['[' EccStr ']']);
if isempty(EccVec)
    EccThrsh = [0 90];
elseif length(EccVec) == 1
    EccThrsh = [0 EccVec];
else
    EccThrsh = EccVec(1:2);
end
set(hObject, 'String', num2str(EccThrsh));
CalculateMapRgb(false);
% Redraw paths
RedrawAllPaths;


function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
