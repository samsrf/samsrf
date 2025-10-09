function varargout = DisplayMaps(varargin)
% 
% ******* SamSrf pRF map display tool ******* 
%
% Run this GUI to display maps and make figures from them.
%
%  If you wish to change default ROI to limit the mesh, 
%  change them below where indicated by '%%% %%% %%%'.
%
% 20/04/2022 - SamSrf 8 version (DSS)

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

SamSrfDefs = LoadSamSrfDefaults;
if ~isfield(SamSrfDefs, 'def_disproi')
    SamSrfDefs.def_disproi = NaN; 
end

%% Default parameters (Feel free to edit)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
RoiName = SamSrfDefs.def_disproi; % ROI name to limit the mesh (without hemisphere prefix, e.g. 'occ') - If this is NaN & Srf.Roi exists, it uses this instead
Pval = 0.0001; % Starting p-value with which to threshold maps
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

%% Version info
[vn,vd] = samsrf_version;

%% Welcome message
samsrf_clrscr; 
samsrf_disp('*** SamSrf Map Display Tool ***');
samsrf_newline;

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
if isscalar(fn)
    samsrf_error('No map file selected!');
end
cd(pn); % Navigate to folder as otherwise anatomicals cannot be loaded
load(fn); % Load Srf 
samsrf_newline;
samsrf_disp(['Loading surface data file: ' pn fn]);
samsrf_newline;
% Is there a pRF model fit?
if exist('Model', 'var')
    Srf.Model = Model;
    clear Model
end
% Backwards compatibility
if ~isfield(Srf, 'Version')
    % Compatible with some(!) files before version field was introduced
    Srf.Version = -Inf; 
end
if Srf.Version < 6
    Srf.Model = struct; % Empty structure
end

% Limit multi-subject Srf?
if size(Srf.Data,3) > 1
    samsrf_disp(['Srf.Data contains ' num2str(size(Srf.Data,3)) ' subjects.']);
    if ~isfield(Srf, 'Subjects')
        Srf.Subjects = {};
        for i = 1:size(Srf.Data,3)
            Srf.Subjects{i} = ['Subject #' num2str(i)];
        end
    end
    SubjNum = listdlg('ListString', Srf.Subjects, 'SelectionMode', 'single', 'PromptString', 'Which subject?');
    if isempty(SubjNum)
        samsrf_disp('You must pick a subject! Selecting #1...');
        SubjNum = 1;
    end
    Srf.Data = Srf.Data(:,:,SubjNum);
    % Limit raw data equally
    if isfield(Srf, Raw_Data)
        Srf.Raw_Data = Srf.Raw_Data(:,:,SubjNum);
    end
end
% In case only raw data has multiple subjects
if isfield(Srf, 'Raw_Data') && size(Srf.Raw_Data,3) > 1
    samsrf_disp(['Srf.Raw_Data contains ' num2str(size(Srf.Raw_Data,3)) ' subjects.']);
    if ~isfield(Srf, 'Subjects')
        Srf.Subjects = {};
        for i = 1:size(Srf.Raw_Data,3)
            Srf.Subjects{i} = ['Subject #' num2str(i)];
        end
    end
    SubjNum = listdlg('ListString', Srf.Subjects, 'SelectionMode', 'single', 'PromptString', 'Which subject?');
    if isempty(SubjNum)
        samsrf_disp('You must pick a subject! Selecting #1...');
        SubjNum = 1;
    end
    Srf.Raw_Data = Srf.Raw_Data(:,:,SubjNum);
end

% Expand if necessary
[Srf,vx] = samsrf_expand_srf(Srf);
% If bilateral ask about hemisphere
if isfield(Srf, 'Nvert_Lhem')
    qh = questdlg('Which hemisphere?', 'Bilateral Srf', 'Left', 'Right', 'Both', 'Both');   
    if isempty(qh) 
        qh = 'Both';
    end 
    if qh(1) == 'L'
        Srf = samsrf_hemi_srfs(Srf); % Remove right hemisphere 
    elseif qh(1) == 'R'
        [~,Srf] = samsrf_hemi_srfs(Srf); % Remove left hemisphere 
    end
end
% What ROI to samsrf_display
if isnan(RoiName)
    % Use ROI from Srf
    samsrf_disp('Using ROI from Srf if it exists');
elseif isempty(RoiName)
    % No ROI defined
    samsrf_disp('No default ROI defined!');
    vx = NaN;
elseif ischar(RoiName) && RoiName(1) == '<' || RoiName(1) == '>' 
    % If ROI defined by coordinates
    if length(RoiName) == 1
        samsrf_error('You must define inflated mesh coordinate in SamSrfDefs.defs_disproi!');
    end
    switch RoiName(2)
        case 'X'
            wv = Srf.Inflated(:,1);
        case 'Y'
            wv = Srf.Inflated(:,2);
        case 'Z'
            wv = Srf.Inflated(:,3);
        otherwise
            samsrf_error('Invalid inflated mesh coordinate specified in SamSrfDefs.defs_disproi!');
    end
    if length(RoiName) < 3
        samsrf_error('You must define inflated mesh cut-off coordinate in SamSrfDefs.defs_disproi!');
    end
    wc = str2double(RoiName(3:end));
    if RoiName(1) == '<'
        vx = find(wv < wc);
    elseif RoiName(1) == '>'
        vx = find(wv > wc);
    end
    samsrf_disp(['Only samsrf_displaying inflated mesh vertices with ' RoiName]);
else
    % If ROI was named, load that instead
    vx = samsrf_loadlabel([pn filesep Srf.Hemisphere '_' RoiName]);
    samsrf_disp(['Only samsrf_displaying ROI ' Srf.Hemisphere '_' RoiName]);   
end

% Add transparent border
if ~isnan(vx)
    Srf.VxAlpha = zeros(size(Srf.Vertices,1),1);
    Srf.VxAlpha(vx) = 1;
    % Loop thru steps
    for i = 1:19
        nb = samsrf_borderpath(Srf, vx); % Vertices surrounding the ROI
        Srf.VxAlpha(nb) = 1 - i/20; % Set alpha for this border
        vx = [vx; nb]; % Add border to ROI
    end
end

% Remove non-ROI vertices if desired
if ~isnan(vx)
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

%% Is this raw data with noise ceiling?
if isfield(Srf, 'Noise_Ceiling')
    Srf.Values(2:end+1) = Srf.Values;
    Srf.Values{1} = 'Noise Ceiling';
    Srf.Data = [Srf.Noise_Ceiling; Srf.Data];
    Srf = rmfield(Srf, 'Noise_Ceiling');
end

%% Update map list
Values = Srf.Values;
% Does file contain retinotopic maps?
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
end
% Set to first item in the list
set(handles.popupmenu2, 'String', Values, 'Value', 1); 

%% Update pRF inspector list
PrfInspList = {'Vertex inspector off'};
if isfield(Srf, 'ConFlds') % Can samsrf_display connective field profiles
    PrfInspList{2} = 'Connective field profiles';
    Tmp = load(Srf.Model.Template); % Load template map
    Tmp = samsrf_expand_srf(Tmp.Srf); % Expand template map
    Tmp = Tmp.Data(2:3,Srf.SeedVx); % Restrict to template seed ROI 
    Srf.TempMap = Tmp;
    PrfInspList{3} = 'Visual CF profiles'; % Can samsrf_display visual CF profiles
end
if isfield(Srf, 'Rmaps') % Can samsrf_display reverse correlation profiles
    PrfInspList{2} = 'Reverse correlation profiles';
end
if isfield(Srf, 'Model') && isfield(Srf.Model, 'Prf_Function') % Only if a pRF model fit
    PrfInspList{end+1} = 'pRF parameter estimates';
end
if isfield(Srf, 'Y') && (isfield(Srf, 'X') || isfield(Srf, 'X_glm')) % Can show both observed & predicted time series
    PrfInspList{end+1} = 'Observed vs predicted time series';
end
if ~isfield(Srf, 'Y') && ~isfield(Srf, 'X') % Raw data file so observed time series only
    if contains(Srf.Functional, 'GLM contrasts')
        PrfInspList{end+1} = 'GLM contrasts';
    else
        PrfInspList{end+1} = 'Observed time series only';
    end
end
% Update & set to first item in the list
set(handles.popupmenu3, 'String', PrfInspList, 'Value', 1); 

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

% pRF inspector mode? 
contents = cellstr(get(handles.popupmenu3,'String'));
PrfInsp = get(handles.popupmenu3,'Value');
PrfInsp = contents{PrfInsp};
if strcmpi(PrfInsp, 'pRF parameter estimates')
    % samsrf_display model fit pRF
    if isfield(Srf, 'Model_')
        Srf.Model = Srf.Model_;
        Srf = rmfield(Srf, 'Model_');
    end
    % Don't samsrf_display reverse correlation profile!
    if isfield(Srf, 'Rmaps')
        Srf.Rmaps_ = Srf.Rmaps;
        Srf = rmfield(Srf, 'Rmaps');
    end
elseif strcmpi(PrfInsp, 'Reverse correlation profiles')
    % samsrf_display reverse correlation profile
    if isfield(Srf, 'Rmaps_')
        Srf.Rmaps = Srf.Rmaps_;
        Srf = rmfield(Srf, 'Rmaps_');
    end
elseif strcmpi(PrfInsp, 'Visual CF profiles')
    % samsrf_display connective field profile
    if isfield(Srf, 'TempMap_')
        Srf.TempMap = Srf.TempMap_;
        Srf = rmfield(Srf, 'TempMap_');
    end
    if isfield(Srf, 'ConFlds')
        Srf.ConFlds_ = Srf.ConFlds;
        Srf = rmfield(Srf, 'ConFlds');
    end
elseif strcmpi(PrfInsp, 'Connective field profiles')
    % samsrf_display connective field profile
    if isfield(Srf, 'ConFlds_')
        Srf.ConFlds = Srf.ConFlds_;
        Srf = rmfield(Srf, 'ConFlds_');
    end
    if isfield(Srf, 'TempMap')
        Srf.TempMap_ = Srf.TempMap;
        Srf = rmfield(Srf, 'TempMap');
    end
elseif strcmpi(PrfInsp, 'Observed vs predicted time series')
    % Plot observed vs predicted time series
    if isfield(Srf, 'Y_')
        Srf.Y = Srf.Y_;
        Srf = rmfield(Srf, 'Y_');
    end
    % Don't samsrf_display model fit pRF
    if isfield(Srf, 'Model')
        Srf.Model_ = Srf.Model;
        Srf = rmfield(Srf, 'Model');
    end
elseif strcmpi(PrfInsp, 'Observed time series only') || strcmpi(PrfInsp, 'GLM contrasts')
    % Plot observed time series
    if isfield(Srf, 'Y')
        Srf.Y_ = Srf.Y;
        Srf = rmfield(Srf, 'Y');
    end
    if isfield(Srf, 'X') && isempty(Srf.X)
        Srf = rmfield(Srf, 'X');
    end
else
    % Not samsrf_displaying pRFs
    if isfield(Srf, 'Model')
        Srf.Model_ = Srf.Model;
        Srf = rmfield(Srf, 'Model');
    end
    if isfield(Srf, 'Rmaps')
        Srf.Rmaps_ = Srf.Rmaps;
        Srf = rmfield(Srf, 'Rmaps');
    end
    if isfield(Srf, 'TempMap')
        Srf.TempMap_ = Srf.TempMap;
        Srf = rmfield(Srf, 'TempMap');
    end
    if isfield(Srf, 'ConFlds')
        Srf.ConFlds_ = Srf.ConFlds;
        Srf = rmfield(Srf, 'ConFlds');
    end
    if ~isfield(Srf, 'X')
        Srf.X = [];
    end
    if isfield(Srf, 'Y')
        Srf.Y_ = Srf.Y;
        Srf = rmfield(Srf, 'Y');
    end
end    

% What R^2 threshold?
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
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
EccRng = eval(['[' get(handles.edit4, 'String') ']']);
if isempty(EccRng)
    EccRng = [0 Inf];
elseif length(EccRng) == 1
    EccRng(2) = Inf;
end
EccRng = EccRng(1:2);
ThrshVec = [ThrshVec EccRng];

% Map transparency level
ThrshVec = [ThrshVec eval(['[' get(handles.edit6, 'String') ']'])];

% Camera angle
if IsNewFig
    FigHdl = figure;
    if ~exist('SamSrfDefs.defs_views', 'var')
        % Focus on early visual cortex
        if Srf.Hemisphere(1) == 'l'
            % Left hemisphere
            CamView = [36 -20 1.8];
        elseif Srf.Hemisphere(1) == 'r'
            % Right hemisphere
            CamView = [-35 -35 2.2];
        else
            % Both hemispheres
            CamView = [4 -30 2.2];
        end
    else
        % Use default camera angle
        if Srf.Hemisphere(1) == 'l'
            % Left hemisphere
            CamView = SamSrfDefs.defs_views(:,1)';
        elseif Srf.Hemisphere(1) == 'r'
            % Right hemisphere
            CamView = SamSrfDefs.defs_views(:,2)';
        else
            % Both hemispheres
            if size(SamSrfDefs.defs_views,2) > 2
                CamView = SamSrfDefs.defs_views(:,3)';
            else
                CamView = [4 -30 2.2];
            end
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
samsrf_newline;
samsrf_disp('SamSrf says:'); 
samsrf_disp(' "I hope you''ll come back to be annoyed by me again soon...');
samsrf_disp('  Mā te wā!"');
samsrf_newline;

clear global Srf RoiName Pval Paths FigHdl R2Thrsh PatchHdl fv
delete(src);


%% Draw ROIs or paths
function pushbutton1_Callback(hObject, eventdata, handles)
global Paths

% Load paths
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
% Draw paths
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
Values = Srf.Values;
% Does file contain retinotopic maps?
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0')
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
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


%% Which map to samsrf_display?
function popupmenu2_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

global Srf

% Which map? 
contents = cellstr(get(handles.popupmenu2,'String'));
MapType = contents{get(handles.popupmenu2,'Value')};
MapNum = find(strcmpi(Srf.Values, MapType));

% Default cut-offs
if strcmpi(MapType, 'R^2') || strcmpi(MapType, 'nR^2') || strcmpi(MapType, 'Noise Ceiling')
    % Goodness-of-fit
    set(handles.edit2, 'String', '0 1');

elseif strcmpi(MapType, 'Polar') || strcmpi(MapType, 'Phase')
    % Angles & phases
    set(handles.edit2, 'String', '0 0');

elseif strcmpi(MapType, 'Eccentricity') 
    % Eccentricity
    MapData = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    set(handles.edit2, 'String', ['0 ' num2str(prctile(MapData(MapData > 0), 95))]);

elseif strcmpi(MapType, 'x0') || strcmpi(MapType, 'y0') || strcmpi(MapType, 'Sigma') || strcmpi(MapType, 'Fwhm') || strcmpi(MapType, 'nSigma') ...
       || strcmpi(MapType, 'Centre') || strcmpi(MapType, 'Surround') || strcmpi(MapType, 'Sigma1') || strcmpi(MapType, 'Sigma2')  || strcmpi(MapType, 'Spread') ...
       || strcmpi(MapType, 'Suppression') || strcmpi(MapType, 'Cmf') || strcmpi(MapType, 'Visual Area') || strcmpi(MapType, 'Surface Area') 
    % Measures that should be positive
    MapData = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    set(handles.edit2, 'String', ['0 ' num2str(prctile(MapData(MapData > 0), 75))]) 

elseif strcmpi(MapType, 'Field Sign')
    % Field sign
    set(handles.edit2, 'String', '0 3'); % Ranges -1 to +1 but this makes it prettier
    
elseif strcmpi(MapType, 'ROI') || strcmpi(MapType, 'ROIs')
    % ROI numbers
    set(handles.edit2, 'String', ['0 ' num2str(nanmax(Srf.Data(MapNum,:)))]);

else
    % Anything else
    set(handles.edit2, 'String', ['0 ' num2str(prctile(Srf.Data(MapNum,:), 99))]);
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


%% Toggles transparent or old-style path colours
function togglebutton3_Callback(hObject, eventdata, handles)
global Paths

% Change button string
if get(hObject, 'Value')
    % Turn on transparent paths
    set(hObject, 'String', 'Old-Style Paths'); % Button is now for going to old-style paths
    if isnan(Paths{end})
        Paths = Paths(1:end-1); % Remove NaN from end if it exists
    end
else
    % Turn on old-style paths
    set(hObject, 'String', 'Transparent Paths'); % Button is now for going back to transparent paths
    Paths{end+1} = NaN;
end
RedrawMaps(handles, false);


%% Map transparency 
function edit6_Callback(hObject, eventdata, handles)

RedrawMaps(handles, false);

%%
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Vertex inspector menu
function popupmenu3_Callback(hObject, eventdata, handles)
global fv 

PrfInsp = get(hObject,'Value');
if PrfInsp > 1
    try 
        figure(fv);
    catch
        fv = figure;
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [.5 .3 .3 .3]);
    end
else
    close(fv);
end
RedrawMaps(handles, false);

%%
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
