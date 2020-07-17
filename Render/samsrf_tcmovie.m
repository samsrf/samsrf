function samsrf_tcmovie(Srf, Mesh, Thrsh, Paths, CamView, MapType, ColorMap)
% samsrf_tcmovie(Srf, Mesh, Thrsh, Paths, CamView, MapType, ColorMap)
%
% Displays a cortical mesh overlayed with a timecourse of the either 
% a) the measured signal or b) the model fit timecourse. Each timepoint is 
% rendered, and captured to form a movie saved to disk.
%
% MapType = 'signal' or 'model'
%
% 17/07/2020 - SamSrf 7 version (DSS)
%

%% Defaults
if nargin < 3
    Thrsh = [];
end
if nargin < 4
    Paths = [];
end
if nargin < 5
    CamView = [];
end
if nargin < 6
    MapType = [];
end
if nargin < 7
    ColorMap = [];
end

%% Load default parameters?
load('SamSrf_defaults.mat');

%% Create figure handle
fh = figure;

%% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

%% Load curvature
Curv = Srf.Curvature;
if exist('def_curv', 'var')
    if strcmpi(def_curv, 'Greyscale')
        % Greyscale for curvature
        CurvGrey = gray(11); % Grey scale colour map
        CurvGrey = CurvGrey(1:10,:); % Remove black & white
    elseif strcmpi(def_curv, 'FreeSurfer')
        % Freesurfer-like curvature
        CurvGrey = gray(4); % Gray scale colour map
        CurvGrey = CurvGrey(2:3, :); % Remove black & white
    else
        % Black & White curvature
        CurvGrey = gray(2); % Black & white colour map
    end
else
    % Greyscale for curvature
    CurvGrey = gray(11); % Grey scale colour map
    CurvGrey = CurvGrey(1:10,:); % Remove black & white
end

% Transform curvature
Curv(isnan(Curv)) = 0;
Curv = -Curv + 0.5;
Curv(Curv <= 0) = 0.000001;
Curv(Curv > 1) = 1;
Curv = ceil(Curv * size(CurvGrey,1));

%% Load mesh
Mesh = lower(Mesh);
Mesh(1) = upper(Mesh(1));
if isfield(Srf, Mesh)
    Vertices = Srf.(Mesh);
else
    Vertices = Srf.Vertices;
end
Vertices(isnan(Srf.Vertices(:,1)),:) = NaN; 
Faces = Srf.Faces;
if length(Curv) == 1
    Curv = repmat(Curv, 1, size(Vertices,1));
end

%% Select data type

% Input
if isempty(MapType)
    Values = {'Signal', 'Model'};
    dt = listdlg('ListString', Values, 'SelectionMode', 'single');
    if isempty(dt)
        return
    end
    MapType = Values{dt};
end

% Pick data out
if strcmpi(MapType, 'Signal')
    X = Srf.Y;
elseif strcmpi(MapType, 'Model')
    X = Srf.X;
end

% Add a dummy timepoint with baseline signal
X = [zeros(1, size(X, 2)); X];

%% Select paths if desired
if isempty(Paths)
    pf = dir('*path*');
    lf = dir('*.label');
    df = dir('del_*.mat');
    pf = [pf; lf; df];
    pf = {pf.name}';
    if ~isempty(pf)
        ps = listdlg('ListString', pf);
        if isempty(ps)
            Paths = {};
        else
            Paths = {pf{ps}}';
        end
    else 
        Paths = {};
    end
end
if ~isempty(Paths) 
    Vs_paths = [];
    for i = 1:length(Paths)
        if strfind(Paths{i}, '.label')
            % Load label
            Label = samsrf_loadlabel(Paths{i}(1:end-6));
            Vs_paths = [Vs_paths; samsrf_borderpath(Srf, Label)];
        else
            % Load paths
            Vs_paths = [Vs_paths; samsrf_loadpaths(Paths{i})];
        end
    end
else
    % No paths needed
    Vs_paths = [];
end

% Colour
PathColour = [1 1 1];

%% Thresholding

% R^2 threshold
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
    r = Srf.Data(1,:) <= 0.01 | isnan(Srf.Data(1,:));
else 
    r = isnan(Srf.Data(1,:));
end

% Input
if isempty(Thrsh)
    pct = prctile(Srf.Data(1, ~r), [2.5 97.5]);
    Thrsh = inputdlg(['R^2 threshold?  95%CI: ' num2str(pct(1)) ', ' num2str(pct(2))]);
    Thrsh = eval(['[' cell2mat(Thrsh) ']']);
    if isempty(Thrsh)
        Thrsh = pct(1);
    end
end

%% Calculate colours for the map 

% Scaling
Pha = round(X / abs(max(X(:))) * 100) + 100;
Pha(Pha==0) = 1;

% Color map
if isempty(ColorMap)
    
    % Hot/cold colour map
    Cmap = [flipud(colormap(winter(100))); colormap(hot(100)); CurvGrey];
else
    Cmap = [ColorMap; CurvGrey];
end

% Loop volumes
for i = 1:size(X, 1)
    Idx = Pha(i,:);
    Idx(:, r|isnan(Idx) | Idx == 100) = 200 + Curv(r | isnan(Idx) | Idx == 100);
    Idx(Idx < 0) = 1;
    Colours(:,:,i) = Cmap(Idx, :);
end

%% Draw paths
Colours(Vs_paths,:) = repmat(PathColour, length(Vs_paths), 1);

%% Initial mesh display

% Title
set(fh, 'name', [MapType ' ( R^2 > ' num2str(Thrsh) ' )'], 'color', 'w');

% Show first frame
patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours(:,:,1), 'FaceColor', 'interp', 'EdgeColor', 'none');

% Tidy
axis off;
daspect([1 1 1]);
ax = gca;
ax.Clipping = 'off';
set(gca, 'projection', 'perspective');

% Apply CamView, if any
if ~isempty(CamView)
    set(gca, 'view', CamView(1:2));
    zoom(CamView(3));
end

%% User input

% Ask to continue
Ok = input('Make movie with this viewpoint? (y/n): ', 's');
if strcmpi(Ok, 'y')
    disp('Capturing...')
else
    disp('Stopping now')
    return;
end

%% Shoot frames

% Make figure invisiible to save on processor
fh.Visible = 'off';

% Open video file
Filename = [Srf.Hemisphere '_' MapType];
while exist([Filename '.avi'], 'file')
    Filename = [Filename '_1'];
end
vidObj = VideoWriter([Filename '.avi']);
vidObj.FrameRate = 5;
open(vidObj);

% Loop timepoints
for i = 1:size(X, 1)

    % Draw frame
    patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours(:,:,1), 'FaceColor', 'interp', 'EdgeColor', 'none');
    patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours(:,:,i), 'FaceColor', 'interp', 'EdgeColor', 'none');

    % Capture
    CurrFrame = getframe(fh);
    
    % Write to file
    writeVideo(vidObj, CurrFrame);
end

% Close video file
close(vidObj);

% Show last frame
fh.Visible = 'on';

% Done
%
