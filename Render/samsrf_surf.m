function PatchHandle = samsrf_surf(Srf, Mesh, Thrsh, Paths, CamView, MapType, PatchHandle)
%
% PatchHandle = samsrf_surf(Srf, Mesh, [Thrsh, Paths, CamView, MapType, PatchHandle])
%
% Displays a cortical mesh (e.g. 'inflated') overlaid with the pRF 
%  statistics from Srf. A GUI is used to select the data type to be displayed. 
%  If the surface mesh can't be located in the Srf, the function simply takes 
%  the WM-GM surface which should already be in the Srf. 
%
% Thrsh defines the thresholds, clipping and transparency parameters.
% All these inputs are optional and if undefined a default is chosen:
%
%   Thrsh(1) is always the minimum R^2. If that doesn't exist this is ignored.
%
%   Thrsh(2:3) are the minimum and maximum point of whatever measure is used.
%    For eccentricity, mu, and sigma-like measures any values outside the range are clipped to that level.
%    For other measures, values below Thrsh(2) are removed, while values above Thrsh(3) are clipped to maximum.
%
%   Thrsh(4:5) restrict the eccentricity range irrespective of what kind of map it is.
%    (Obviously, this only works if there is actually a 2D visual map present in the data).
%   
%   Thrhs(6) defines the proporion of the range beyond which the map is scaled to be transparent.
%    If R^2 is present, this proportion is relative to the range between Thrsh(1) and 1.
%    If no R^2 is present, this proportion refers to the range between Thrsh(2) and Thrsh(3).
%   If Thrsh(6) is negative, the same transparency level is used uniformly for the whole map.
%   To turn off transparency, set Thrsh(6) to zero.
%
% Paths is a cell array that defines the filenames of the paths to be displayed. If this is 
%   undefined a dialog box opens allowing you to select the file (close it if none needed).
%  If the last entry in this array is a 1x3 vector, this defines the path colour.
%  If the last entry in this array is NaN, then a default path colour is used.
%  If this is empty, or the last entry is a filename, then the path colour
%   is automatically defined as the opposite polarity of the underlying colour.
%
% If Srf.Data contains more than one subject in the third dimension then
%  another dialog box is opened to select the subject you want to display.
%
% CamView is a two-element vector defining the camera position. Uses the
%  default values defined in SamSrf_defaults. If they aren't included there, 
%  it defaults to pointing to early visual cortex (assuming a FreeSurfer mesh).
%
% MapType is a string containing the name of the data entry to display
%  (e.g. 'Polar', 'Eccentricity', 'R^2', etc.). It can also be a scalar in
%  which case it will use this as the index for Srf.Data and Srf.Values.
%
% PatchHandle is the handle to the patch with the mesh returned by the function.
%  By adding this to the input arguments the figure will simply update the
%  maps instead of opening a whole new figure.
%
% The colour schemes for maps must be defined as strings in SamSrf_defaults.mat.
%
% 10/08/2018 - SamSrf 6 version (DSS)
% 18/06/2019 - Values below Thrsh(2) are now removed unless it's sigma or eccentricity (DSS)
% 22/06/2019 - Changed how colours are assigned (DSS)
%              Added transparency option (DSS)
% 03/07/2019 - Path colour is now by default complementary colour (DSS)
%              Added option to define path colour in Paths (DSS)
% 14/07/2019 - Added option to have uniform transparency level (DSS)
%

%% Create global variables
global Vertices Type Data

%% Load default parameters?
load('SamSrf_defaults.mat');

%% Default thresholds
if nargin < 3
    Thrsh = 0;
end
for i = length(Thrsh)+1:5
    if mod(i,2)
        % Odd thresholds set to Inf
        Thrsh(i) = Inf;
    else
        % Even thresholds set to 0
        Thrsh(i) = 0;
    end
end
% Default transparency
if length(Thrsh) < 6
    Thrsh(6) = 0.3;
end
% If uniform alpha desired
if Thrsh(6) < 0
    % All good data have equal alpha
    UniformAlpha = true;
    Thrsh(6) = 1 + Thrsh(6); % Convert opaqueness to transparency
else
    % By default alphas depend on R^2 or value range
    UniformAlpha = false; 
end

%% Figure handle
fh = gcf;

%% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

%% Is this a group map file?
if size(Srf.Data,3) > 1
    SubjNum = inputdlg(['Which subject? (N=' num2str(size(Srf.Data,3)) ')']);
    SubjNum = cell2mat(SubjNum);
    SubjNum = str2double(SubjNum);
    Srf.Data = Srf.Data(:,:,SubjNum);
end

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
elseif strcmpi(Mesh, 'Fake-Flat')
    Vertices = Srf.Sphere;
    Vertices(:,2) = 0;
else
    error(['Unknown mesh ' Mesh ' specified!']);
end
Vertices(isnan(Srf.Vertices(:,1)),:) = NaN; 
Faces = Srf.Faces;
if length(Curv) == 1
    Curv = repmat(Curv, 1, size(Vertices,1));
end

%% Restrict eccentricities
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
    % Only if both x0 & y0 exist eccentricity can be calculated
    Srf.Data(1,sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2) < Thrsh(4)) = NaN;
    Srf.Data(1,sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2) > Thrsh(5)) = NaN;
end

%% Remove rubbish
if strcmpi(Srf.Values{1}, 'R^2')
    r = Srf.Data(1,:) <= Thrsh(1) | isnan(Srf.Data(1,:));
    Alpha = CalcAlphas(Srf.Data(1,:), [Thrsh(1) Thrsh(1) + (1-Thrsh(1))*Thrsh(6)]); % Transparency based on R^2
    % If uniform alpha desired
    if UniformAlpha
        Alpha(Alpha > 0) = Thrsh(6);
    end
else 
    r = isnan(Srf.Data(1,:));
end

%% Select data type
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
    Values = Srf.Values;
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
else
    Values = Srf.Values;
end
if nargin < 6
    dt = listdlg('ListString', Values, 'SelectionMode', 'single');
    if isempty(dt)
        return
    end
else
    if isscalar(MapType)
        dt = MapType;
    else
        dt = find(strcmpi(Values, MapType));
    end
end
if isempty(dt)
    close(gcf);
    error('Invalid map type specified!');
end
Type = Values{dt};

%% Select paths if desired
if nargin < 4
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

%% Load paths
if ~iscell(Paths)
    Paths = {Paths};
end
PathColour = Inf; % Default path colour is opposite polarity 
if ~isempty(Paths) 
    Vs_paths = [];
    for i = 1:length(Paths)
        % Is a path colour defined?
        if i == length(Paths) && ~ischar(Paths{i})
            if isnan(Paths{i})
                PathColour = NaN;
            else
                PathColour = Paths{i};
            end
        else % No path colour defined
            if strfind(Paths{i}, '.label')
                % Load label
                Label = samsrf_loadlabel(Paths{i}(1:end-6));
                Vs_paths = [Vs_paths; samsrf_borderpath(Srf, Label)];
            else
                % Load paths
                Vs_paths = [Vs_paths; samsrf_loadpaths(Paths{i})];
            end
        end
    end
else
    % No paths needed
    Vs_paths = [];
end

%% Calculate colours for the map 
if strcmpi(Type, 'Polar') 
    % Polar map
    Pha = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    Data = Pha;
    if Srf.Hemisphere(1) == 'l'
        Pha = mod(ceil(Pha + 270), 360) + 1;
    else
        Pha = -Pha;
        Pha = mod(ceil(Pha + 90), 360) + 1;
    end
    Pha(Pha == 0) = 360;
    Pha(r) = 360;
    
    % Colourmap
    cstr = ['colormap(' def_cmap_angle '(360));'];
    Cmap = eval(cstr);        
    
    % Determine colours
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [1 1 1];
    end
    
elseif strcmpi(Type, 'Phase') || strcmpi(Type, 'Phi') 
    % Phase map
    Pha = Srf.Data(dt,:);
    Data = Pha;
    if Srf.Hemisphere(1) == 'r'
        Pha = -Pha;
    end
    Pha = mod(ceil(Pha + 270), 360) + 1;
    Pha(Pha == 0) = 360;
    Pha(r) = 360;
    
    % If no R^2 present
    if ~strcmpi(Srf.Values{1}, 'R^2')
        Alpha = CalcAlphas(Pha, [Thrsh(2) Thrsh(2) + (Thrsh(3)-Thrsh(2))*Thrsh(6)]); % Transparency based on Pha
        % If uniform alpha desired
        if UniformAlpha
            Alpha(Alpha > 0) = Thrsh(6);
        end
    end
    
    % Colourmap
    cstr = ['[colormap(' def_cmap_angle '(360)); CurvGrey];'];
    Cmap = eval(cstr);        
    
    % Determine colours
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    Thrsh = [Thrsh Inf];
    if isnan(PathColour) 
        PathColour = [1 1 1];
    end
    
elseif strcmpi(Type, 'Eccentricity')
    % Eccentricity map
    Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    Data = Rho;
    Rho(r) = 0;
        
    % Set all below minimum to minimum
    Rho(Rho < Thrsh(2)) = Thrsh(2);
    % Adjust minimum
    AdjThr = Thrsh(3) - Thrsh(2);  
    Rho = Rho - Thrsh(2);
    Rho(Rho < 0) = 0;
    % Set all above maximum to maximum
    Rho = Rho / AdjThr;
    Rho(Rho > 1) = 1;   
    % Convert to integers
    Pha = round(Rho * 360);

    % Colourmap
    cstr = ['[colormap(' def_cmap_eccen '(360)); CurvGrey];'];
    Cmap = eval(cstr);        
    
    % Determine colours
    Pha = mod(Pha, 360);
    Pha(Pha==0) = 360;
    Pha(r|isnan(Pha)) = 360;
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [0 0 0];
    end
    
elseif strcmpi(Type, 'Mu') 
    % Mu map
    Mu = Srf.Data(dt,:);
    Data = Mu;
    Mu(r) = 0;
    
    % Set all below minimum to minimum
    Mu(Mu>0 & Mu<+Thrsh(2)) = +Thrsh(2);
    Mu(Mu<0 & Mu>-Thrsh(2)) = -Thrsh(2);
    % Adjust minimum
    AdjThr = Thrsh(3) - Thrsh(2);  
    Mu(Mu>0) = Mu(Mu>0) - Thrsh(2);
    Mu(Mu<0) = Mu(Mu<0) + Thrsh(2);
    % Set all above maximum to maximum
    Mu(Mu>0 & Mu>+Thrsh(3)) = +AdjThr;
    Mu(Mu<0 & Mu<-Thrsh(3)) = -AdjThr;
    % Convert to integers
    Pha = round(Mu / AdjThr * 100) + 100;
    
    % Colourmap
    cstr = ['colormap(' def_cmap_angle '(200));'];
    Cmap = eval(cstr);        
        
    % Determine colours
    Pha(Pha==0) = 1;
    Pha(r|isnan(Pha)|isinf(Pha)) = 200;
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [0 0 0];  
    end
    
elseif strcmpi(Type, 'Sigma') || strcmpi(Type, 'Fwhm') || strcmpi(Type, 'Visual Area') || strcmpi(Type, 'Spread') ...
        || strcmpi(Type, 'Centre') || strcmpi(Type, 'Surround') || strcmpi(Type, 'Sigma1') || strcmpi(Type, 'Sigma2')
    % pRF size map
    Sigma = Srf.Data(dt,:);
    Data = Sigma;
    Sigma(r) = 0;
    
    % Set all below minimum to minimum
    Sigma(Sigma < Thrsh(2)) = Thrsh(2);
    % Adjust minimum
    AdjThr = Thrsh(3) - Thrsh(2);
    Sigma = Sigma - Thrsh(2);
    Sigma(Sigma < 0) = 0;
    % Set all above maximum to maximum
    Sigma = Sigma / AdjThr;
    % Convert to integers
    Pha = round(Sigma * 200);
    
    % Colourmap
    cstr = ['[colormap(' def_cmap_sigma '(200)); CurvGrey];'];
    Cmap = eval(cstr);        
    
    % Determine colours
    Pha(Pha > 200) = 200;
    Pha(Pha <= 0) = 1;
    Pha(r|isnan(Pha)) = 200 + Curv(r|isnan(Pha));
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [1 0 1];
    end
    
else
    % Any generic data
    X = Srf.Data(dt,:);
    % Ensure Field Sign is binary
    if strcmpi(Type, 'Field Sign')
        X = sign(X);
        X(X == 0) = -1;
    end
    % If CMF use logarithm
    if strcmpi(Type, 'Cmf') 
        X = log(X);
        X(X < 0) = 0;
    end
    Data = X;
    X(r) = 0;
    
    % If no R^2 present
    if ~strcmpi(Srf.Values{1}, 'R^2')
        Alpha = CalcAlphas(X, [Thrsh(2) Thrsh(2) + (Thrsh(3)-Thrsh(2))*Thrsh(6)]); % Transparency based on X
        % If uniform alpha desired
        if UniformAlpha
            Alpha(Alpha > 0) = Thrsh(6);
        end
    end
    
    % Set all below minimum to minimum
    X(X>0 & X<+Thrsh(2)) = NaN;
    X(X<0 & X>-Thrsh(2)) = NaN;
    % Adjust minimum
    AdjThr = Thrsh(3) - Thrsh(2);
    X(X>0) = X(X>0) - Thrsh(2);
    X(X<0) = X(X<0) + Thrsh(2);
    % Set all above maximum to maximum
    X(X>0 & X>+AdjThr) = +AdjThr;
    X(X<0 & X<-AdjThr) = -AdjThr;
    % Convert to integers
    Pha = round(X / AdjThr * 100) + 100;
    
    % Colormap
    cstr = ['colormap(' def_cmap_other '(200));'];
    Cmap = eval(cstr);
    
    % Determine colours
    Pha(Pha==0) = 1;
    Pha(r|isnan(Pha)|isinf(Pha)) = 100;
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [0 1 1];
    end
end

%% Draw paths
if isinf(PathColour)
    Colours(Vs_paths,:) = 1 - Colours(Vs_paths,:);
else
    Colours(Vs_paths,:) = repmat(PathColour, length(Vs_paths), 1);
end

%% Display the mesh
set(fh, 'name', [Type ' (' num2str(Thrsh(2)) ' -> ' num2str(Thrsh(3)) ')'], 'color', 'w');
if nargin < 7
    % Draw a new patch
    PatchHandle = patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours, 'FaceColor', 'interp', 'EdgeColor', 'none');
else
    % Simply replace the colours
    set(PatchHandle, 'FaceVertexCData', Colours); 
end
axis off;
ax = gca;
ax.Clipping = 'off';
if nargin < 5 || isempty(CamView)
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
end
set(gca, 'view', CamView(1:2));
zoom(CamView(3));
set(gca, 'projection', 'perspective');
daspect([1 1 1]); % Correct aspect ratio

%% Vertex selection hack
dcm_obj = datacursormode(fh);
set(dcm_obj, 'UpdateFcn', @getvertexfun)

%% Close request function
crfcn = @closereq;
set(fh, 'CloseRequestFcn', crfcn);


% Return transparency levels 
function Alpha = CalcAlphas(VxData, Levels)
Alpha = abs(VxData) - Levels(1); % Set opaqueness below threshold to 0
Alpha = Alpha / (Levels(2)-Levels(1)); % Clip upper level to 1
Alpha(Alpha > 1) = 1;
Alpha(Alpha < 0) = 0;
Alpha(isnan(Alpha)) = 0;
% Replicate into matrix
Alpha = repmat(Alpha',1,3);


function txt = getvertexfun(empt, event_obj)
global Vertices Type Data

pos = get(event_obj, 'Position');
v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1);
txt = {['Vertex: ', num2str(v)];...
       [Type ': ' num2str(Data(v))]};
return
   

function closereq(src, evnt)
% Clear global variables when figure is closed
clear global Vertices Type Data
delete(src);