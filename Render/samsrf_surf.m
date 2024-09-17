function PatchHandle = samsrf_surf(Srf, Mesh, Thrsh, Paths, CamView, MapType, PatchHandle)
%
% PatchHandle = samsrf_surf(Srf, Mesh, [Thrsh, Paths, CamView, MapType, PatchHandle])
%
% Displays a cortical surface Mesh (e.g. 'Inflated') overlaid with the data in Srf. 
%   If the surface mesh isn't found in the Srf, the WM-GM surface is displayed.
%
% If the Srf contains a reverse correlation or connective field data, a second figure 
%   is opened to display the relevant data for the currently selected vertex.
%   However, this functionality only works when using the DisplayMaps tool.
%
% Note: This function is primarily designed to be called internally by the DisplayMaps tool!
%       While it can be used directly, it is overly complex. You may want to create 
%       a streamlined, simplified version for displaying maps directly.
%
% Optional input arguments:
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
%   Thrsh(6) defines the proportion of the range beyond which the map is scaled to be transparent.
%    If R^2 is present, this proportion is relative to the range between Thrsh(1) and 1.
%    If no R^2 is present, this proportion refers to the range between Thrsh(2) and Thrsh(3).
%   If Thrsh(6) is negative, the same transparency level is used uniformly for the whole map.
%   To turn off transparency, set Thrsh(6) to zero. This is the default.
%
% Paths is a cell array that defines the filenames of the paths to be displayed. If this is 
%   undefined a dialog box opens allowing you to select the file (close it if none needed).
%  If the last entry is NaN followed by a 1x3 vector, this defines the path colour.
%  If the last entry is simply NaN, then a default path colour is used.
%  If Paths is empty, or the last entry is a filename, then the path colour
%   is automatically defined as the opposite polarity of the underlying colour.
%
% If Srf.Data contains more than one subject in the third dimension then
%  another dialog box is opened to select the subject you want to display.
%
% CamView is a three-element vector. The first two values define the camera position. 
%  The third value is the zoom factor. If undefined or empty, it uses the default values 
%  defined in SamSrf_defaults. If they aren't included there, it defaults to pointing to 
%  early visual cortex (assuming a FreeSurfer mesh).
%
% MapType determines which map is being displayed. This can be a...
%   - String containing the name of the data entry in Srf.Values to display
%      (e.g. 'R^2', 'Sigma', etc.). If this is a normal pRF map with x0 and y0 
%      coordinates, you can also display 'Polar' and 'Eccentricity' maps.
%   - Scalar (e.g. 5) defining the index for Srf.Data and Srf.Values.
%   - Vector containing the actual map data to display. 
%      This must have the same number of columns as Srf.Data. 
%   - If undefined, the function opens a selection dialog to select from Srf.Values. 
%
% PatchHandle is the handle to the patch with the mesh returned by the function.
%  By adding this to the input arguments the figure will simply update the
%  maps instead of opening a whole new figure. This is mostly for internal
%  use by the DisplayMaps tool.
%
% The colour schemes for maps must be defined as strings in SamSrf_defaults.mat.
%
% 13/03/2022 - Now reports which default parameter file it's loading (DSS)
%              Fixed bug when trying to display white-matter surface (DSS)
% 14/03/2022 - Vertex inspector no longer scales amplitude for reverse correlation profiles (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 29/07/2022 - Turned off transparency by default (DSS)
%              Fixed bug with eccentricity bound clipping (DSS)
% 10/08/2022 - Now ensures no NaNs in displayed maps to prevent crash (DSS)
% 12/08/2022 - Vertex inspector can display visual CF profiles now (DSS)
% 27/09/2022 - Plots of CFs in visual space are now zoomed in (DSS)
% 19/06/2023 - Added option for transparency in surface mesh (DSS)
% 03/10/2023 - Removed overly verbose defaults message (DSS)
%              R^2 maps now use eccentricity colour scheme (DSS)
% 12/02/2024 - Fixed bug when using 32bit data with CF maps (DSS)  
%

%% Create global variables
global Vertices Type Data CurvGrey fh fv pv

%% Load default parameters?
load('SamSrf_defaults.mat');
% Ensure colour maps have sign
if def_cmap_angle(1) ~= '-' && def_cmap_angle(1) ~= '+'
    def_cmap_angle = ['+' def_cmap_angle];
end
if def_cmap_eccen(1) ~= '-' && def_cmap_eccen(1) ~= '+'
    def_cmap_eccen = ['+' def_cmap_eccen];
end
if def_cmap_other(1) ~= '-' && def_cmap_other(1) ~= '+'
    def_cmap_other = ['+' def_cmap_other];
end
if def_cmap_sigma(1) ~= '-' && def_cmap_sigma(1) ~= '+'
    def_cmap_sigma = ['+' def_cmap_sigma];
end
if ~exist('def_logmaps', 'var')
    def_logmaps = false;
end  

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
    Thrsh(6) = 0;
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

%% Expand (if necessary) & denoise Srf 
Srf = samsrf_expand_srf(Srf);

%% Multi-subject Srf?
if size(Srf.Data,3) > 1
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
end

%% Load curvature
Curv = Srf.Curvature;
if exist('def_curv', 'var')
    if strcmpi(def_curv, 'Greyscale')
        % Greyscale for curvature
        CurvGrey = gray(11); % Grey scale colour map
        CurvGrey = CurvGrey(1:10,:); % Remove white to avoid transparency issues 
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
    CurvGrey = CurvGrey(1:10,:); % Remove white to avoid transparency issues
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
elseif strcmpi(Mesh, 'Inflated-50%')
    Vertices = Srf.Vertices + (Srf.Inflated-Srf.Vertices)*0.5;
elseif strcmpi(Mesh, 'Pial-50%')
    Vertices = Srf.Pial + (Srf.Inflated-Srf.Pial)*0.5;
else
    if strcmpi(Mesh, 'White')
        Vertices = Srf.Vertices;
    else
        samsrf_error(['Unknown mesh ' Mesh ' specified!']);
    end
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
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
    r = Srf.Data(1,:) <= Thrsh(1) | isnan(Srf.Data(1,:));
    Srf.Data(1,r) = NaN; % Set rubbish to NaN
    Alpha = CalcAlphas(Srf.Data(1,:), [Thrsh(1) Thrsh(1) + (1-Thrsh(1))*Thrsh(6)]); % Transparency based on R^2
    % If uniform alpha desired
    if UniformAlpha
        Alpha(Alpha > 0) = Thrsh(6);
    end
else 
    r = isnan(Srf.Data(1,:));
end

%% Select data type
Values = Srf.Values;
% If desired, open second figure for data inspector
if (isfield(Srf, 'ConFlds') || isfield(Srf, 'Rmaps') || isfield(Srf, 'Model')) && nargin < 7
    CallingFuncs = dbstack;
    % Only if using DisplayMaps tool
    if ~strcmpi(CallingFuncs(end).file, 'samsrf_surf.m')
        fv = figure; % Visualization figure handle
        pv = []; % Visualisation patch handle
        figure(fh);
    end
end
% Does file contain retinotopic maps?
if length(Srf.Values) > 1 && strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
    Values(end+1:end+2) = {'Polar'; 'Eccentricity'};
end
if nargin < 6
    dt = listdlg('ListString', Values, 'SelectionMode', 'single');
    if isempty(dt)
        return
    end
else
    if isvector(MapType) && length(MapType) == size(Srf.Data,2)
        % Map provided as data vector
        Type = 'Custom';
        dt = NaN;
    elseif isscalar(MapType)
        % Map type given by number
        dt = MapType;
    elseif ischar(MapType)
        % Map type given by name
        dt = find(strcmpi(Values, MapType));
    else
        samsrf_error('Invalid map selected!');
    end
end
if isempty(dt)
    close(gcf);
    samsrf_error('Invalid map type specified!');
end
if ~isnan(dt)
    Type = Values{dt};
end

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
PathColour = Inf; % Default path colour is opposite polarity 
% If paths defined
if ~isempty(Paths) 
    % Ensure cell array
    if ~iscell(Paths)
        Paths = {Paths};
    end
    % Vector of path vertices
    Vs_paths = [];
    % Loop thru paths
    for i = 1:length(Paths)
        if ~ischar(Paths{i}) && ~isnan(Paths{i}(1))            
            % If paths contain vertices
            Vs_paths = [Vs_paths; Paths{i}];
        else
            % If paths contain strings
            if i == length(Paths) && ~ischar(Paths{i})
                % Is a path colour defined?
                if isnan(Paths{i}(1)) && length(Paths{i}) == 1
                    PathColour = NaN;
                elseif isnan(Paths{i}(1))
                    PathColour = Paths{i}(2:4);
                end
            else
                % No path colour defined
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
    end
    % Default path thickness defined?
    if exist('def_pathwidth', 'var')
        for i = 1:def_pathwidth-1
            Vs_paths = [Vs_paths; samsrf_neighbours(Vs_paths, Srf.Faces)];
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
    Pha(isnan(Pha)) = 0;
    Data = Pha;
    if Srf.Hemisphere(1) == 'l'
        % Left hemisphere
        Pha = mod(ceil(Pha + 270), 360) + 1;
    elseif Srf.Hemisphere(1) == 'b'
        % Bilateral maps
        Pha(1:Srf.Nvert_Lhem) = mod(ceil(Pha(1:Srf.Nvert_Lhem) + 270), 360) + 1;
        Pha(Srf.Nvert_Lhem+1:end) = -Pha(Srf.Nvert_Lhem+1:end);
        Pha(Srf.Nvert_Lhem+1:end) = mod(ceil(Pha(Srf.Nvert_Lhem+1:end) + 90), 360) + 1;
    else
        % Right hemisphere
        Pha = -Pha;
        Pha = mod(ceil(Pha + 90), 360) + 1;
    end
    Pha(Pha == 0) = 360;
    Pha(r) = 360;
    
    % Colourmap
    cstr = ['colormap(' def_cmap_angle(2:end) '(360));'];
    Cmap = eval(cstr);        
    if def_cmap_angle(1) == '-'
        Cmap = flipud(Cmap);
    end
    
    % Determine colours
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [1 1 1];
    end
    
elseif strcmpi(Type, 'Phase') || strcmpi(Type, 'Phi') 
    % Phase map
    Pha = Srf.Data(dt,:);
    Pha(isnan(Pha)) = 0;
    Data = Pha;
    if Srf.Hemisphere(1) == 'r'
        % Right hemisphere
        Pha = -Pha;
    elseif Srf.Hemisphere(1) == 'b'
        % Bilateral data
        Pha(Srf.Nvert_Lhem+1:end) = -Pha(Srf.Nvert_Lhem+1:end);
    end
    Pha = mod(ceil(Pha + 270), 360) + 1;
    Pha(Pha == 0) = 360;
    Pha(r) = 360;
    
    % If no R^2 present
    if ~strcmpi(Srf.Values{1}, 'R^2') && ~strcmpi(Srf.Values{1}, 'nR^2')
        Alpha = CalcAlphas(Pha, [Thrsh(2) Thrsh(2) + (Thrsh(3)-Thrsh(2))*Thrsh(6)]); % Transparency based on Pha
        % If uniform alpha desired
        if UniformAlpha
            Alpha(Alpha > 0) = Thrsh(6);
        end
    end
    
    % Colourmap
    cstr = ['[colormap(' def_cmap_angle(2:end) '(360)); CurvGrey];'];
    Cmap = eval(cstr);        
    if def_cmap_angle(1) == '-'
        Cmap = flipud(Cmap);
    end
    
    % Determine colours
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    Thrsh = [Thrsh Inf];
    if isnan(PathColour) 
        PathColour = [1 1 1];
    end
    
elseif strcmpi(Type, 'Eccentricity')
    % Eccentricity map
    Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    Rho(isnan(Rho)) = 0;
    Data = Rho;
    Rho(r) = 0;
 
    % If logarithmic maps
    if def_logmaps
    	Rho = log2(Rho);
    end	
        
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
    cstr = ['[colormap(' def_cmap_eccen(2:end) '(360)); CurvGrey];'];
    Cmap = eval(cstr);        
    if def_cmap_eccen(1) == '-'
        Cmap = flipud(Cmap);
    end
    
    % Determine colours
    Pha(Pha==0) = 1;%360;
    Pha(r|isnan(Pha)) = 1;
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [0 0 0];
    end
    
elseif strcmpi(Type, 'Mu') 
    % Mu map
    Mu = Srf.Data(dt,:);
    Mu(isnan(Mu)) = 0;
    Data = Mu;
    Mu(r) = 0;

    % If logarithmic maps
    if def_logmaps
    	Mu = log2(Mu);
    end	
    
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
    cstr = ['colormap(' def_cmap_angle(2:end) '(200));'];
    Cmap = eval(cstr);        
    if def_cmap_angle(1) == '-'
        Cmap = flipud(Cmap);
    end
        
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
    Sigma(isnan(Sigma)) = 0;
    Data = Sigma;
    Sigma(r) = 0;

    % If logarithmic maps
    if def_logmaps
    	Sigma = log2(Sigma);
    end	
    
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
    cstr = ['[colormap(' def_cmap_sigma(2:end) '(200)); CurvGrey];'];
    Cmap = eval(cstr);        
    if def_cmap_sigma(1) == '-'
        Cmap = flipud(Cmap);
    end
    
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
    if strcmpi(Type, 'Custom')
        % Using a custom data vector as map
        X = MapType;
    else
        % Any named or index field in Srf.Data
        X = Srf.Data(dt,:);
    end
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
    X(isnan(X)) = 0;
    Data = X;
    X(r) = 0;
    
    % If no R^2 present
    if ~strcmpi(Srf.Values{1}, 'R^2') && ~strcmpi(Srf.Values{1}, 'nR^2')
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
    if strcmpi(Type, 'R^2') || strcmpi(Type, 'nR^2') 
        cstr = ['colormap(' def_cmap_eccen(2:end) '(200));'];
    else
        cstr = ['colormap(' def_cmap_other(2:end) '(200));'];
    end
    Cmap = eval(cstr);
    if def_cmap_other(1) == '-'
        Cmap = flipud(Cmap);
    end
    
    % Determine colours
    Pha(Pha==0) = 1;
    Pha(r|isnan(Pha)|isinf(Pha)) = 100;
    Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
    if isnan(PathColour) 
        PathColour = [0 1 1];
    end
end

%% Draw paths
if ~isnan(Vs_paths) % Only if not a NaN
    if isinf(PathColour)
        Colours(Vs_paths,:) = 1 - Colours(Vs_paths,:);
    else
        Colours(Vs_paths,:) = repmat(PathColour, length(Vs_paths), 1);
    end
end

%% Display the mesh
if isfield(Srf, 'VxAlpha')
    % Surface has transparency
    SurfAlpha = 'interp';
    VxAlpha = Srf.VxAlpha;
else
    % Surface has no transparency
    SurfAlpha = 1;
    VxAlpha = ones(size(Vertices,1),1);
end
set(fh, 'name', [Type ' (' num2str(Thrsh(2)) ' -> ' num2str(Thrsh(3)) ')'], 'color', 'w');
if nargin < 7
    % Draw a new patch
    PatchHandle = patch('Vertices', Vertices, 'Faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceVertexAlphaData', VxAlpha, 'FaceAlpha', SurfAlpha);
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
            % Left hemisphere
            CamView = [15 -30 1.8];
        elseif Srf.Hemisphere(1) == 'r'
            % Right hemisphere
            CamView = [-13 -38 1.8];
        else
            % Both hemispheres
            CamView = [4 -30 2.2];
        end
    else
        % Use default camera angle
        if Srf.Hemisphere(1) == 'l'
            % Left hemisphere
            CamView = def_views(:,1)';
        elseif Srf.Hemisphere(1) == 'r'
            % Right hemisphere
            CamView = def_views(:,2)';
        else
            % Both hemispheres
            if size(def_views,2) > 2
                CamView = def_views(:,3)';
            else
                CamView = [4 -30 2.2];
            end
        end
    end
end
% Don't update view if simply redrawing 
if nargin < 7
    set(gca, 'view', CamView(1:2));
    zoom(CamView(3));
    set(gca, 'projection', 'perspective');
    daspect([1 1 1]); % Correct aspect ratio
end

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
global Vertices Type Data CurvGrey fh fv pv Srf 

% Datatip update
pos = get(event_obj, 'Position'); % Selected position
v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1); % Selected vertex
txt = {['Vertex: ', num2str(v)];...
       [Type ': ' num2str(Data(v))]}; % Datatip message

% Update visualisation figure?   
if ~isempty(Srf)
    % Which data to display?
    if isfield(Srf, 'Rmaps')
        % Reverse correlation profile
        figure(fv);
        if ~isnan(Srf.Rmaps)
            % pRF profiles stored in Srf
            if var(Srf.Rmaps(:,v)) == 0
                Srf.Rmaps(1,v)=.001;
            end
        end
        samsrf_showprf(Srf, v);
        legend off
        title(['Vertex: ', num2str(v)]);
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [.6 .1 .3 .3]);
        set(fv, 'name', 'Reverse correlation profile');
        figure(fh);
        
    elseif isfield(Srf, 'Model')
        % pRF profile based on fit parameters
        if isfield(Srf.Model, 'Prf_Function')
            figure(fv);
            samsrf_showprf(Srf, v, Srf.Model);
            legend off
            caxis([-1 1]);
            title(['Vertex: ', num2str(v)]);
            set(gcf, 'Units', 'normalized');
            set(gcf, 'Position', [.6 .1 .3 .3]);
            set(fv, 'name', 'pRF model fit');
            figure(fh);
        end
        
    elseif isfield(Srf, 'Y')
        % Predicted & observed time series
        if isfield(Srf, 'Model_')
            figure(fv);
            samsrf_fitvsobs(Srf, Srf.Model_, v);
            set(fv, 'name', 'Observed vs predicted time series');
            figure(fh);
        end
        
    elseif ~isfield(Srf, 'Y') && ~isfield(Srf, 'X') && ~isfield(Srf, 'Y_')
        CurrData = Srf.Data(:,v);
        ts = ['Vertex: ', num2str(v)];
        if strcmpi(Srf.Values{1}, 'Noise Ceiling')
            % If file includes noise ceiling in row one extract that
            ts = [ts ', Noise ceiling: ' num2str(CurrData(1))];
            CurrData = CurrData(2:end);
        end
        figure(fv);
        hold off
        if contains(Srf.Functional, 'GLM contrasts')
            % Plot GLM contrasts
            plot(Srf.Data(:,v), 'o-k', 'linewidth', 2);
            xlim([.5 size(Srf.Data,1)+.5]);
            ylim([nanmin(Srf.Data(:)) nanmax(Srf.Data(:))]);
            line(xlim, [0 0], 'color', [1 1 1]/2, 'linewidth', 2);
            set(gca, 'fontsize', 12, 'xtick', 1:size(Srf.Data,1), 'xticklabel', Srf.Values);
            title(ts);
            xlabel('Contrast');
            ylabel('Differential response');
            set(gcf, 'Units', 'normalized');
            set(gcf, 'Position', [.1 .1 .8 .4]);
            set(fv, 'name', 'GLM contrasts');
            
        else
            % Plotting observed time series
            if isfield(Srf, 'Raw_Data')
                plot(Srf.Raw_Data(:,v), 'color', [.6 .6 1], 'linewidth', 2);
                hold on
            end
            plot(CurrData, 'color', [0 0 1], 'linewidth', 2);
            hold on
            xlim([1 length(CurrData)]);
            line(xlim, [0 0], 'color', [1 1 1]/2, 'linewidth', 2);
            if isfield(Srf, 'Raw_Data')
                legend({'Smooth' 'Raw'});
            end
            title(ts);
            set(gcf, 'Units', 'normalized');
            set(gcf, 'Position', [.1 .1 .8 .4]);
            set(gca, 'fontsize', 12);
            xlabel('Volumes (#)');
            ylabel('Response (z)');
            set(fv, 'name', 'Observed time series');
        end
        grid on
        figure(fh);
        
    elseif isfield(Srf, 'TempMap')
        % Connective field visual profile
        figure(fv); % Change to this figure
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [.5 .3 .3 .3]);
        hold off
        pv = []; 
        % Seed ROI correlations
        if isnan(Srf.ConFlds_)
            % Recompute CF correlation profile
            R = corr(Srf.Y_(:,v), Srf.Y_(:,Srf.SeedVx));
        else
            % CF correlation profile saved in file
            R = Srf.ConFlds_(:,v);
        end
        [txy,u] = unique(Srf.TempMap', 'rows'); % Unique coordinates in seed map 
        txy = double(txy); % Ensure double format
        R = R(u); % Remove duplicate correlations
        R(isnan(R)) = 0; % Remove NaNs
        mR = max(R); % Maximum after redundant vertices removed
        % Plot correlation profile
        polarpatch(txy(:,1), txy(:,2), R); % Plot profile
        hold on
        scatter(Srf.Data(2,v), Srf.Data(3,v), 'go'); % Plot centroid
        axis([-1 1 -1 1] * (sqrt(sum(Srf.Data(2:3,v).^2)) + Srf.Data(4,v))*4); % Restrict axes
        xlabel('Horizontal position (deg)');
        ylabel('Vertical position (deg)');
        colorbar
        colormap berlin
        caxis([-1 1] * sqrt(Srf.Data(1,v)));
        cblabel('Correlation');
        title(['Selected vertex: #' num2str(v)]);
        % Plot CF outline
        warning off
        % Region growing for positive subfield
        tri = delaunay(txy); % Delaunay triangles
        pos = samsrf_clusterroi(find(R==mR,1), R, mR/2, tri); % Points in positive subfield
        % Quantify positive subfield
        ps = polyshape(txy(pos,:)); % Points in positive subfield
        ps = ps.convhull; % Convex hull of positive subfield
        % Plot outline
        for i = 1:size(ps.Vertices,1)
            if i == size(ps.Vertices,1)
                line(ps.Vertices([i 1],1), ps.Vertices([i 1],2), 'color', 'g');
            else
                line(ps.Vertices(i:i+1,1), ps.Vertices(i:i+1,2), 'color', 'g');
            end
        end
        % Quantify negative subfield
        nR = min(R); % Minimum correlation
        neg = R < nR/2; % Points in negative subfield
        ps = polyshape(txy(neg,:)); % Points in negative subfield
        ps = ps.convhull; % Convex hull of negative subfield
        % Plot outline
        for i = 1:size(ps.Vertices,1)
            if i == size(ps.Vertices,1)
                line(ps.Vertices([i 1],1), ps.Vertices([i 1],2), 'color', 'g', 'linestyle', '--');
            else
                line(ps.Vertices(i:i+1,1), ps.Vertices(i:i+1,2), 'color', 'g', 'linestyle', '--');
            end
        end    
        % Back to main figure
        hold off
        figure(fh);
        
    elseif isfield(Srf, 'ConFlds')
        % Connective field profile
        curcam = get(gca, 'view');
        figure(fv);
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [.5 .3 .3 .3]);
        zoom(gca, 'reset');
        % Curvature
        Curv = Srf.Curvature;
        Curv(isnan(Curv)) = 0;
        Curv = -Curv + 0.5;
        Curv(Curv <= 0) = 0.000001;
        Curv(Curv > 1) = 1;
        Curv = ceil(Curv * size(CurvGrey,1));
        % Seed ROI correlations
        if isnan(Srf.ConFlds)
            % Recompute CF correlation profile
            Cf = corr(Srf.Y_(:,v), Srf.Y_(:,Srf.SeedVx));
        else
            % CF correlation profile saved in file
            Cf = Srf.ConFlds(:,v);
        end
        CfThrsh = [.01 max(Cf(:))];
        X = NaN(1,size(Vertices,1));
        X(1,Srf.SeedVx) = Cf;
        % Set all below minimum to minimum
        X(X>0 & X<+CfThrsh(1)) = NaN;
        X(X<0 & X>-CfThrsh(1)) = NaN;
        % Adjust minimum
        AdjThr = CfThrsh(2) - CfThrsh(1);
        X(X>0) = X(X>0) - CfThrsh(1);
        X(X<0) = X(X<0) + CfThrsh(1);
        % Set all above maximum to maximum
        X(X>0 & X>+AdjThr) = +AdjThr;
        X(X<0 & X<-AdjThr) = -AdjThr;
        % Convert to integers
        Pha = round(X / AdjThr * 100) + 100;    
        Pha(Pha==0) = 1;
        Pha(isnan(Pha)|isinf(Pha)) = 100;  
        % Transparency
        Alpha = zeros(size(Vertices,1),3);
        Alpha(Srf.SeedVx,:) = repmat(X(Srf.SeedVx)'/AdjThr * 0.5 + 0.5, [1 3]);
        Alpha(isnan(Alpha)) = 0.5; % To remove white artifacts
        % Colour map
        Cmap = hotcold(200);
        Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
        % Draw surface
        set(gca, 'view', curcam);
        if isempty(pv) || ~isvalid(pv)
            % New plot
            cla(gca);
            pv = patch('vertices', Vertices, 'faces', Srf.Faces(:,[1 3 2]), 'FaceVertexCData', Colours, 'FaceColor', 'interp', 'EdgeColor', 'none');
            colorbar off
            zoom(3);
            axis off;
            ax = gca;
            ax.Clipping = 'off';
            set(gca, 'projection', 'perspective');
            daspect([1 1 1]); % Correct aspect ratio    
        else
            % Simply redraw
            set(pv, 'FaceVertexCData', Colours); 
        end
        title(['Vertex: ', num2str(v)]);
        set(fv, 'name', 'Connective field profile');
        figure(fh);
    end
end
return
   

function closereq(src, evnt)
% Clear global variables when figure is closed
global fv 
% Just in case it's already closed
try
    close(fv);
catch
    true; 
end
clear global Vertices Type Data CurvGrey fh fv pv
delete(src);
