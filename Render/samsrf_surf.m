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
%   Thrsh(6) defines the proportion of the range beyond which the map is scaled to be transparent.
%    If R^2 is present, this proportion is relative to the range between Thrsh(1) and 1.
%    If no R^2 is present, this proportion refers to the range between Thrsh(2) and Thrsh(3).
%   If Thrsh(6) is negative, the same transparency level is used uniformly for the whole map.
%   To turn off transparency, set Thrsh(6) to zero.
%
% Paths is a cell array that defines the filenames of the paths to be displayed. If this is 
%   undefined a dialog box opens allowing you to select the file (close it if none needed).
%  If the last entry in this array is a 1x3 vector, this defines the path colour.
%  If the last entry in this array is NaN, then a default path colour is used.
%  If Paths is empty, or the last entry is a filename, then the path colour
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
% 17/07/2020 - SamSrf 7 version (DSS)
% 07/08/2020 - Changed default camera view for right hemisphere (DSS)
%              Added support for displaying ROI numbers (DSS)
% 16/10/2020 - Fixed bug with determining transparency when no goodness-of-fit exists (DSS)
% 29/10/2020 - Polar/phase colour schemes now account for bilateral data files (DSS)
% 26/11/2020 - Added default camera angle for bilateral data files (DSS)  
%

%% Create global variables
global Vertices Type Data

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
    Thrsh(6) = 0.1;
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
        disp('You must pick a subject! Selecting #1...');
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
% Does file contain connective field profiles?
if isfield(Srf, 'ConFlds')
    Values{end+1} = 'Connective Field';
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
    CfVx = 0;
else
    if isscalar(MapType)
        % Map type given by number
        dt = MapType;
        % CF profiles present?
        if isfield(Srf, 'ConFlds')
            CfVx = 0;
        end
    else
        % Map type given by name
        if contains(MapType, 'Connective Field')
            hash = strfind(MapType,'#');
            if ~isempty(hash)
                CfVx = str2double(MapType(hash+1:end));
                if isnan(CfVx)
                    error('Invalid vertex number provided!');
                end
            else
                CfVx = 0;
            end
            MapType = MapType(1:16);
        end
        dt = find(strcmpi(Values, MapType));
    end
end
if isempty(dt)
    close(gcf);
    error('Invalid map type specified!');
end
Type = Values{dt};

%% If CF profile selected
if strcmpi(Type, 'Connective Field')
    % Select vertex?
    if CfVx == 0
        CfVx = inputdlg('Vertex #', ['Which vertex? (1-' num2str(size(Srf.Data,2)) ')']);
        CfVx = str2double(cell2mat(CfVx));
        if isnan(CfVx)
            error('Invalid vertex number provided!');
        end
    end
    % Overlay connective field profile onto cortex
    Srf.Data(end+1,:) = zeros(1,size(Srf.Data,2)); % Empty surface
    Srf.Data(end,Srf.SeedVx) = Srf.ConFlds(:,CfVx)'; % Fill in seed ROI with CF profile
    Srf.Values{1} = ''; % Remove R^2 value label
    Srf.Values{end+1} = 'Connective Field'; % Remove values field
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
                if isnan(Paths{i})
                    PathColour = NaN;
                else
                    PathColour = Paths{i};
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
    cstr = ['[colormap(' def_cmap_eccen(2:end) '(360)); CurvGrey];'];
    Cmap = eval(cstr);        
    if def_cmap_eccen(1) == '-'
        Cmap = flipud(Cmap);
    end
    
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
    
elseif strcmpi(Type, 'Sigma') || strcmpi(Type, 'Fwhm') || strcmpi(Type, 'Visual Area') || strcmpi(Type, 'Spread') || strcmpi(Type, 'ROI') ...
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
    cstr = ['colormap(' def_cmap_other(2:end) '(200));'];
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
global Vertices Type Data Srf

pos = get(event_obj, 'Position');
v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1);
txt = {['Vertex: ', num2str(v)];...
       [Type ': ' num2str(Data(v))]};
return
   

function closereq(src, evnt)
% Clear global variables when figure is closed
clear global Vertices Type Data
delete(src);