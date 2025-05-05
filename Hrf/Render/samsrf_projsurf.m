function OutSrf = samsrf_projsurf(Srf, Mesh, Img, Eccen, Thrsh, CamView)
%
% OutSrf = samsrf_surf(Srf, Mesh, Img, Eccen, [Thrsh=0.01, CamView])
%
% Displays a cortical mesh (e.g. 'inflated') overlaid with image Img as it
% would appear on the cortical surface using the pRF parameters in Srf 
% (assumes that Srf came from a typical pRF analysis). Img is an RGB image.
% If Img is empty, it assumes the overlay was saved in Srf.Data(2:4,:).
% Eccen defines the maximal eccentricity of the image i.e. the edges on the
% shorter side of the image are assumed to have eccentricity Eccen.
% The optional Thrsh defines the R^2 threshold of vertices to include.
%
% The function returns the Srf with the colours in Srf.Data(2:4,:).
%
% This projection does not use pRF size but only position. It also simply 
% uses the pixel nearest to the pRF position rather than interpolating. 
%
% 20/09/2024 - Fixed for JSON file for default parameters (DSS)
%

%% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

%% Output structure
OutSrf = Srf;

%% Default threshold?
if nargin < 6
    Thrsh = 0.01;
end

%% Create figure handle
fh = figure;

%% Load curvature
Curv = Srf.Curvature;

%% Load mesh
Mesh = lower(Mesh);
Mesh(1) = upper(Mesh(1));
if isfield(Srf, Mesh)
    Vertices = getfield(Srf, Mesh);
elseif exist([MeshFolder filesep Srf.Hemisphere '.' lower(Mesh)], 'file')
    Vertices = fs_read_surf([MeshFolder filesep Srf.Hemisphere '.' lower(Mesh)]);
else
    Vertices = Srf.Vertices;
end
Faces = Srf.Faces;
if length(Curv) == 1
    Curv = repmat(Curv, 1, size(Vertices,1));
end

%% Remove rubbish
r = Srf.Data(1,:) <= Thrsh | isnan(Srf.Data(1,:));

%% Overlay image 
if isempty(Img)
    % Use saved overlay
    Colours = Srf.Data(2:4,:)';
else
    % Calculate overlay 
    if isa(Img,'uint8') 
        Img = double(Img);
        Img = Img / 255;
    end
    dims = size(Img);
    if dims(3) == 1
        % If greyscale convert to RGB
        repmat(Img, [1 1 3]);
    end
    % Image space grid
    [iy ix] = ndgrid(1:dims(1), 1:dims(2));
    iy = flipud(iy); % Flip Y because of matlab convention
    % Relative to centre in pixels   
    ix = ix - dims(2)/2; 
    iy = iy - dims(1)/2;
    ihw = min(dims(1:2))/2; % Image half width
    % Relative to eccentricity
    ix = ix / ihw * Eccen;
    iy = iy / ihw * Eccen;

    % Calculate overlay
    Colours = NaN(size(Vertices,3),3);
    samsrf_disp('Projecting pixels to cortical surface...');
    samsrf_disp(' Please stand by...');
    parfor v = 1:length(r)
        if ~r(v)
            % Convert image space into visual space
            vx = Srf.Data(2,v); % X-position of pRF
            vy = Srf.Data(3,v); % Y-position of pRF
            ed = sqrt((ix-vx).^2+(iy-vy).^2); % Euclidean distance of pixels from pRF centre
            [ir,ic] = find(ed == min(ed(:))); % Position of pRF in image matrix
            if ir > 0 && ir <= size(Img,1) && ic > 0 && ic <= size(Img,2)
                % Inside the image so gets a colour
                Colours(v,:) = squeeze(Img(ir,ic,:))';
            else
                % Outside the image so set to curvature
                if Curv(v) >= 0
                    Colours(v,:) = [.2 .2 .2];
                else
                    Colours(v,:) = [.8 .8 .8];
                end
            end
        end
    end
end
% Set bad vertices to curvature only
Colours(r&Curv<0,:) = repmat([.8 .8 .8],sum(r&Curv<0),1);
Colours(r&Curv>=0,:) = repmat([.2 .2 .2],sum(r&Curv>=0),1);
% Store in OutSrf
OutSrf.Data = [OutSrf.Data(1,:); Colours'];
if isfield(OutSrf, 'Raw_Data')
    OutSrf = rmfield(OutSrf, 'Raw_Data');
end
OutSrf.Values = {'R^2'; 'Red'; 'Green'; 'Blue'}; 

%% Display the mesh
set(fh, 'name', 'Cortical image', 'color', 'w');
patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), 'FaceVertexCData', Colours, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis off;
ax = gca;
ax.Clipping = 'off';
if nargin < 7
    SamSrfDefs = LoadSamSrfDefaults;
    if ~isfield(SamSrfDefs, 'def_views')
        samsrf_disp('WARNING: def_views not defined in SamSrf_defaults.json');
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
            CamView = SamSrfDefs.def_views(:,1)';
        elseif Srf.Hemisphere(1) == 'r'
            % Right hemisphere
            CamView = SamSrfDefs.def_views(:,2)';
        else
            % Both hemispheres
            if size(SamSrfDefs.def_views,2) > 2
                CamView = SamSrfDefs.def_views(:,3)';
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
samsrf_lighting('on');

