function samsrf_benson2srf(gii, surfdir, anatpath)
% 
% samsrf_benson2srf(gii, surfdir, [anatpath])
%
% Converts a Benson-style map prediction GII file into a SamSrf surface file.
%   (Need to select the 'all' file containing polar, eccentricity, and ROIs)
% Also saves the ROI labels for V1-V3 in a subfolder called ROIs_Benson.
%
% This is using the original Benson templates that only contained V1-V3.
% A future version could adjust this to use a newer template with a larger number 
% of ROIs as well as predictions for pRF size.
%
%   gii:        Name of GII file (without extension)
%   surfdir:    Folder containing the subject's surface data
%   anatpath:   The folder where anatomical surfaces to stored
%               (See also other surface projection functions)
%               Defaults to '' so anatomy is not split off!
%
% Note that this function requires gifti functionality from SPM12.
%
% 18/09/2024 - Now expects maps in GII instead of MGH format (DSS)
%

if nargin < 3
    anatpath = '';
end

%% Determine file parts 
[pn, gii] = fileparts(gii);  % Split file name into components
sl = strfind(gii, '_');
hemis = gii(1:sl-1);
gii = gii(4:end);
if isempty(pn)
    pn = '.';
end

%% Data labels
valstrs = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'ROI'};

%% Load surface vertices
[V0 F] = fs_read_surf([surfdir filesep hemis '.white']); % Grey-white surface
P = fs_read_surf([surfdir filesep hemis '.pial']); % Pial surface
I = fs_read_surf([surfdir filesep hemis '.inflated']); % Inflated surface
S = fs_read_surf([surfdir filesep hemis '.sphere']); % Spherical surface
C = fs_read_curv([surfdir filesep hemis '.curv']); % Cortical curvature 
A = fs_read_curv([surfdir filesep hemis '.area']); % Cortical surface area
T = fs_read_curv([surfdir filesep hemis '.thickness']); % Cortical thickness
N = P - V0; % Cortical vectors for each vertex 

%% Use surf folder name as Structural field
strimg = surfdir;
samsrf_disp(['Saving Benson maps for ' strimg]);

%% Load GII data
D = gifti([pn filesep hemis '_' gii '.gii']);
D = D.cdata;
D(:,1) = -D(:,1) + 90;
% Convert to Cartesian coordinates
x0 = D(:,2) .* cos(D(:,1)/180*pi);
y0 = D(:,2) .* sin(D(:,1)/180*pi);
rois = D(:,3);
% If right hemisphere flip sign
if strcmp(hemis, 'rh')
    x0 = -x0;
end

%% Create surface structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Structural = strimg;
Srf.Functional = gii;
Srf.Hemisphere = hemis;
Srf.Cortex_Steps = NaN;
Srf.Vertices = V0;
Srf.Pial = P;
Srf.Inflated = I;
Srf.Sphere = S;
Srf.Faces = F;
Srf.Normals = N;
Srf.Curvature = C';
Srf.Area = A';
Srf.Thickness = T';
% Add map data
Srf.Data = NaN(5, size(V0,1));
Srf.Data(1,:) = zeros(1,size(V0,1));
Srf.Data(1,rois~=0) = ones(1, sum(rois~=0)); % Label predicted maps as R^2 = 1
Srf.Data(2,:) = x0';
Srf.Data(3,:) = y0';
Srf.Data(5,:) = rois';
Srf.Rule = 'X';
Srf.Values = valstrs;

%% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);

%% Save surface data
save([hemis '_' gii '.mat'], 'Srf', '-v7.3');
samsrf_disp(['Saved ' hemis '_' gii '.mat']);
samsrf_anatomy_srf([hemis '_' gii], anatpath);
samsrf_newline;
% Save ROI labels
mkdir('ROIs_Benson');
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V1'], 4, find(abs(rois) == 1));
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V2'], 4, find(abs(rois) == 2));
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V3'], 4, find(abs(rois) == 3));
