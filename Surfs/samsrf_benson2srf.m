function samsrf_benson2srf(mghimg, surfdir)
% 
% samsrf_benson2srf(mghimg, surfdir)
%
% Converts a Benson-style map prediction MGH file into a SamSrf surface file.
%   (Need to select the 'all' file containing polar, eccentricity, and ROIs)
% Also saves the ROI labels for V1-V3 in a subfolder called ROIs_Benson.
%
%   mghimg:     Name of MGH file (without extension)
%   surfdir:    Folder containing the surface data
%
% Note that this function requires fs_load_mgh.m, a function that was 
% modified from code on the FreeSurfer website.
%
% The function automatically finds the NII file for the T1 in <SubjectID>/mri
% and adds this to Srf.Structural.
%
% 07/08/2018 - SamSrf 6 version (DSS)
%

%% Determine file parts 
[pn, hemis, mghimg] = fileparts(mghimg);  % Split file name into components
mghimg = mghimg(2:end);  % Remove dot from beginning

%% Data labels
valstrs = {'R^2'; 'x0'; 'y0'; 'ROIs'};

%% Load surface vertices
[V0 F] = fs_read_surf([surfdir filesep hemis '.white']); % Grey-white surface
C = fs_read_curv([surfdir filesep hemis '.curv']); % Curvature 
P = fs_read_surf([surfdir filesep hemis '.pial']); % Pial surface
I = fs_read_surf([surfdir filesep hemis '.inflated']); % Inflated surface
S = fs_read_surf([surfdir filesep hemis '.sphere']); % Spherical surface
C = fs_read_curv([surfdir filesep hemis '.curv']); % Cortical curvature 
A = fs_read_curv([surfdir filesep hemis '.area']); % Cortical surface area
T = fs_read_curv([surfdir filesep hemis '.thickness']); % Cortical thickness
N = P - V0; % Cortical vectors for each vertex 

%% Find T1 used for reconstruction
strimg = dir([surfdir filesep '..' filesep 'mri' filesep '*.nii']);
if isempty(strimg)
    strimg = dir([surfdir filesep '..' filesep 'mri' filesep 'orig' filesep '*.nii']);
end
strimg = strimg(1).name;
if isempty(strimg)
    error('Couldn''t find any NII file for the T1!');
end
disp(['Saving Benson maps for ' strimg]);

%% Load MGH data
D = fs_load_mgh([pn filesep hemis '.' mghimg '.mgh']);
D = [squeeze(D(:,:,:,1)), squeeze(D(:,:,:,2)), squeeze(D(:,:,:,3))];
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
Srf.Structural = strimg;
Srf.Functional = mghimg;
Srf.Hemisphere = hemis;
Srf.Cortex_Steps = NaN;
Srf.Vertices = V0;
Srf.Curvature = C';
Srf.Pial = P;
Srf.Inflated = I;
Srf.Sphere = S;
Srf.Faces = F;
Srf.Normals = N;
Srf.Curvature = C';
Srf.Area = A';
Srf.Thickness = T';
% Add map data
Srf.Data = NaN(4, size(V0,1));
Srf.Data(1,:) = zeros(1,size(V0,1));
Srf.Data(1,rois~=0) = ones(1, sum(rois~=0)); % Label predicted maps as R^2 = 1
Srf.Data(2,:) = x0';
Srf.Data(3,:) = y0';
Srf.Data(4,:) = rois';
Srf.Rule = 'X';
Srf.Values = valstrs;

%% Save surface data
save([hemis '_' mghimg '.mat'], 'Srf', '-v7.3');
disp(['Saved ' hemis '_' mghimg '.mat']);
samsrf_anatomy_srf([hemis '_' mghimg]);
new_line;
% Save ROI labels
mkdir('ROIs_Benson');
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V1'], 4, find(abs(rois) == 1));
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V2'], 4, find(abs(rois) == 2));
samsrf_srf2label(Srf, ['ROIs_Benson' filesep hemis '_V3'], 4, find(abs(rois) == 3));
