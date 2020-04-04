function samsrf_label2nii(labelfile, funimg, strimg, hemsurf, ctxsteps, scalar)
% 
% samsrf_label2nii(labelfile, funimg, strimg, hemsurf, [ctxsteps], [scalar])
%
% Converts a freesurfer label into NII binary mask.
%
%   labelfile:  name of label file (without extension)
%   funimg:     name of functional NII file (without extension)
%   strimg:     name of structural NII file (without extension)
%   hemsurf:    hemisphere of surfaces (& folder if needed)
%   ctxsteps:   vector with proportional steps through the grey matter 
%                   to check which functional voxel contains a vertex
%                   (Optional, defaults to [0.5] but [0.1] is closer to what FreeSurfer does)
%   scalar:     save scalar values stored in labelfile (true/false)
%
% 09/08/2018 - SamSrf 6 version (DSS) 
% 21/02/2020 - Added Matlab-native NIfTI support (IA)
% 11/03/2020 - Removed MatLab-native NIfTI support again because too complex (DSS)
% 01/04/2020 - IMPORTANT UPDATE: Removed the need (I hope) for Coregistration.txt!!! 
%              Fixed bug with native reader still being present for functional image *sigh* (DSS)
%              Corrected the help section in this function (DSS)
% 02/04/2020 - Removed inconsequential line from code (DSS)
%

if nargin < 5
    ctxsteps = 0.5;
end
if nargin < 6
    scalar = false;
end

%% Load structural header 
if exist('spm', 'file')
    hdr = spm_vol([strimg '.nii']);
    hdr = hdr(1);
    % Origin in the actual structural
    nii_orig = hdr.mat(1:3,4);
    % Origin in Freesurfer space (1/2 dimensions)
    fs_orig = hdr.dim' / 2;
    fs_orig = fs_orig([3 1 2]) .* sign(nii_orig);
else
    error('Sorry but I need SPM to load NII files :(');
end
    
%% Load functional image
fhdr = spm_vol([funimg '.nii']);
fhdr = fhdr(1);
% Empty functional image
fimg = zeros(fhdr.dim);
% Image matrix dimensions
funcdim = fhdr.dim(1:3);

%% Transformation matrices
fs_mat = [hdr.mat(1:3,1:3)' hdr.dim'/2; 0 0 0 1]; % FreeSurfer matrix
mat = fhdr.mat; % Functional matrix

%% Load surface vertices
V0 = fs_read_surf([hemsurf '.white']); % Grey-white surface
P = fs_read_surf([hemsurf '.pial']); % Pial surface
N = P - V0; % Cortical vectors for each vertex 

% Load ROI label
labeldata = Read_FreeSurfer([labelfile '.label']);
R = samsrf_loadlabel(labelfile); % Should be equal to labeldata(:,1)+1
V0 = V0(R,:);
N = N(R,:);

% What value to store
if scalar
    M = double(labeldata(:,5));
else
    M = ones([size(labeldata, 1), 1]);
end

%% Transform the vertices 
for cl = ctxsteps
    % Step through cortex layers
    V = V0 + N*cl;
    
    % Transformation into voxel space
    tV = fs_mat * [V'; ones(1,size(V,1))]; % Transform FreeSurfer into T1 voxel space
    % Transform into functional space?
    if ~strcmpi(funimg, strimg)
        tV = hdr.mat * tV; % Transform with T1 matrix
        tV = mat \ tV; % Now transform with functional matrix
    end
    tV = round(tV);
    tV = tV(1:3,:)';
    tV(tV(:,1) <= 0,:) = [];
    tV(tV(:,2) <= 0,:) = [];
    tV(tV(:,3) <= 0,:) = [];
    tV(tV(:,1) > funcdim(1),:) = [];
    tV(tV(:,2) > funcdim(2),:) = [];
    tV(tV(:,3) > funcdim(3),:) = [];

    % Mark voxels for each vertex
    for i = 1:size(tV,1) 
        fimg(tV(i,1), tV(i,2), tV(i,3)) = M(i); 
    end
end

%% Save new image 
% Enforce data type
if fhdr.dt(1) < 16
    fhdr.dt = [16 0]; % Enforce float32 type
end
% Write
fhdr.fname = [labelfile '.nii'];
spm_write_vol(fhdr, fimg);  
% Message
disp(['Saved ' labelfile '.nii.']);

