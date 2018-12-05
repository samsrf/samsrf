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
% A text file called Coregistration.txt must be present in the surface
% folder. If it isn't the transformation between native space and
% FreeSurfer may not work. At BUCNI or in the FIL this has not been an
% issue but it may be an issue in other centres and it will be a problem if
% your structural NII is not in standard 1mm^3 resolution.
%
% 09/08/2018 - SamSrf 6 version (DSS) 

if nargin < 5
    ctxsteps = 0.5;
end
if nargin < 6
    scalar = false;
end

%% If using the registration matrices from FreeSurfer
if exist([hemsurf(1:end-2) 'Coregistration.txt'], 'file')
    fs2nii = dlmread([hemsurf(1:end-2) 'Coregistration.txt']);
    Tmov = fs2nii(1:4,:);
    Reg = fs2nii(5:8,:);
    useRegDat = true;
else
    useRegDat = false;
end

%% Load structural header 
hdr = spm_vol([strimg '.nii']);
hdr = hdr(1);
% Origin in the actual structural
nii_orig = hdr.mat(1:3,4);
% Origin in Freesurfer space (1/2 dimensions)
fs_orig = hdr.dim' / 2;
fs_orig = fs_orig([3 1 2]) .* sign(nii_orig);

%% Load functional image
fhdr = spm_vol([funimg '.nii']);
fhdr = fhdr(1);
% Empty functional image
fimg = zeros(fhdr.dim);

%% Adjust transformation matrix
mov = nii_orig - fs_orig;
mat = fhdr.mat;
if useRegDat
    smat = hdr.mat;
else
    mat(1:3,4) = mat(1:3,4) - mov;
end

%% Load surface vertices
V0 = fs_read_surf([hemsurf '.white']); % Grey-white surface
P = fs_read_surf([hemsurf '.pial']); % Pial surface
N = P - V0; % Cortical vectors for each vertex 

% Load ROI label
labeldata = Read_FreeSurfer([labelfile '.label']);
V = labeldata(:,1)+1;
R = samsrf_loadlabel(labelfile);
V0 = V0(R,:);
N = N(R,:);

% What value to store
if scalar
    M = double(labeldata(:,5));
    if fhdr.dt(1) < 16
        fhdr.dt = [16 0]; % Enforce float32 type
    end
else
    M = ones([size(labeldata, 1), 1]);
end

%% Transform the vertices 
for cl = ctxsteps
    % Step through cortex layers
    V = V0 + N*cl;
    
    % Transformation into voxel space
    if useRegDat
        tV = Tmov \ Reg * [V'; ones(1,size(V,1))]; 
        if ~strcmpi(funimg, strimg)
            tV = smat * tV;
            tV = mat \ tV; 
        end
        tV = round(tV);
    else
        tV = round(mat \ [V'; ones(1,size(V,1))]); 
    end
    tV = tV(1:3,:)';
    tV(tV(:,1) <= 0,:) = [];
    tV(tV(:,2) <= 0,:) = [];
    tV(tV(:,3) <= 0,:) = [];
    tV(tV(:,1) > fhdr.dim(1),:) = [];
    tV(tV(:,2) > fhdr.dim(2),:) = [];
    tV(tV(:,3) > fhdr.dim(3),:) = [];

    % Mark voxels for each vertex
    for i = 1:size(tV,1) 
        fimg(tV(i,1), tV(i,2), tV(i,3)) = M(i); 
    end
end

%% Save new image 
fhdr.fname = [labelfile '.nii'];
spm_write_vol(fhdr, fimg);
disp(['Saved ' labelfile '.nii.']);
