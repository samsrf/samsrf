function Srf = samsrf_vol2mat(funimg, roi, nrmls)
% Srf = samsrf_vol2mat(funimg, [mask, nrmls=true])
%
% Converts a NII functional volume to a SamSrf-compatible Matlab structure.
% The data file is prefixed 'vol_'
%
%   funimg:     name of functional NII file (without extension)
%                 This can be a cell array if more than one file is to be
%                 averaged. In that case you should probably normalise!
%   roi:        name of binary mask, in NII format (without extension)
%   nrmls:      if true, it will detrend & z-score the time series in each voxel.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 29/10/2020 - Added prefix for these files (DSS)
% 13/03/2022 - Now reports which default parameter file it's loading (DSS)
%              Changes error message when NII loading fails (DSS)
%

%% Default parameters
if exist('SamSrf_defaults.mat', 'file')
    disp(['Using defaults in: ' which('SamSrf_defaults.mat')]);
    load('SamSrf_defaults.mat');
end
if nargin < 2
    roi = [];
end
if nargin < 3
    nrmls = true;
end

% If input functional is a string, turn into cell array
if isa(funimg, 'char')
    funimg = {funimg};
end
% Trim file names if neccesary
for fi = 1:length(funimg) 
    if strcmpi(funimg{fi}(end-3:end), '.nii')
        funimg{fi} = funimg{fi}(1:end-4);
    end
end

% Trim ROI file name if necessary
if strcmpi(roi(end-3:end), '.nii')
    roi = roi(1:end-4);
end

%% Load functional image
if exist('spm', 'file') % Use SPM 
    fhdr = spm_vol([funimg{1} '.nii']);
    fimg = NaN([fhdr(1).dim length(fhdr) length(funimg)]);
    for fi = 1:length(funimg)
        fhdr = spm_vol([funimg{fi} '.nii']);
        fimg(:,:,:,:,fi) = spm_read_vols(fhdr);
    end
else 
    error('Sorry but I need SPM or nifti-patch to load NII files :(');
end

%% Load Roi
if ~isempty(roi)
    % Use SPM
    mhdr = spm_vol([roi '.nii']);
    mimg = spm_read_vols(mhdr);
    
    % Enforce binary mask
    mimg = logical(mimg);
end

%% Create SamSrf structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Functional = funimg;
Srf.Hemisphere = 'vol';
Srf.NiiHeader = fhdr;
Srf.Roi = [];

%% Indexing
% Size of the input matrix
Dims = size(fimg);
if length(Dims) < 5
    Dims(5) = 1;
end

% Reshape functional data to time x voxel x run
y = reshape(fimg, numel(fimg(:,:,:,1,1)), Dims(4), Dims(5));
Srf.Data = permute(y, [2 1 3]);

% Voxel order
% Rows -> Columns -> Slice
[xi, yi, zi] = meshgrid(1:Dims(1), 1:Dims(2), 1:Dims(3));

% Restrict to mask
if ~isempty(roi)
    % Reshape mask to matching 1 x voxel
    m = reshape(mimg, numel(mimg), 1);
    Srf.Roi = m;
    
    % Trim data, voxel indices
    Srf.Data = Srf.Data(:, m, :);
end

% Voxel dimensions
Srf.VoxDim = Dims(1:3);

% Fake vertices
Srf.Vertices = 1:size(Srf.Data, 2);
Srf.Vertices = Srf.Vertices';

%% Normalise
if nrmls
    if size(Srf.Data, 1) > 1       
        % Normalise each run
        for fi = 1:length(funimg)            
            Srf.Data(:,:,fi) = detrend(Srf.Data(:,:,fi)); % Linear detrend to remove drift
            Srf.Data(:,:,fi) = zscore(Srf.Data(:,:,fi)); % Normalise timeseries to z-score
        end
    else
        disp('Only one volume so won''t perform normalization!');
    end
end

%% Average separate runs 
if length(funimg) > 1
    Srf.Data = nanmean(Srf.Data, 3);
end
% Srf.Data: (time, voxel)
Srf.Data = squeeze(Srf.Data);

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Save data
[~, f, ~] = fileparts(funimg{1});
save(['vol_' f], 'Srf', '-v7.3');
disp(['Saved vol_' f '.mat']);
new_line;
