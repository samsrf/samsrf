function Srf = samsrf_vol2mat(funimg, roi, nrmls, avrgd, nsceil)
%
% Srf = samsrf_vol2mat(funimg, [roi, nrmls=1, avrgd=true, nsceil=true])
%
% Converts a NII functional volume to a SamSrf-compatible Matlab structure.
% If no output argument is provided, the Srf is saved as a file prefixed 'vol_'
%
% NOTE: This is really a hack trying to make volumetric data work with the
%       standard Srf variable structure. Not all functionality is supported.
%
%   funimg:     Name of functional NII file (without extension)
%                 This can be a cell array if more than one file is to be
%                 averaged. In that case you should probably normalise!
%   roi:        Name of binary mask, in NII format (without extension)
%   nrmls:      If true, it will detrend & normalise the time series in each vertex.
%                 If positive, it will use z-normalisation.
%                 If negative, it will only detrend but not z-normalise.
%   avrgd:      If true, runs will be averaged into one SamSrf structure (default).
%               If false, runs will be concatenated into one SamSrf structure.
%   nsceil:     If true, calculates the noise ceiling by splitting data into odd and even runs.
%                 The noise ceiling is stored in the vector Srf.Noise_Ceiling.
%                 This option only works when averaging runs - otherwise it is ignored 
%                   (this may change in future versions)
%
% 14/09/2024 - Saving the file is now optional (DSS)  
% 08/10/2025 - Now accepts wildcard input for GII files (DSS)
% 			   Adapted for compiled command line analysis (DSS)
%

%% Default parameters
if nargin < 2
    roi = [];
end
if nargin < 3
    nrmls = 1;
end
if nargin < 4
    avrgd = true;
end
if nargin < 5
    nsceil = true;
end

% If input functional is a string, turn into cell array
if isa(funimg, 'char')
    if contains(funimg, '*')
    	funimg = dir([funimg '.nii']);
        funimg = {funimg.name}';
        samsrf_disp('Found NII files:');
        samsrf_disp(funimg);
    else
    	funimg = {funimg};
    end
end
% Trim file names if neccesary
for fi = 1:length(funimg) 
    if strcmpi(funimg{fi}(end-3:end), '.nii')
        funimg{fi} = funimg{fi}(1:end-4);
    end
end

% Trim ROI file name if necessary
if ~isempty(roi)
    if strcmpi(roi(end-3:end), '.nii')
        roi = roi(1:end-4);
    end
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
    samsrf_error('Sorry but I need SPM or NIfTI-patch to load NII files :(');
end

%% Load Roi
if ~isempty(roi)
    % Use SPM
    mhdr = spm_vol([roi '.nii']);
    mimg = spm_read_vols(mhdr);
    
    % Enforce binary mask
    mimg = logical(mimg);
else
    samsrf_disp('WARNING: No ROI specified! This is probably unwise...');
end

%% Create SamSrf structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Functional = funimg;
Srf.Hemisphere = 'vol';
Srf.NiiHeader = fhdr(1);
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
            if sign(nrmls) > 0
                % Z-normalisation
                Srf.Data(:,:,fi) = zscore(Srf.Data(:,:,fi)); 
            end
        end
    else
        samsrf_disp('Only one volume so won''t perform normalization!');
    end
end

%% Average separate runs 
if length(funimg) > 1
    if avrgd
        % Average runs 
        samsrf_disp('Averaging runs...');
        % Calculate noise ceiling?
        if nsceil && size(Srf.Data, 3) > 1
            samsrf_disp('Calculating noise ceiling...');
            OddRuns = nanmean(Srf.Data(:,:,1:2:end), 3);
            EvenRuns = nanmean(Srf.Data(:,:,2:2:end), 3);
            % Loop thru vertices
            Srf.Noise_Ceiling = NaN(1, size(Srf.Data,2));
            for v = 1:size(Srf.Data, 2)
                Rho_xxp = corr(OddRuns(:,v), EvenRuns(:,v)); % Correlation between odd & even runs
                Srf.Noise_Ceiling(v) = (2*Rho_xxp) / (1+Rho_xxp); % Spearman-Brown prediction formula
                % For observable correlation we would need to take 
                % square root so without the square root this is R^2!
            end
        end
        % Calculate mean across runs
        Srf.Data = nanmean(Srf.Data, 3);
    else
        % Concatenate runs
        samsrf_disp('Concatenating runs...');
        if size(Srf.Data, 3) > 1 % If individual runs contained only one row this is unnecessary because they have already been squeezed
            nSrf = Srf;
            nSrf.Data = [];
            for fi = 1:length(funimg)
                nSrf.Data = [nSrf.Data; Srf.Data(:,:,fi)];
            end
            Srf = nSrf;
            clear nSrf
        end
    end
end
% Srf.Data: (time, voxel)
Srf.Data = squeeze(Srf.Data);

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);

%% Save data?
if nargout == 0
    [~, f, ~] = fileparts(funimg{1});
    save(['vol_' f], 'Srf', '-v7.3');
    samsrf_disp(['Saved vol_' f '.mat']);
    samsrf_newline;
end
