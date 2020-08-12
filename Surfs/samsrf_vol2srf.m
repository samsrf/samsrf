function samsrf_vol2srf(funimg, strimg, hemsurf, ctxsteps, rule, nrmls, avrgd, nsceil, anatpath)
% 
% samsrf_vol2srf(funimg, strimg, hemsurf, [ctxsteps=0.5, rule='Mean', nrmls=true, avrgd=true, nsceil=true, anatpath='../anatomy/'])
%
% Converts NII functional files to a SamSrf surface file and saves it.
%
%   funimg:     Name of functional NII files (without extension!)
%                 If this a cell array, files are averaged or concatenated (see avgconsep) 
%                 In that case you should probably normalise! (see nrmls)
%   strimg:     Name of structural NII file (without extension!)
%   hemsurf:    Hemisphere of surfaces (& folder if needed)
%   ctxsteps:   Vector with proportional steps through the grey matter 
%                   to check which functional voxel contains a vertex
%                   (Optional, defaults to [0.5])
% 	rule:       Method how multiple voxels per vertex are interpreted:
%                   'Mean':     Arithmetic mean (default) 
%                   'Maximum':  Maximum
%                   'Minimum':  Minimum
%                   'Median':   Median
%                   'Sum':      Sum
%                   'Geomean':  Geometric mean
%               (anything else will retain the individual steps - this is not advised for large data sets!)
%   nrmls:      If true, it will detrend & z-score the time series in each vertex.
%   avrgd:      If true, runs will be averaged into one SamSrf file (default).
%               If false, runs will be concatenated into one SamSrf file.
%   nsceil:     If true, calculates the noise ceiling by splitting data into odd and even runs.
%                 The noise ceiling is stored in the vector Srf.Noise_Ceiling.
%                 This option only works when averaging runs - otherwise it is ignored 
%                   (this may change in future versions)
%   anatpath:   Defines path where anatomy meshes are stored. Defaults to '../anatomy/'
%
% If multiple ctxsteps are requested and no collapsing rule is specified, 
% final output will consists of timepoints x vertices x ctxsteps.
% Note that normalisation is applied independently to each cortex step.
%
% IMPORTANT DIFFERENCE FROM PRIOR VERSIONS OF SAMSRF:
% The anatomical meshes are automatically split from the functional data
% and are then stored in ../anatomy/ (see samsrf_anatomy_srf) unless you 
% change the default path. If this file with anatomical meshes already exists 
% this step is skipped. 
% IT IS YOUR RESPONSIBILITY TO CHECK YOU'RE USING THE RIGHT ANATOMICAL DATA!
%
% 29/06/2020 - SamSrf 7 version (DSS)
% 22/07/2020 - Added progress reports but still no parallel processing (DSS)
%

%% Default parameters
if nargin < 4
    ctxsteps = 0.5;
end
if nargin < 5
    rule = 'Mean';
end
if nargin < 6
    nrmls = true;
end
if nargin < 7
    avrgd = true;
end
if nargin < 8
    nsceil = true;
end
if nargin < 9
    anatpath = ['..' filesep 'anatomy' filesep];
end

% If input functional is a string, turn into cell array
if isa(funimg, 'char')
    funimg = {funimg};
end
% Trim functional file names if neccesary
for fi = 1:length(funimg) 
    if strcmpi(funimg{fi}(end-3:end), '.nii')
        funimg{fi} = funimg{fi}(1:end-4);
    end
end

% Trim T1 file name if neccesary
if strcmpi(strimg(end-3:end), '.nii') 
    strimg = strimg(1:end-4);
end

%% Load structural header
if exist('spm', 'file') % Use SPM
    hdr = spm_vol([strimg '.nii']);
    % Origin in the actual structural
    nii_orig = hdr.mat(1:3,4);
    % Origin in Freesurfer space (1/2 dimensions)
    fs_orig = hdr.dim' / 2;
    fs_orig = fs_orig([3 1 2]) .* sign(nii_orig); 
else % Sadly no way to load NIIs
    error('Sorry but I need SPM to load NII files :(');
end
    
%% Load functional image
fhdr = spm_vol([funimg{1} '.nii']);
fimg = NaN([fhdr(1).dim length(fhdr) length(funimg)]);
for fi = 1:length(funimg)
    fhdr = spm_vol([funimg{fi} '.nii']);
    fimg(:,:,:,:,fi) = spm_read_vols(fhdr);
end
fhdr(1).dim(4) = length(fhdr);

%% Transformation matrices
fs_mat = [hdr.mat(1:3,1:3)' hdr.dim'/2; 0 0 0 1]; % FreeSurfer matrix
mat = fhdr.mat; % Functional matrix

%% Load surface vertices
[V0 F] = fs_read_surf([hemsurf '.white']); % Grey-white surface
P = fs_read_surf([hemsurf '.pial']); % Pial surface
I = fs_read_surf([hemsurf '.inflated']); % Inflated surface
S = fs_read_surf([hemsurf '.sphere']); % Spherical surface
C = fs_read_curv([hemsurf '.curv']); % Cortical curvature 
A = fs_read_curv([hemsurf '.area']); % Cortical surface area
T = fs_read_curv([hemsurf '.thickness']); % Cortical thickness
N = P - V0; % Cortical vectors for each vertex 
[~,hemsurf] = fileparts(hemsurf);   % Remove folder from hemsurf

%% Create surface structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Structural = strimg;
Srf.Functional = funimg;
Srf.Hemisphere = hemsurf;
Srf.Cortex_Steps = ctxsteps;
Srf.Vertices = V0;
Srf.Pial = P;
Srf.Inflated = I;
Srf.Sphere = S;
Srf.Faces = F;
Srf.Normals = N;
Srf.Curvature = C';
Srf.Area = A';
Srf.Thickness = T';
Srf.Data = NaN(length(ctxsteps), fhdr(1).dim(4), size(V0,1), length(funimg));
Srf.Rule = rule;

%% Transform the vertices 
disp('Running surface projection...');
cs = 0;
for cl = ctxsteps
    disp([' Cortical sampling step: ' num2str(cl)]);
    cs = cs + 1;
    % Step through cortex layers
    V = V0 + N*cl;
    
    % Transformation into voxel space
    tV = fs_mat * [V'; ones(1,size(V,1))]; % Transform FreeSurfer into T1 voxel space
    % Transform into functional space?
    if ~strcmpi(funimg{1}, strimg)
        tV = hdr.mat * tV; % Transform with T1 matrix
        tV = mat \ tV; % Now transform with functional matrix
    end
    tV = round(tV);
    tV = tV(1:3,:)';

    % Find voxels for each vertex
    for i = 1:size(tV,1) 
        for fi = 1:length(funimg)
            if tV(i,1)>0 && tV(i,2)>0 && tV(i,3)>0 && tV(i,1)<fhdr(1).dim(1) && tV(i,2)<fhdr(1).dim(2) && tV(i,3)<fhdr(1).dim(3) 
                Srf.Data(cs, :, i, fi) = squeeze(fimg(tV(i,1), tV(i,2), tV(i,3), :, fi))'; 
            else
                Srf.Data(cs, :, i, fi) = NaN; 
            end
        end
        % Progress report
        if mod(i,25000) == 0
            disp(['  ' num2str(round(i/size(tV,1)*100)) '% complete']);
        end
    end       
end

%% Calculate one value per vertex
if strcmpi(rule, 'Mean')
    disp('Using mean of steps through cortex');
    Srf.Data = squeeze(nanmean(Srf.Data,1));
elseif strcmpi(rule, 'Maximum')
    disp('Using maximum of steps through cortex');
    Srf.Data = squeeze(nanmax(Srf.Data,[],1));
elseif strcmpi(rule, 'Minimum')
    disp('Using minimum of steps through cortex');
    Srf.Data = squeeze(nanmin(Srf.Data,[],1));
elseif strcmpi(rule, 'Median')
    disp('Using median of steps through cortex');
    Srf.Data = squeeze(nanmedian(Srf.Data,1));
elseif strcmpi(rule, 'Sum')
    disp('Using sum of steps through cortex');
    Srf.Data = squeeze(nansum(Srf.Data,1));
elseif strcmpi(rule, 'Geomean')
    disp('Using geometric mean of steps through cortex');
    Srf.Data = squeeze(exp(nanmean(log(Srf.Data,1))));
else
    % Unless this happens
    disp('Retaining individual steps through cortex');
    Srf.Data = permute(Srf.Data, [2 3 4 1]);
end

% Far fewer columns than rows? Something went wrong above & vertices are in rows 
if size(Srf.Data, 2) < size(Srf.Data,1) 
    Srf.Data = Srf.Data';
end

%% Normalise data
if nrmls
    % Only if more than 1 timepoint
    if size(Srf.Data, 1) > 1

        % Number of cortical steps
        if ndims(Srf.Data) < 4
            cs = 1;
        else
            cs = length(ctxsteps);
        end
        
        % Normalise each image, at each cortical step
        for fi = 1:length(funimg)
            for cl = 1:cs
                Srf.Data(:,:,fi,cl) = detrend(Srf.Data(:,:,fi,cl)); % Linear detrending to remove drift
                Srf.Data(:,:,fi,cl) = zscore(Srf.Data(:,:,fi,cl)); % Normalize time series to z-score
            end
        end
    else
        disp('Only one volume so won''t perform normalization!');
    end
end

%% Combine separate runs into one file
if length(funimg) > 1
    if avrgd
        % Average runs 
        disp('Averaging runs...');
        % Calculate noise ceiling?
        if nsceil && size(Srf.Data, 3) > 1
            disp('Calculating noise ceiling...');
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
        disp('Concatenating runs...');
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
% Srf.Data: (volume, vertex, [ctxsteps, if >1])
Srf.Data = squeeze(Srf.Data);

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Save surface data
[p f e] = fileparts(funimg{1});
save([hemsurf '_' f], 'Srf', '-v7.3');
disp(['Saved ' hemsurf '_' f '.mat']);
samsrf_anatomy_srf([hemsurf '_' f], anatpath);
new_line;

