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
%   nrmls:      If true, it will detrend & normalise the time series in each vertex.
%                 If positive, it will use z-normalisation.
%                 If negative, it will only detrend but not z-normalise.
%   avrgd:      If true, runs will be averaged into one SamSrf file (default).
%               If false, runs will be concatenated into one SamSrf file.
%   nsceil:     If true, calculates the noise ceiling by splitting data into odd and even runs.
%                 The noise ceiling is stored in the vector Srf.Noise_Ceiling.
%                 This option only works when averaging runs - otherwise it is ignored 
%                   (this may change in future versions)
%   anatpath:   Defines path where anatomy meshes are stored. Defaults to '../anatomy/'
%                 If this is empty, the anatomy is not split off!
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
% 13/03/2022 - Ensures now that random files aren't loaded from path (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 05/10/2022 - Can now also detrend without z-normalisation (DSS)
% 15/10/2022 - Default normalisation is now 1 instead of true (DSS)
% 29/06/2023 - Added conversion to 32 bit (single) data (DSS)
% 19/12/2023 - Now has option not to split off anatomy (DSS)  
%

%% Default parameters
if nargin < 4
    ctxsteps = 0.5;
end
if nargin < 5
    rule = 'Mean';
end
if nargin < 6
    nrmls = 1;
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
    hdr = spm_vol([EnsurePath(strimg) '.nii']);
    % Origin in the actual structural
    nii_orig = hdr.mat(1:3,4);
    % Origin in Freesurfer space (1/2 dimensions)
    fs_orig = hdr.dim' / 2;
    fs_orig = fs_orig([3 1 2]) .* sign(nii_orig); 
else % Sadly no way to load NIIs
    samsrf_error('Sorry but I need SPM to load NII files :(');
end
    
%% Load functional image
fhdr = spm_vol([EnsurePath(funimg{1}) '.nii']);
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
try
    C = fs_read_curv([hemsurf '.curv']); % Cortical curvature 
    A = fs_read_curv([hemsurf '.area']); % Cortical surface area
    T = fs_read_curv([hemsurf '.thickness']); % Cortical thickness
catch
    % If no binary files available, hope there are ASC files, and use those
    C = Read_FreeSurfer([hemsurf '.curv.asc']); % Cortical curvature 
    A = Read_FreeSurfer([hemsurf '.area.asc']); % Cortical surface area
    T = Read_FreeSurfer([hemsurf '.thickness.asc']); % Cortical thickness
    C = C(:,5);
    A = A(:,5);
    T = T(:,5);
end
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
samsrf_disp('Running surface projection...');
samsrf_newline;
cs = 0;
for cl = ctxsteps
    samsrf_disp(['Cortical sampling step: ' num2str(cl)]);
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
    samsrf_progbar(0);
    for i = 1:size(tV,1) 
        for fi = 1:length(funimg)
            if tV(i,1)>0 && tV(i,2)>0 && tV(i,3)>0 && tV(i,1)<fhdr(1).dim(1) && tV(i,2)<fhdr(1).dim(2) && tV(i,3)<fhdr(1).dim(3) 
                Srf.Data(cs, :, i, fi) = squeeze(fimg(tV(i,1), tV(i,2), tV(i,3), :, fi))'; 
            else
                Srf.Data(cs, :, i, fi) = NaN; 
            end
        end
        % Progress report
        samsrf_progbar(i/size(tV,1));            
    end
    samsrf_newline;
end

%% Calculate one value per vertex
if strcmpi(rule, 'Mean')
    samsrf_disp('Using mean of steps through cortex');
    Srf.Data = squeeze(nanmean(Srf.Data,1));
elseif strcmpi(rule, 'Maximum')
    samsrf_disp('Using maximum of steps through cortex');
    Srf.Data = squeeze(nanmax(Srf.Data,[],1));
elseif strcmpi(rule, 'Minimum')
    samsrf_disp('Using minimum of steps through cortex');
    Srf.Data = squeeze(nanmin(Srf.Data,[],1));
elseif strcmpi(rule, 'Median')
    samsrf_disp('Using median of steps through cortex');
    Srf.Data = squeeze(nanmedian(Srf.Data,1));
elseif strcmpi(rule, 'Sum')
    samsrf_disp('Using sum of steps through cortex');
    Srf.Data = squeeze(nansum(Srf.Data,1));
elseif strcmpi(rule, 'Geomean')
    samsrf_disp('Using geometric mean of steps through cortex');
    Srf.Data = squeeze(exp(nanmean(log(Srf.Data,1))));
else
    % Unless this happens
    samsrf_disp('Retaining individual steps through cortex');
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
                if sign(nrmls) > 0
                    % Z-normalisation
                    Srf.Data(:,:,fi,cl) = zscore(Srf.Data(:,:,fi,cl)); 
                end
            end
        end
    else
        samsrf_disp('Only one volume so won''t perform normalization!');
    end
end

%% Combine separate runs into one file
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
% Srf.Data: (volume, vertex, [ctxsteps, if >1])
Srf.Data = squeeze(Srf.Data);

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);

%% Save surface data
[p f e] = fileparts(funimg{1});
save([hemsurf '_' f], 'Srf', '-v7.3');
samsrf_disp(['Saved ' hemsurf '_' f '.mat']);
samsrf_anatomy_srf([hemsurf '_' f], anatpath);
samsrf_newline;

