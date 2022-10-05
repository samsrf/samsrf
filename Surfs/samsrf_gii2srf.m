function samsrf_gii2srf(funimg, hemsurf, nrmls, avrgd, nsceil, anatpath)
% 
% samsrf_gii2srf(funimg, hemsurf, [nrmls=true, avrgd=true, nsceil=true, anatpath='../anatomy/'])
%
% Converts a GIfTI surface in GII format to a SamSrf surface file and saves it.
% You can use this instead of the MGH format when using FreeSurfer to project functional data to the surface and then analysing it in SamSrf.
%
% NOTE: This function requires SPM12 for GIfTI functionality.
%
%   funimg:     Name of functional GII files (without extension!)
%                 If this a cell array, files are averaged or concatenated (see avgconsep) 
%                 In that case you should probably normalise! (see nrmls)
%   hemsurf:    Hemisphere of surfaces (& folder if needed)
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
%
% 14/05/2022 - Written (DSS)
% 14/09/2022 - Edited help section to declare SPM12 dependency (DSS)  
% 05/10/2022 - Can now also detrend without z-normalisation (DSS)
%

%% Default parameters
if nargin < 3
    nrmls = true;
end
if nargin < 4
    avrgd = true;
end
if nargin < 5
    nsceil = true;
end
if nargin < 6
    anatpath = ['..' filesep 'anatomy' filesep];
end

% If input functional is a string, turn into cell array
if isa(funimg, 'char')
    funimg = {funimg};
end
% Trim functional GIfTI file names if neccesary
for fi = 1:length(funimg) 
    if strcmpi(funimg{fi}(end-3:end), '.gii')
        funimg{fi} = funimg{fi}(1:end-4);
    end
end
    
%% Load functional image
run1 = gifti([EnsurePath(funimg{1}) '.gii']);
run1 = run1.cdata';
fimg = NaN([size(run1) length(funimg)]);
for fi = 1:length(funimg)
    cur_run = gifti([funimg{fi} '.gii']);
    cur_run = cur_run.cdata';
    fimg(:,:,fi) = cur_run;
end

%% Load surface vertices
[V0 F] = fs_read_surf([hemsurf '.white']); % Grey-white surface
P = fs_read_surf([hemsurf '.pial']); % Pial surface
I = fs_read_surf([hemsurf '.inflated']); % Inflated surface
S = fs_read_surf([hemsurf '.sphere']); % Spherical surface
C = fs_read_curv([hemsurf '.curv']); % Cortical curvature 
A = fs_read_curv([hemsurf '.area']); % Cortical surface area
T = fs_read_curv([hemsurf '.thickness']); % Cortical thickness
N = P - V0; % Cortical vectors for each vertex 
[strimg,hemsurf] = fileparts(hemsurf);   % Remove folder from hemsurf

%% Create surface structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Structural = strimg;
Srf.Functional = funimg;
Srf.Hemisphere = hemsurf;
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
Srf.Data = fimg;
Srf.Rule = 'mri_vol2surf';

%% Normalise data
if nrmls
    % Only if more than 1 timepoint
    if size(Srf.Data, 1) > 1        
        % Normalise image for each run
        for fi = 1:length(funimg)
            Srf.Data(:,:,fi) = detrend(Srf.Data(:,:,fi)); % Linear detrending to remove drift
            if sign(nrmls) > 0
                % Z-normalisation
                Srf.Data(:,:,fi) = zscore(Srf.Data(:,:,fi)); 
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
% Srf.Data should now be only (volume, vertex)

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Save surface data
[p f e] = fileparts(funimg{1});
save(f, 'Srf', '-v7.3');
disp(['Saved ' f '.mat']);
samsrf_anatomy_srf(f, anatpath);
new_line;
