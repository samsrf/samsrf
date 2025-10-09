function Srf = samsrf_gii2srf(funimg, hemsurf, nrmls, avrgd, nsceil, anatpath)
% 
% Srf = samsrf_gii2srf(funimg, [hemsurf='', nrmls=true, avrgd=true, nsceil=true, anatpath='../anatomy/'])
%
% Converts a GIfTI surface in GII format to a SamSrf surface structure. If no output argument is provided, it saves the Srf as a file.
% You can use this instead of the MGH format when using FreeSurfer to project functional data to the surface and then analysing it in SamSrf.
%
% IMPORTANT: This function requires SPM12 for GIfTI functionality!
%
%   funimg:     Name of functional GII files (without extension!)
%                 If this a cell array, files are averaged or concatenated (see avrgd) 
%                 In that case you should probably normalise! (see nrmls)
%				  You can instead provide a wildcard to search for files (e.g. 'lh_*') 
%   hemsurf:    Hemisphere of surfaces & folder (e.g. '../surf/lh')
%				  This will produce a standard file for use with SamSrf.
%				Alternatively, you can provide the GII file for -one- surface so they can be saved together.
%				  This will save that surface in Srf.Vertices so you wont know which model this is. 
%				  It will contain dummy values for the other surface meshes so they cannot be used!
%				  For such files anatpath (see below) is automatically set to ''
%				Providing anatomical info is -not- necessary - without this it will create a Srf without any surface data.
%				  Such files are probably fine for most purposes, for example they can be saved out as FreeView-readable GIIs.
%				  But obviously, such a file cannot be rendered as a surface model in DisplayMaps or samsrf_surf.
%				  In this case, the Srf.Hemisphere prefix will be set to 'srf'.
%				  This is in fact the default setting so hemsurf defaults to ''
%				  For such files anatpath (see below) is automatically set to ''
%   nrmls:      If true, it will detrend & normalise the time series in each vertex.
%                 If positive, it will use z-normalisation.
%                 If negative, it will only detrend but not z-normalise.
%   avrgd:      If true, runs will be averaged into one SamSrf structure (default).
%               If false, runs will be concatenated into one SamSrf structure.
%   nsceil:     If true, calculates the noise ceiling by splitting data into odd and even runs.
%                 The noise ceiling is stored in the vector Srf.Noise_Ceiling.
%                 This option only works when averaging runs - otherwise it is ignored 
%                   (this may change in future versions)
%   anatpath:   Defines relative path where anatomy meshes are stored. Defaults to '../anatomy/'
%                 If this is empty, the anatomy is not split off!
%				  Note that in standalone mode anatomy is never split off!
%
% 04/09/2024 - Clarifications in help section (DSS)
%              Now allows returning an Srf instead of saving it. (DSS)
% 19/09/2024 - Srf.Structural is now the subject's surf folder (DSS)
% 10/06/2025 - Added check if SPM is on path for GII functionality (DSS)
% 07/10/2025 - Now accepts wildcard input for GII files (DSS)
% 08/10/2025 - Adapted for compiled command line analysis (DSS)
% 09/10/2025 - Can now read in GII files without providing anatomical info (DSS)
%			   Default for hemsurf is now also '' to not include anatomical info (DSS) 		   
%

%% Default parameters
if nargin < 2
	hemsurf = '';
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
if nargin < 6
    anatpath = ['..' filesep 'anatomy' filesep];
end

% If input functional is a string, turn into cell array
if isa(funimg, 'char')
    if contains(funimg, '*')
    	funimg = dir([funimg '.gii']);
        funimg = {funimg.name}';
        samsrf_disp('Found GII files:');
        samsrf_disp(funimg);
    else
    	funimg = {funimg};
    end
end
% Trim functional GIfTI file names if neccesary
for fi = 1:length(funimg) 
    if strcmpi(funimg{fi}(end-3:end), '.gii')
        funimg{fi} = funimg{fi}(1:end-4);
    end
end
    
%% Load functional image
if exist('spm', 'file') % Use SPM 
    run1 = gifti([EnsurePath(funimg{1}) '.gii']);
    run1 = run1.cdata';
    fimg = NaN([size(run1) length(funimg)]);
    for fi = 1:length(funimg)
        cur_run = gifti([funimg{fi} '.gii']);
        cur_run = cur_run.cdata';
        fimg(:,:,fi) = cur_run;
    end
else
    samsrf_error('Requires SPM for reading GII files!');
end

%% Load surface vertices
if ~isempty(hemsurf)
    if contains(lower(hemsurf), '.gii')
        % Read surface data from GII file
        AnatSurf = gifti(EnsurePath(hemsurf));
        strimg = hemsurf;
        hemsurf = 'srf';
	    anatpath = ''; % Don't split off non-existing anatomy info
        V0 = double(AnatSurf.vertices);
        F = double(AnatSurf.faces);
    else
        % Read surface data from FreeSurfer/surf folder
        [V0 F] = fs_read_surf([hemsurf '.white']); % Grey-white surface
        P = fs_read_surf([hemsurf '.pial']); % Pial surface
        I = fs_read_surf([hemsurf '.inflated']); % Inflated surface
        S = fs_read_surf([hemsurf '.sphere']); % Spherical surface
        C = fs_read_curv([hemsurf '.curv']); % Cortical curvature 
        A = fs_read_curv([hemsurf '.area']); % Cortical surface area
        T = fs_read_curv([hemsurf '.thickness']); % Cortical thickness
        N = P - V0; % Cortical vectors for each vertex 
        [strimg,hemsurf] = fileparts(hemsurf);   % Remove folder from hemsurf
    end
else
    % No surface data provided
    samsrf_disp('Reading in GII file without any anatomical info');
    strimg = 'n/a'; % No structural provided as name could be anything
    hemsurf = 'srf'; % Indicates this is a surface without anatomy info
    anatpath = ''; % Don't split off non-existing anatomy info
    V0 = 1:size(fimg,2); % Vertices provided by indeces
    V0 = V0'; % Ensure column vector
    F = NaN; % No faces here
end

%% Create surface structure
Srf = struct;
Srf.Version = samsrf_version;
Srf.Structural = strimg;
Srf.Functional = funimg;
Srf.Hemisphere = hemsurf;
Srf.Cortex_Steps = NaN;
Srf.Vertices = V0;
Srf.Faces = F;
if ~strcmpi(hemsurf, 'srf')
    Srf.Pial = P;
    Srf.Inflated = I;
    Srf.Sphere = S;
    Srf.Normals = N;
    Srf.Curvature = C';
    Srf.Area = A';
    Srf.Thickness = T';
end
Srf.Data = fimg;
Srf.Rule = 'gii';

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
% Srf.Data should now be only (volume, vertex)

%% Add values field
Srf.Values = {};
for r = 1:size(Srf.Data,1)
    Srf.Values{r,1} = ['Volume #' num2str(r)];
end

%% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);

%% Save surface data?
if nargout == 0 
    [p f e] = fileparts(funimg{1});
    save([f '.mat'], 'Srf', '-v7.3');
    samsrf_disp(['Saved ' f '.mat']);
    samsrf_anatomy_srf(f, anatpath);
    samsrf_newline;
end