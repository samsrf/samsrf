function OutFile = samsrf_fit_cf(Model, SrfFile, Roi)
%
% OutFile = samsrf_fit_cf(Model, SrfFile, [Roi = ''])
%
% Fits a CF model using a fast grid search fitting procedure.
% For each vertex in the region of interest Roi, determines the best CF parameters
% for which the time series in the seed region Model.SeedRoi can predict its response.
%
%   Model:        Contains the parameters defining the analysis.
%
%   SrfFile:      Data to be analysed. You can provide a Srf structure directly or define the 
%                 Srf data filename is a char array without .mat extension (e.g. 'lh_Bars1').
%
%                 NOTE: Unlike in previous SamSrf versions, you cannot provide a file list 
%                       for concatenation! If you want concatenated runs, do this using the 
%                       surface projection tools available and provide the input as Srf struct.                     
%
%   Roi:          ROI label to restrict the analysis (default = '') 
%                     Optional, but without this the analysis can take forever.
%
% Returns the name of the map file it saved.
%
% 04/09/2024 - Instead of a list of files, you can now specify a Srf as input (DSS)
% 13/09/2024 - Removed option to provife data file list for concatenation (DSS)
%              Fixed help descriptions (DSS)
% 15/09/2024 - Fixed bug with undefined output filename when providing Srf data (DSS)
% 07/10/2025 - Models can now be provided as JSON files (DSS)
%

%% Load Model as JSON?
if ischar(Model) || isstring(Model)
    try
        Model = LoadModelJson(Model);
    catch
        samsrf_error(['Cannot load model file ' Model '.json!']);
    end
end

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end
% Define search space name
SearchspaceFile = ['src_' Model.Name '.mat'];  

%% Default model parameters
Model = samsrf_model_defaults('samsrf_fit_cf', Model);

%% Start time of analysis
t0 = tic; samsrf_newline;  
samsrf_disp('*** SamSrf CF model fitting ***');
[vn, vd] = samsrf_version;
samsrf_disp([' Version ' num2str(vn) ' - ' vd]);
samsrf_newline;
samsrf_disp('Current working directory:');
samsrf_disp([' ' pwd]);
samsrf_newline;

%% Load images 
% Were data provided directly?
if isstruct(SrfFile)
    samsrf_disp('Surface data input...');
    Srf = SrfFile; % Srf data provided directly
    clear SrfFile % To avoid duplicating massive variable
else
    samsrf_disp('Reading data file...')
    load([pwd filesep SrfFile]);
    Srf = samsrf_expand_srf(Srf);
    samsrf_disp([' Loading ' SrfFile ': ' num2str(size(Srf.Vertices,1)) ' vertices & ' num2str(size(Srf.Data,1)) ' volumes']);
end
Tc = Srf.Data; % Time course matrix 
OutFile = [Srf.Hemisphere '_' Model.Name];

%% Load ROI mask
if isempty(Roi)
    samsrf_newline; samsrf_disp('Using all vertices (''tis gonna take forever!)...');
    mver = 1:size(Srf.Vertices,1);
else
    samsrf_newline; samsrf_disp('Reading ROI mask...')
    mver = samsrf_loadlabel(Roi);
    samsrf_disp([' Loading ' Roi ': ' num2str(size(mver,1)) ' vertices']);
end

%% Correct by global mean signal?
if Model.Global_Signal_Correction
    samsrf_newline; samsrf_disp('Applying global signal correction...');
    Srf.Data = Tc;
    Srf = samsrf_removenoise(Srf, nanmean(Srf.Data,2), mver);
    Tc = Srf.Data;
    Srf.Data = [];
end
samsrf_newline; 

%% Limit data due to noise ceiling?
if isfield(Srf, 'Noise_Ceiling')
    if Model.Noise_Ceiling_Threshold > 0
        mver = mver(Srf.Noise_Ceiling(mver) > Model.Noise_Ceiling_Threshold);
        samsrf_disp(['Limiting analysis to ' num2str(length(mver)) ' vertices above noise ceiling ' num2str(Model.Noise_Ceiling_Threshold)]);
        samsrf_newline;
    end
end

%% Load seed ROI 
samsrf_disp(['Loading seed ROI: ' Model.SeedRoi]);
svx = samsrf_loadlabel(Model.SeedRoi);
Srf.SeedVx = svx; % Store for posterity

%% Load template map
samsrf_disp(['Loading template map: ' Model.Template]);
Temp = load(EnsurePath(Model.Template));
Temp.Srf = samsrf_expand_srf(Temp.Srf);
samsrf_newline;

%% Add version number
Srf.Version = samsrf_version;

%% Smooth for coarse fit?
if Model.Smoothing > 0
    Srf = samsrf_smooth_sphere(Srf, Model.Smoothing, Roi, 0); % Smooth within ROI
    Srf = rmfield(Srf, 'Raw_Data'); % Remove raw data field
    Tc = Srf.Data; % Put smoothed data back into time course variable
end

%% Generate prediction matrix
if ~exist([pwd filesep SearchspaceFile], 'file') 
    samsrf_disp('Generating predictions...');
    if isfield(Model, 'Polar')
        samsrf_disp(' Using polar angle patches as CFs');
        Polar = atan2(Temp.Srf.Data(3,svx), Temp.Srf.Data(2,svx)) / pi * 180;
        X = NaN(size(Tc,1), length(Model.Polar)-1); % Time courses
        S = NaN(1,length(Model.Polar)-1); % Patch polar angle
        for p = 1:length(Model.Polar)-1
            pvx = svx(Polar >= Model.Polar(p) & Polar < Model.Polar(p+1)); % Patch vertices in seed ROI
            pX = roi_1steigvar(Tc(:,pvx)); % 1st eigenvariate of patch
            S(p) = mean(Model.Polar(p:p+1)); % Polar angle of patch
            X(:,p) = pX; % Store extracted time course
        end
    elseif isfield(Model, 'Eccentricity')
        samsrf_disp(' Using eccentricity patches as CFs');
        Eccentricity = sqrt(Temp.Srf.Data(2,svx).^2 + Temp.Srf.Data(3,svx).^2);
        X = NaN(size(Tc,1), length(Model.Eccentricity)-1); % Time courses
        S = NaN(1,length(Model.Eccentricity)-1); % Patch eccentricity
        for p = 1:length(Model.Eccentricity)-1
            pvx = svx(Eccentricity >= Model.Eccentricity(p) & Eccentricity < Model.Eccentricity(p+1)); % Patch vertices in seed ROI
            pX = roi_1steigvar(Tc(:,pvx)); % 1st eigenvariate of patch
            S(p) = mean(Model.Eccentricity(p:p+1)); % Polar angle of patch
            X(:,p) = pX; % Store extracted time course
        end
    elseif isfield(Model, 'Sizes')
        samsrf_disp(' Using vertex-wise search space with circular CFs');
        [X,S] = cf_generate_searchspace(Srf, svx, Model.Sizes); % Store time courses for each seed vertex
    end
    save(SearchspaceFile, 'X', 'S', '-v7.3');
    t1 = toc(t0); 
    samsrf_disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
else
    samsrf_disp('Loading pre-defined predictions...');
    load([pwd filesep SearchspaceFile]);
    samsrf_disp([' Loading ' SearchspaceFile]);
    % Does number of grid points match model?
    if isfield(Model, 'Sizes')
        if size(S,2) ~= length(svx) * length(Model.Sizes)
            samsrf_error('Mismatch between saved search space and model definition!');
        end
    end
    % Does length of search parameter & prediction matrix match?
    if size(S,2) ~= size(X,2)
        % Shouldn't ever happen unless someone screwed with the search space file
        samsrf_error('Search space is corrupt! Mismatch between number of parameters & predictions!');
    end
end
samsrf_disp(['Using search space with ' num2str(size(S,2)) ' grid points.']);
samsrf_newline; 

%% Preprocess data
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = [];  % Clear data field

%% Coarse fit 
samsrf_disp('Coarse fitting...');
Xfits = zeros(size(Tc,1), length(mver));  % Matrix with predictions
Rimg = zeros(1, length(mver)); % R^2 map
Vimg = zeros(1, length(mver)); % Fitted vertex number/patch parameter map
if isfield(Model, 'Sizes')
    Simg = zeros(1, length(mver)); % Fitted size parameter map
    Pimg = zeros(2, length(mver)); % Template parameter maps
end

% Loop through mask vertices 
samsrf_disp(' Please stand by...');
parfor v = 1:length(mver)
    % Current vertex index
    vx = mver(v);
    % Find best prediction
    Pv = samsrf_georoi(vx, Model.Patch_Size, Srf.Vertices, Srf.Faces); % Circular patch
    Y = roi_1steigvar(Tc(:,Pv));
    if ~isempty(Y)
        R = corr(Y,X); % Mean corrected correlation 
        mR = max(R).^2; % Find best fit & square now
        R = R.^2; % Now turn others into R^2 too
        rx = find(R == mR,1); % Matrix position of best prediction
        if ~isempty(rx)
            % Store prediction
            Xfits(:,v) = X(:,rx);  % Best fitting prediction
            % Store parameters
            Rimg(1,v) = mR;  % Variance explained
            if isfield(Model, 'Sizes')
                Vimg(1,v) = S(1,rx); % Add vertex number
                Simg(1,v) = S(2,rx); % Add size parameter in geodesic steps
                Pimg(:,v) = Temp.Srf.Data(2:3, S(1,rx)); % Retrieve template coordinates 
            else
                Vimg(1,v) = S(rx); % Patch parameter (whichever it may be)
            end
        end
    end
end
t2 = toc(t0); 
samsrf_disp(['Coarse fitting completed in ' num2str(t2/60) ' minutes.']);
samsrf_newline;

%% Store coarse fit parameters
Srf.X = zeros(size(Tc)); % Matrix with best-fitting predictions
Srf.X(:,mver) = Xfits; % Fill in best-fitting predictions 

% Prepare surface structure
Srf.Functional = 'Connective Field Fast-Fit'; % pRF function name
if isfield(Model, 'Sizes')
    Data = [Rimg; Pimg; Simg; Vimg];
    Srf.Values = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'Vx'}; % Add parameter names
else
    Data = [Rimg; Vimg];
    Srf.Values = {'R^2'; 'Phase'}; % Add parameter name
end
Srf.Data = zeros(size(Data,1), size(Srf.Vertices,1));
Srf.Data(:,mver) = Data; % Add parameter maps into full data matrix
% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
    Srf = samsrf_normr2(Srf); % Calculate normalised R^2
end
samsrf_newline;

% Compress to save space
Srf = samsrf_compress_srf(Srf, mver);
% Save map files
OutFile = [OutFile '_Fwd'];
samsrf_disp('Saving CF fitting results...');
save(OutFile, 'Model', 'Srf', '-v7.3');
samsrf_disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60/60);
samsrf_newline; samsrf_disp(['Whole analysis completed in ' EndTime ' hours.']);
samsrf_disp('******************************************************************');
samsrf_newline; samsrf_newline;

