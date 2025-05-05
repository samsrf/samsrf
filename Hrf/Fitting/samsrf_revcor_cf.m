function OutFile = samsrf_revcor_cf(Model, SrfFile, Roi)
%
% OutFile = samsrf_revcor_cf(Model, SrfFile, [Roi])
%
% Calculates the connective field profiles using reverse correlation:
% For each vertex in the region of interest Roi, it calculates the 
% correlation with all the time series in the seed region Model.SeedRoi.
% The most correlated vertices are the connective field of this vertex.
%
% Note that this procedure does not estimate CF size. In order to do this,
% you will need to run a post-processing function to fit CF models to the
% reverse correlation profiles.
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
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end

%% Default model parameters
Model = samsrf_model_defaults('samsrf_revcor_cf', Model);

%% Start time of analysis
t0 = tic; samsrf_newline;  
samsrf_disp('*** SamSrf connective field analysis ***');
[vn, vd] = samsrf_version;
samsrf_disp([' Version ' num2str(vn) ' - ' vd]);
samsrf_newline;
samsrf_disp('Current working directory:');
samsrf_disp([' ' pwd]);
samsrf_newline;
% Are we also fitting pRF model?
if Model.Fit_pRF == 1
    % Which optimisation algorithm is used?
    if isfield(Model, 'Hooke_Jeeves_Steps')
        % Hooke-Jeeves algorithm
        samsrf_disp('Using Hooke-Jeeves pattern search algorithm')
        hjs = [' with step sizes: '];
        for p = 1:length(Model.Hooke_Jeeves_Steps)
            hjs = [hjs num2str(Model.Hooke_Jeeves_Steps(p))];
            if p < length(Model.Hooke_Jeeves_Steps)
                hjs = [hjs ', '];
            end
        end
        samsrf_disp(hjs);
        if isfield(Model, 'Nelder_Mead_Tolerance')
            samsrf_disp('WARNING: (Nelder-Mead parameter tolerance was also defined but not used...)');
        end
    else
        % Nelder-Mead algorithm
        samsrf_disp('Using Nelder-Mead (fminsearch) algorithm');
        if isfield(Model, 'Nelder_Mead_Tolerance')
            samsrf_disp([' with parameter tolerance: ' num2str(Model.Nelder_Mead_Tolerance)]);
        else
            samsrf_disp(' with default parameter tolerance');
        end
    end
else
    if Model.Fit_pRF == 0
        % Convex hull algorithm instead of model fitting
        samsrf_disp('Using convex hull algorithm instead of model fitting')    
    elseif Model.Fit_pRF == -1
        % Summary statistics instead of model fitting
        samsrf_disp('Using summary statistics instead of model fitting')    
    else
        samsrf_error('Invalid value chosen for Model.Fit_pRF!');
    end
end
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
Srf.SeedVx = samsrf_loadlabel(Model.SeedRoi); % Store in Srf for posterity
X = Tc(:, Srf.SeedVx); % Time courses in seed ROI

%% Calculate cortical distances within seed ROI
if Model.Smoothing > 0 % No point when not smoothing
    Ds = samsrf_geomatrix(Srf.Vertices, Srf.Faces, Srf.SeedVx); % Cortical distance matrix
    Weights = exp(-(Ds.^2)/(2*Model.Smoothing.^2)); % Smoothing weight matrix
    samsrf_disp(['Distance matrix computed in ' num2str(toc(t0)/60) ' minutes.']);
else 
    Weights = [];
end

%% Load template map
samsrf_disp(['Loading template map: ' Model.Template]);
Temp = load(EnsurePath(Model.Template));
Temp.Srf = samsrf_expand_srf(Temp.Srf);
samsrf_newline;

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = []; % Needed for later

%% Determine CF parameters
% CF parameters
if Model.Fit_pRF == 1
    Srf.Data = zeros(8,size(Srf.Vertices,1)); % Output data when fitting pRFs
elseif Model.Fit_pRF == 0
    Srf.Data = zeros(9,size(Srf.Vertices,1)); % Output data when using convex hull algorithm 
elseif Model.Fit_pRF == -1
    Srf.Data = zeros(6,size(Srf.Vertices,1)); % Output data when calculating summary stats only
else
    samsrf_error('Invalid Fit_pRF value provided!');
end
samsrf_disp('Reverse correlation & CF parameter estimation...');
if Model.Fit_pRF == 1
    samsrf_disp(' Fitting pRF parameters to template pRF coordinates.');
    % Which fitting algorithm?
    if isfield(Model, 'Hooke_Jeeves_Steps')
        % Use Hooke-Jeeves algorithm
        AlgorithmParam = Model.Hooke_Jeeves_Steps; % Beta step size is pre-set 
    else
        % Use Nelder-Mead algorithm
        if isfield(Model, 'Nelder_Mead_Tolerance')
            AlgorithmParam = [NaN Model.Nelder_Mead_Tolerance]; % Define parameter tolerance
        else
            AlgorithmParam = NaN;
        end
    end    
else
    if Model.Fit_pRF == 0
        samsrf_disp(' Using convex hull algorithm to estimate parameters.');
        AlgorithmParam = [];
    elseif Model.Fit_pRF == -1
        samsrf_disp(' Using summary statistics to estimate parameters.');
        AlgorithmParam = Inf;
    end
end
% Run estimation loop
[fVimg, fXimg, fYimg, fWimg, fRimg, fSimg, fBimg, Rmaps] = samsrf_cfparam_loop(Tc(:,mver), X, Srf.Area, Srf.SeedVx, Temp.Srf.Data, Weights, AlgorithmParam, Model.Save_Rmaps);
t2 = toc(t0); 
samsrf_disp(['Parameter estimates completed in ' num2str(t2/60/60) ' hours.']);
samsrf_newline;

% Save as surface structure
Srf.Functional = 'Connective field';
if Model.Fit_pRF == 1
    % Fitting parameters
    Data = [fRimg; fXimg; fYimg; fSimg; fBimg; fWimg; fVimg];
    Srf.Values = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'Beta'; 'Baseline'; 'Fwhm'; 'Vx'};
elseif Model.Fit_pRF == 0
    % Convex hull estimation of parameters
    Data = [fRimg; fXimg; fYimg; fSimg; fBimg; fWimg; fVimg];
    Srf.Values = {'R^2'; 'x0'; 'y0'; 'Centre'; 'Surround'; 'Beta'; 'Suppression'; 'Fwhm'; 'Vx'};
elseif Model.Fit_pRF == -1
    % Using summary stats only
    Data = [fRimg; fXimg; fYimg; fSimg; fWimg; fVimg];
    Srf.Values = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'Fwhm'; 'Vx'};
end
Srf.Data(:,mver) = Data;

% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
end

% Are we saving correlation profiles?
if Model.Save_Rmaps
    samsrf_disp('Saving CF profiles in data file.');
    Srf.ConFlds = NaN(length(Srf.SeedVx), size(Srf.Vertices,1));
    Srf.ConFlds(:,mver) = Rmaps;
else
    samsrf_disp('Not saving CF profiles...');
    Srf.ConFlds = Rmaps; % No Rmaps to save space
end

% Save map files
samsrf_disp('Saving CF results...');
Srf = samsrf_compress_srf(Srf, mver);
save(OutFile, 'Model', 'Srf', '-v7.3');
samsrf_disp(['Saved ' OutFile '.mat']); 

% End time
t3 = toc(t0); 
EndTime = num2str(t3/60/60);
samsrf_newline; samsrf_disp(['Whole analysis completed in ' EndTime ' hours.']);
samsrf_disp('******************************************************************');
samsrf_newline; samsrf_newline;

