function OutFile = samsrf_revcor_cf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_revcor_cf(Model, SrfFiles, [Roi])
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
%   Model:          Contains the parameters defining the analysis.
%   SrfFiles:       Cell array of file names without extension (e.g. {'lh_Bars1' 'lh_Bars2'})
%                       Files will be concatenated in this order.
%   Roi:            ROI label to restrict the analysis (default = '') 
%                       Optional, but without this the analysis can take forever.
%
% Returns the name of the map file it saved.
%
% 07/04/2022 - pRF fitting now thresholds correlations by half-maximum (DSS)
% 08/04/2022 - Improved algorithm to home in on pRF size estimates (DSS)
% 15/04/2022 - Warns if both Hooke-Jeeves steps & Nelder-Mead tolerance are defined (DSS)
%              Outsourced check for default parameters so no longer needs to check these (DSS)
%              Fixed error with duration report for generating distance matrix (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 06/05/2022 - Now can estimate CF parameters using convex hull algorithm (DSS)
%              Fixed inconsequential error with command window reports (DSS)
% 16/05/2022 - Streamlined correlation & parameter estimation for efficiency (DSS)
%              Added option to use summary statistics for parameter estimation (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end

%% Default model parameters
Model = samsrf_model_defaults('samsrf_revcor_cf', Model);

%% Start time of analysis
t0 = tic; new_line;  
disp('*** SamSrf connective field analysis ***');
[vn, vd] = samsrf_version;
disp([' Version ' num2str(vn) ' - ' vd]);
new_line;
disp('Current working directory:');
disp([' ' pwd]);
new_line;
% Are we also fitting pRF model?
if Model.Fit_pRF == 1
    % Which optimisation algorithm is used?
    if isfield(Model, 'Hooke_Jeeves_Steps')
        % Hooke-Jeeves algorithm
        disp('Using Hooke-Jeeves pattern search algorithm')
        hjs = [' with step sizes: '];
        for p = 1:length(Model.Hooke_Jeeves_Steps)
            hjs = [hjs num2str(Model.Hooke_Jeeves_Steps(p))];
            if p < length(Model.Hooke_Jeeves_Steps)
                hjs = [hjs ', '];
            end
        end
        disp(hjs);
        if isfield(Model, 'Nelder_Mead_Tolerance')
            warning('(Nelder-Mead parameter tolerance was also defined but isn''t used...)');
        end
    else
        % Nelder-Mead algorithm
        disp('Using Nelder-Mead (fminsearch) algorithm');
        if isfield(Model, 'Nelder_Mead_Tolerance')
            disp([' with parameter tolerance: ' num2str(Model.Nelder_Mead_Tolerance)]);
        else
            disp(' with default parameter tolerance');
        end
    end
else
    if Model.Fit_pRF == 0
        % Convex hull algorithm instead of model fitting
        disp('Using convex hull algorithm instead of model fitting')    
    elseif Model.Fit_pRF == -1
        % Summary statistics instead of model fitting
        disp('Using summary statistics instead of model fitting')    
    else
        error('Invalid value chosen for Model.Fit_pRF!');
    end
end
new_line;

%% Load images 
if ischar(SrfFiles)
    SrfFiles = {SrfFiles};
end
disp('Reading surface images...')
Tc = [];
for f = 1:length(SrfFiles)
    % Load surface image from each run
    load([pwd filesep SrfFiles{f}]);
    Srf = samsrf_expand_srf(Srf);
    if f == 1
        OutFile = [Srf.Hemisphere '_' Model.Name];
    end
    disp([' Loading ' SrfFiles{f} ': ' num2str(size(Srf.Vertices,1)) ' vertices & ' num2str(size(Srf.Data,1)) ' volumes']);
    Tc = [Tc; Srf.Data]; % Add run to time course
end

%% Correct by global mean signal?
if Model.Global_Signal_Correction
    new_line; disp('Applying global signal correction...');
    Srf.Data = Tc;
    Srf = samsrf_removenoise(Srf, nanmean(Srf.Data,2));
    Tc = Srf.Data;
    Srf.Data = [];
end

%% Load ROI mask
if isempty(Roi)
    new_line; disp('Using all vertices (''tis gonna take forever!)...');
    mver = 1:size(Srf.Vertices,1);
else
    new_line; disp('Reading ROI mask...')
    mver = samsrf_loadlabel(Roi);
    disp([' Loading ' Roi ': ' num2str(size(mver,1)) ' vertices']);
end
new_line; 

%% Limit data due to noise ceiling?
if isfield(Srf, 'Noise_Ceiling')
    if Model.Noise_Ceiling_Threshold > 0
        mver = mver(Srf.Noise_Ceiling(mver) > Model.Noise_Ceiling_Threshold);
        disp(['Limiting analysis to ' num2str(length(mver)) ' vertices above noise ceiling ' num2str(Model.Noise_Ceiling_Threshold)]);
        new_line;
    end
end

%% Load seed ROI 
disp(['Loading seed ROI: ' Model.SeedRoi]);
Srf.SeedVx = samsrf_loadlabel(Model.SeedRoi); % Store in Srf for posterity
X = Tc(:, Srf.SeedVx); % Time courses in seed ROI

%% Calculate cortical distances within seed ROI
if Model.Smoothing > 0 % No point when not smoothing
    Ds = samsrf_geomatrix(Srf.Vertices, Srf.Faces, Srf.SeedVx); % Cortical distance matrix
    Weights = exp(-(Ds.^2)/(2*Model.Smoothing.^2)); % Smoothing weight matrix
    disp(['Distance matrix computed in ' num2str(toc(t0)/60) ' minutes.']);
else 
    Weights = [];
end

%% Load template map
disp(['Loading template map: ' Model.Template]);
Temp = load(EnsurePath(Model.Template));
Temp.Srf = samsrf_expand_srf(Temp.Srf);
new_line;

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = []; % Needed for later

%% Determine CF parameters
% CF parameters
if Model.Fit_pRF > 0
    Srf.Data = zeros(8,size(Srf.Vertices,1)); % Output data when fitting pRFs
else
    Srf.Data = zeros(6,size(Srf.Vertices,1)); % Output data when not fitting pRFs
end
disp('Reverse correlation & CF parameter estimation...');
if Model.Fit_pRF == 1
    disp(' Fitting pRF parameters to template pRF coordinates.');
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
        disp(' Using convex hull algorithm to estimate parameters.');
        AlgorithmParam = [];
    elseif Model.Fit_pRF == -1
        disp(' Using summary statistics to estimate parameters.');
        AlgorithmParam = Inf;
    end
end
% Run estimation loop
[fVimg, fXimg, fYimg, fWimg, fRimg, fSimg, fBimg, Rmaps] = samsrf_cfparam_loop(Tc(:,mver), X, Srf.Area, Srf.SeedVx, Temp.Srf.Data, Weights, AlgorithmParam, Model.Save_Rmaps);
t2 = toc(t0); 
disp(['Parameter estimates completed in ' num2str(t2/60/60) ' hours.']);
new_line;

% Save as surface structure
Srf.Functional = 'Connective field';
if Model.Fit_pRF == 1
    % Fitting pRF parameters
    Data = [fRimg; fXimg; fYimg; fSimg; fBimg; fWimg; fVimg];
    Srf.Values = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'Beta'; 'Baseline'; 'Fwhm'; 'Vx'};
else
    % Not fitting pRF parameters
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
    disp('Saving CF profiles in data file.');
    Srf.ConFlds = NaN(length(Srf.SeedVx), size(Srf.Vertices,1));
    Srf.ConFlds(:,mver) = Rmaps;
else
    disp('Not saving CF profiles...');
    Srf.ConFlds = Rmaps; % No Rmaps to save space
end

% Save map files
disp('Saving CF results...');
Srf = samsrf_compress_srf(Srf, mver);
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t3 = toc(t0); 
EndTime = num2str(t3/60/60);
new_line; disp(['Whole analysis completed in ' EndTime ' hours.']);
disp('******************************************************************');
new_line; new_line;

