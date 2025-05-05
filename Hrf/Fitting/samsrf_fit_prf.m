function OutFile = samsrf_fit_prf(Model, SrfFile, Roi)
%
% OutFile = samsrf_fit_prf(Model, SrfFile, [Roi = ''])
%
% Fits a pRF model using the slow fine-fitting procedure (fminsearch in Optimization toolbox).
%
% IMPORTANT NOTE: In SamSrf 7 the precision of the model fitting algorithm was improved!
%                 Model fits using this version are by default -NOT- compatible with previous versions!
%                 Under some conditions previous versions sometimes reverted to the search grid.
%                   
%                 In SamSrf 8, the option to use the Hooke-Jeeves algorithm was added.
%                 Moreover, you can now define parameter tolerance for Nelder-Mead algorithm,
%                 but note that both of these can result in imprecise parameter esstimates.
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
%                     Optional, but without this the analysis can take forever depending on how much data there is.
%
% Returns the name of the map file it saved.
%
% 04/09/2024 - Now automatically vectorises apertures if necessary (DSS)  
%              Instead of a list of files, you can now specify a Srf as input (DSS)
% 11/09/2024 - Fixed bug with Hooke-Jeeves algorithm step sizes not being scaled (DSS)
%              Changed selection of separate HRF file (DSS)
% 12/09/2024 - Added Hrf=0 option for using SPM canonical HRF (DSS)
% 13/09/2024 - Removed option to provife data file list for concatenation (DSS)
%              Fixed help descriptions (DSS)
% 15/09/2024 - Fixed bug with undefined output filename when providing Srf data (DSS)
% 18/09/2024 - Fixed bug with pRF-from-CF analysis using incorrect search space (DSS)
% 20/10/2024 - Fixed bug when using SPM canonical HRF (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end
% Define search space name
SearchspaceFile = ['src_' Model.Name '.mat'];  
% pRF function name
PrfFcnName = char(Model.Prf_Function);
bc = strfind(PrfFcnName, ')');
bo = strfind(PrfFcnName, '(');
PrfFcnName = PrfFcnName(bc(1)+1:bo(2)-1); % Remove rubbish

%% Default model parameters
Model = samsrf_model_defaults('samsrf_fit_prf', Model);
% Estimating HRF in fine-fit?
if isinf(Model.Hrf)
    % Flexible HRF fitting
    FittingHrf = true;
    Model.Hrf = samsrf_doublegamma(Model.TR); % Use SPM's canonical for coarse-fit 
else
    % Predefined HRF
    FittingHrf = false;
end

%% MatLab R2012a or higher can do fast coarse-fit
if verLessThan('matlab','7.13')
    cfvb = 1;
else
    % Coarse-fit vertex block size
    cfvb = Model.Coarse_Fit_Block_Size;
end

%% Start time of analysis
t0 = tic; samsrf_newline;  
samsrf_disp('*** SamSrf pRF model fitting ***');
[vn, vd] = samsrf_version;
samsrf_disp([' Version ' num2str(vn) ' - ' vd]);
samsrf_disp([' pRF model: ' PrfFcnName]);
samsrf_newline;
samsrf_disp('Current working directory:');
samsrf_disp([' ' pwd]);
samsrf_newline;
% Are we fine-fitting?
if ~Model.Coarse_Fit_Only
    % Which optimisation algorithm is used?
    if isfield(Model, 'Hooke_Jeeves_Steps') 
        % Hooke-Jeeves algorithm
        samsrf_disp('Using Hooke-Jeeves pattern search algorithm')
        hjs = [' with step sizes: '];
        for p = 1:length(Model.Hooke_Jeeves_Steps)
            % Scale step size if necessary
            if Model.Scaled_Param(p)
                Model.Hooke_Jeeves_Steps(p) = Model.Hooke_Jeeves_Steps(p) * Model.Scaling_Factor;
            end
            % Output step sizes
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
end
samsrf_newline;

%% Prepare apertures
if ~isfield(Model, 'SeedRoi')
    % If loading apertures for conventional pRF analysis
    samsrf_disp('Load stimulus apertures...');
    load(EnsurePath(Model.Aperture_File));  % Supposed to loads variables called ApFrm & ApXY
    samsrf_disp([' Loading ' Model.Aperture_File]); 
    if sum(ApFrm(:)<0) > 0
        samsrf_disp(' Warning: Apertures contain negative values!');
    end

    % Apertures vectorised already?
    if ~exist('ApXY', 'var') 
        samsrf_disp('Aperture pixel coordinates undefined. Vectorising apertures...');
        VectoriseApertures(EnsurePath(Model.Aperture_File));
        Model.Aperture_File = [Model.Aperture_File '_vec']; % Use vectorised apertures instead
        load(EnsurePath(Model.Aperture_File));  % Loads variables called ApFrm & ApXY
    end

    % Rescale apertures
    ApXY = ApXY / max(abs(ApXY(:))); % Normalise scale to maximum
    ApXY = ApXY * Model.Scaling_Factor; % Rescale to current scaling factor

    % How many volumes?
    samsrf_disp([' There are ' num2str(size(ApFrm,2)) ' volumes']);
else
    % If backprojecting seed activity for pRF-from-CF analysis
    samsrf_disp('Backprojecting CF apertures...');
    if isfield(Model, 'Template') && ~isempty(Model.Template)
        [ApName, ApMax] = BackprojAps(Model.Template, Model.SeedRoi, SrfFile);
    else
        samsrf_error('No template map defined for pRF-from-CF analysis!');
    end
    % Set defaults for pRF-from-CF analysis
    Model.Aperture_File = ApName;
    Model.Scaling_Factor = ApMax;
    clear ApName ApMax
    % Rescale search space as needed
    if Model.Polar_Search_Space
        p = 2;
    else
        p = 1;
    end
    % Loop thru parameters
    for i = p:length(Model.Param_Names)
        if Model.Scaled_Param(i)
            Model.(['Param' num2str(i)]) = Model.(['Param' num2str(i)]) * Model.Scaling_Factor;
        end
    end    
    % Load backprojected apertures
    load(EnsurePath(Model.Aperture_File));  % Supposed to loads variables called ApFrm & ApXY
end

% Which biophysical model to use
if Model.Aperture_Mean_Response
    samsrf_disp(' Using Dumoulin & Wandell 2008 biophysical model (% aperture response)');
    ApXY = [ApXY; NaN NaN]; % Conventional model is marked by NaN in final pixel of pRF profile
else
    samsrf_disp(' Using SamSrf biophysical model (% pRF overlap)');    
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
samsrf_newline; 

%% Limit data due to noise ceiling?
if isfield(Srf, 'Noise_Ceiling')
    if Model.Noise_Ceiling_Threshold > 0
        mver = mver(Srf.Noise_Ceiling(mver) > Model.Noise_Ceiling_Threshold);
        samsrf_disp(['Limiting analysis to ' num2str(length(mver)) ' vertices above noise ceiling ' num2str(Model.Noise_Ceiling_Threshold)]);
        samsrf_newline;
    end
end

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = [];  % Clear data field
% Do aperture & data length match?
if size(Tc,1)*Model.Downsample_Predictions ~= size(ApFrm,2)
    samsrf_error('Mismatch between length of apertures and data!');
end

%% Load or generate HRF
samsrf_disp('Haemodynamic response function...')
if isempty(Model.Hrf) 
    samsrf_disp(' Using de Haas canonical HRF');
    Model.Hrf = samsrf_hrf(Model.TR);
elseif isscalar(Model.Hrf) && Model.Hrf == 0
    samsrf_disp(' Using SPM canonical HRF');
    Model.Hrf = samsrf_doublegamma(Model.TR, [6 16 1 1 6 0 32]);
elseif FittingHrf
    samsrf_disp(' Estimating HRF during fine-fit');
elseif isscalar(Model.Hrf) && Model.Hrf == 1
    samsrf_disp(' No HRF used!');
else
    if ischar(Model.Hrf)
        samsrf_disp([' Subject-specific HRF: ' Model.Hrf]);
        load(Model.Hrf);
         % HRF based on loaded parameters but TR defined here
        Model.Hrf = samsrf_doublegamma(Model.TR, [fP(1:2) 1 1 fP(3) 0 32])' * fP(4);
    else
        samsrf_disp(' Using Subject-specific HRF provided');
    end
end
samsrf_newline; 

%% Generate prediction matrix
if isempty(Model.Seed_Fine_Fit) % Only if running coarse fit
    if ~exist([pwd filesep SearchspaceFile], 'file') 
        samsrf_disp('Generating predictions...');
        [X,S] = prf_generate_searchspace(Model.Prf_Function, ApFrm, ApXY, Model.Param1, Model.Param2, Model.Param3, Model.Param4, Model.Param5, ...
                                                                          Model.Param6, Model.Param7, Model.Param8, Model.Param9, Model.Param10, Model.Polar_Search_Space);    
        save(SearchspaceFile, 'X', 'S', '-v7.3');
        t1 = toc(t0); 
        samsrf_disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
    else
        samsrf_disp('Loading pre-defined predictions...');
        load([pwd filesep SearchspaceFile]);
        samsrf_disp([' Loading ' SearchspaceFile]);
        % Does number of grid points match model?
        if size(S,2) ~= length(Model.Param1) * length(Model.Param2) * length(Model.Param3) * length(Model.Param4) * length(Model.Param5) ...
                      * length(Model.Param6) * length(Model.Param7) * length(Model.Param8) * length(Model.Param9) * length(Model.Param10)
            samsrf_error('Mismatch between saved search space and model definition!');
        end
        % Does length of search space match apertures?
        if size(ApFrm,2) ~= size(X,1)
            samsrf_error('Mismatch between length of saved search space and apertures!');
        end
        % Does length of search parameter & prediction matrix match?
        if size(S,2) ~= size(X,2)
            % Shouldn't ever happen unless someone screwed with the search space file
            samsrf_error('Search space is corrupt! Mismatch between number of parameters & predictions!');
        end
    end
    samsrf_disp(['Using search space with ' num2str(size(S,2)) ' grid points.']);
    samsrf_newline; 

    %% Convolution with HRF
    samsrf_disp('Convolving predictions with HRF...');
    cX = NaN(size(Tc,1),size(X,2)); % Convolved X (has lower number of volumes than X, if downsampling) 
    for p = 1:size(X,2)        
        cX(:,p) = prf_convolve_hrf(X(:,p), Model.Hrf, Model.Downsample_Predictions); % Convolve each prediction & downsample if desired
    end
    X = cX; % Replace X with convolution
    samsrf_newline; 
end

%% Smooth for coarse fit?
if Model.Smoothed_Coarse_Fit > 0
    Srf.Data = Tc; % Store full time course in Srf
    Srf = samsrf_smooth_sphere(Srf, Model.Smoothed_Coarse_Fit, Roi, 0); % Smooth within ROI
    Srf = rmfield(Srf, 'Raw_Data'); % Remove raw data field
    Tc = Srf.Data; % Put smoothed data back into time course variable
    Srf.Data = []; % Clear data field
end

%% Coarse fit / Load seed map
if ~isempty(Model.Seed_Fine_Fit)
  % Load a previous map as seeds for fine fit
  samsrf_disp(['Loading ' Model.Seed_Fine_Fit ' to seed fine fit...']);
  SeedMap = load(EnsurePath(Model.Seed_Fine_Fit));
  SeedMap.Srf = samsrf_expand_srf(SeedMap.Srf);
  Pimg = SeedMap.Srf.Data(2:length(Model.Param_Names)+1,:); % Fitted parameter maps
  Rimg = SeedMap.Srf.Data(1,:); % R^2 map
  % In case saved as singles
  Pimg = double(Pimg);
  Rimg = double(Rimg);
else
  % Coarse fitting procedure
  samsrf_disp('Coarse fitting...');
  Srf.X = zeros(size(Tc,1), size(Tc,2));  % Matrix with predictions
  Pimg = zeros(length(Model.Param_Names), size(Srf.Vertices,1)); % Fitted parameter maps
  Rimg = zeros(1, size(Srf.Vertices,1)); % R^2 map
  % If only running coarse fit
  if Model.Coarse_Fit_Only
      Bimg = zeros(2,size(Srf.Vertices,1)); % Beta map
  end
  
  % Which percentile of correlations to include in parameter estimate
  if Model.Coarse_Fit_Percentile == 100
    samsrf_disp(' Using only the maximal coarse fit R^2 as parameter estimate');
  else
    samsrf_disp([' Including top ' num2str(100-Model.Coarse_Fit_Percentile) '% of coarse fit R^2s in parameter estimates']);
  end
  
  % Loop through mask vertices (in blocks if Matlab R2012a or higher)
  samsrf_disp([' Block size: ' num2str(cfvb) ' vertices']);
  samsrf_progbar(0);
  for vs = 1:cfvb:length(mver)
      % Starting index of current vertex block
      ve = vs + cfvb - 1;
      if ve > length(mver)
          ve = length(mver);
      end
      vx = mver(vs:ve);
      % Find best prediction
      Y = Tc(:,vx);  % Time course of current vertex
      if Model.Only_Positive_Coarse_Fits
      	 R = corr(Y,X); % Best correlating prediction 
      	 mR = max(R,[],2).^2; % Find best fit & square now
      	 R = R.^2; % Now turn others into R^2 too
      else
      	 R = corr(Y,X).^2; % Best correlating prediction (squared to allow for negative betas!)
      	 mR = max(R,[],2); % Find best fit
      end
      for v = 1:length(vx)
          rx = R(v,:) >= prctile(R(v,:), Model.Coarse_Fit_Percentile); % Predictions with correlation above percentile
          % Store parameters
          Pimg(:,vx(v)) = mean(S(1:length(Model.Param_Names),rx),2); % Mean parameter across top percentile predictions

          % Ensure vector isn't all false
          if sum(rx) > 0
              % Store prediction
              if Model.Coarse_Fit_Percentile == 100
                % Only take maximum
                Srf.X(:,vx(v)) = X(:,find(rx,1));  % Best fitting convolved prediction
                Rimg(1,vx(v)) = mR(v);  % Variance explained at maximum
              else
                % Top percentile of predictions
                Rfp = Model.Prf_Function(Pimg(:,vx(v))', ApXY); % New pRF profile replicated by number of volumes
                nX = prf_predict_timecourse(Rfp, ApFrm); % New predicted time courase
                nX = prf_convolve_hrf(nX, Model.Hrf, Model.Downsample_Predictions); % Convolve new time course with HRF
                Rimg(1,vx(v)) = corr(Y(:,v),nX)^2; % New goodness of fit
                Srf.X(:,vx(v)) = nX;  % Store new prediction with convolution
              end              

              % If running coarse fit only determine betas now
              if Model.Coarse_Fit_Only 
                % Fit betas for amplitude & intercept
                warning off % In case of rank deficient GLM
                B = [ones(length(Y(:,v)),1) Srf.X(:,vx(v))] \ Y(:,v); % GLM fit 
                warning on
                Bimg(1,vx(v)) = B(2); % Amplitude
                Bimg(2,vx(v)) = B(1); % Intercept
              end           
          end
          samsrf_progbar((vs+v-1)/length(mver));
      end
  end
  t2 = toc(t0); 
  samsrf_disp(['Coarse fitting completed in ' num2str(t2/60) ' minutes.']);
end
samsrf_newline;

%% Run fine-fit or break now?
if Model.Coarse_Fit_Only 
    %% Store coarse fit parameters
    samsrf_disp('Only running coarse fit!');
    fPimg = Pimg(:,mver); % Coarse fitted parameter maps
    fRimg = Rimg(1,mver); % R^2 map
    fBimg = Bimg(:,mver); % Beta map
    OutFile = [OutFile '_CrsFit']; % Suffix to indicate coarse fit only 
    % If fitting HRF restore flag in Model
    if FittingHrf
        Model.Hrf = Inf;
    end
else
    %% Run fine fit for each vertex
    samsrf_disp('Fine fitting...');
    % Restrict matrices to mask
    Tc = Srf.Y(:,mver); % Always use unsmoothed data for fine fit
    Rimg = Rimg(1,mver);
    Pimg = Pimg(:,mver);
    
    % Are we modelling compressive nonlinearity?
    if Model.Compressive_Nonlinearity
        samsrf_disp(' Modelling neural response with compressive spatial summation');
        Model.Param_Names{end+1} = 'Exponent'; % Add name for CSS exponent
        Model.Scaled_Param(end+1) = 0; % Exponent isn't scaled
        Model.Only_Positive(end+1) = 1; % Exponent cannot be negative
        Pimg = [Pimg; ones(1,size(Pimg,2))]; % Add seed parameters for exponent
        OutFile = [OutFile '_Css']; % Suffix to indicate compressive nonlinearity 
    end
    
    % Are we estimating HRF parameters?
    if FittingHrf
        HrfParams = {'RLat' 'ULat' 'RDisp' 'UDisp' 'R/U'};
        for p = 1:5
            Model.Param_Names{end+1} = HrfParams{p}; % Add name for CSS exponent
            Model.Scaled_Param(end+1) = 0; % HRF parameters aren't scaled
            Model.Only_Positive(end+1) = 1; % HRF parameters cannot be negative
        end
        Pimg = [Pimg; repmat([6 16 1 1 6]', 1, size(Pimg,2))]; % Seed parameters for HRF fit
        Model.Hrf = -Model.TR; % Negative TR indicates HRF fit
    end
    
    % Run fine fit (must be separate function for parallel computing toolbox
    [fPimg, fRimg] = samsrf_prfoptim_loop(Model, Tc, ApFrm, ApXY, Rimg, Pimg);
    if FittingHrf
        Model.Hrf = Inf;
    end
    t3 = toc(t0);
    samsrf_disp(['Fine fitting completed in ' num2str(t3/60/60) ' hours.']);
    samsrf_newline;
    
    samsrf_disp('Fitting betas & storing parameter estimates...');   
    % Additional data fields
    fBimg = zeros(2,length(mver)); % Beta maps
    Srf.X = zeros(size(Tc,1)*Model.Downsample_Predictions, size(Srf.Vertices,1)); % Matrix with unconvolved predictions
    
    % Process & fit betas for mask vertices
    samsrf_progbar(0);
    for v = 1:length(mver)
        % Only accept non-zero fits
        if fRimg(v) > 0
            IsGoodFit = true;
        else
            IsGoodFit = false;
        end
        
        if IsGoodFit
            % Loop thru fitted parameters
            for p = 1:size(fPimg,1)
                % If scaled parameter is out of bounds
                if Model.Scaled_Param(p) && abs(fPimg(p,v)) > 2 * Model.Scaling_Factor
                    IsGoodFit = false;
                end
                % If a parameter is negative but shouldn't be
                if Model.Only_Positive(p) && fPimg(p,v) <= 0
                    IsGoodFit = false;
                end
            end
        end

        % Store parameters & predicted time course 
        if IsGoodFit 
            % Only keep good parameters
            Rfp = Model.Prf_Function(fPimg(:,v)', ApXY); % Generate predicted pRF
            fX = prf_predict_timecourse(Rfp, ApFrm); % Generate predicted time course 
            if Model.Compressive_Nonlinearity
                fX = fX.^fPimg(end,v); % Compressive spatial summation
            end
            % Store prediction without HRF convolution
            Srf.X(:,mver(v)) = fX;  % Best fitting prediction
        else
            % Replace bad fits with coarse fit?
            if Model.Replace_Bad_Fits
                % Generate predicted time course
                Rfp = Model.Prf_Function(Pimg(:,v)', ApXY);
                fX = prf_predict_timecourse(Rfp, ApFrm);
                % Store prediction without HRF convolution
                Srf.X(:,mver(v)) = fX;  % Best fitting prediction
                % Replace bad fit with coarse fit
                fPimg(:,v) = Pimg(:,v); % Parameters
                fRimg(1,v) = Rimg(1,v); % Goodness-of-fit
            else
                % Just set all to zero
                fRimg(1,v) = 0;
                fPimg(:,v) = 0;
            end
        end

        % Convolve with HRF & fit betas
        if IsGoodFit || Model.Replace_Bad_Fits
            % Time course of current vertex
            Y = Tc(:,v);  
            % Convolve predicted time course with HRF & downsample if desired
            if FittingHrf
                fX = prf_convolve_hrf(fX, samsrf_doublegamma(Model.TR, fPimg(end-4:end,v)), Model.Downsample_Predictions);
            else
                fX = prf_convolve_hrf(fX, Model.Hrf, Model.Downsample_Predictions);
            end
            % Fit betas for amplitude & intercept
            warning off % In case of rank deficient GLM
            fB = [ones(length(Y),1) fX] \ Y; % GLM fit
            warning on
            fBimg(:,v) = fB([2 1]); % Amplitude & intercept
        end
        
        % Progress report
        samsrf_progbar(v/length(mver));
    end
end

samsrf_disp('Tidying up final results structure...');
% Prepare surface structure
Srf.Functional = PrfFcnName; % pRF function name
Data = [fRimg; fPimg; fBimg];
Srf.Data = zeros(size(Data,1), size(Srf.Vertices,1));
Srf.Data(:,mver) = Data; % Add parameter maps into full data matrix
% Add parameter names
Srf.Values = {'R^2'}; 
for p = 1:length(Model.Param_Names)
    Srf.Values{p+1} = Model.Param_Names{p}; 
end
Srf.Values{end+1} = 'Beta'; % Amplitude
Srf.Values{end+1} = 'Baseline'; % Intercept
% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
    Srf = samsrf_normr2(Srf); % Calculate normalised R^2
end
Srf.Values = Srf.Values'; % So it is the same as Srf.Data
samsrf_newline;

% Compress to save space
Srf = samsrf_compress_srf(Srf, mver);
% Save fine fit map files
samsrf_disp('Saving pRF fitting results...');
save(OutFile, 'Model', 'Srf', '-v7.3');
samsrf_disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60/60);
samsrf_newline; samsrf_disp(['Whole analysis completed in ' EndTime ' hours.']);
samsrf_disp('******************************************************************');
samsrf_newline; samsrf_newline;

