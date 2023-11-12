function OutFile = samsrf_fit_prf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_fit_prf(Model, SrfFiles, [Roi = ''])
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
%   Model:        Contains the parameters defining the eccentricity/scaling factor of the space.
%   SrfFiles:     Cell array of file names without extension (e.g. {'lh_Bars1' 'lh_Bars2'})
%                     Files will be concatenated in this order.
%   Roi:          ROI label to restrict the analysis (default = '') 
%                     Optional, but without this the analysis can take forever.
%
% Returns the name of the map file it saved.
%
% 06/07/2022 - New faster version using vector apertures instead of movies (DSS)
% 07/07/2022 - Added option for 10 free parameters (excluding CSS exponent) (DSS)
%              Fixed small bug with multiple maximal correlations in coarse fit (DSS)
% 10/08/2022 - Fixed bug with coarse fit when time series is flat (DSS)
% 24/01/2023 - Added more info to error when apertures aren't vectorised (DSS) 
% 17/03/2023 - Added option for conventional Dumoulin & Wandell approach to predict neural response (DSS)
% 24/10/2023 - Bugfix when using 32bit data for seeding fine-fit (DSS)
%              Bugfix for fitting betas when fits are bad (DSS)
% 13/11/2023 - Added option to estimate HRF parameters during fine-fitting (DSS)
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
t0 = tic; new_line;  
disp('*** SamSrf pRF model fitting ***');
[vn, vd] = samsrf_version;
disp([' Version ' num2str(vn) ' - ' vd]);
disp([' pRF model: ' PrfFcnName]);
new_line;
disp('Current working directory:');
disp([' ' pwd]);
new_line;
% Are we fine-fitting?
if ~Model.Coarse_Fit_Only
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
end
new_line;

%% Load apertures
disp('Load stimulus apertures...');
load(EnsurePath(Model.Aperture_File));  % Loads a variable called ApFrm
disp([' Loading ' Model.Aperture_File ': ' num2str(size(ApFrm,2)) ' volumes']);
if sum(ApFrm(:)<0) > 0
    disp(' Warning: Apertures contain negative values!');
end
if ~exist('ApXY', 'var') 
    error('Aperture pixel coordinates undefined! Did you use VectoriseApertures?');
end
if Model.Aperture_Mean_Response
    disp(' Using conventional Dumoulin & Wandell 2008 biophysical model to predict neural responses.');
    ApXY = [ApXY; NaN NaN]; % Conventional model is marked by NaN in final pixel of pRF profile
else
    disp(' Using default biophysical model to predict neural responses from percent pRF overlap.');    
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

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = [];  % Clear data field
% Do aperture & data length match?
if size(Tc,1)*Model.Downsample_Predictions ~= size(ApFrm,2)
    error('Mismatch between length of apertures and data!');
end

%% Load or generate HRF
disp('Haemodynamic response function...')
if isempty(Model.Hrf)
    disp(' Using canonical HRF');
    Model.Hrf = samsrf_hrf(Model.TR);
elseif FittingHrf
    disp(' Estimating HRF during fine-fit');
elseif isscalar(Model.Hrf) && Model.Hrf == 1
    disp(' No HRF used!');
else
    if ischar(Model.Hrf)
        disp([' Subject-specific HRF: ' Model.Hrf]);
        load([pwd filesep Model.Hrf]);
         % HRF based on loaded parameters but TR defined here
        Model.Hrf = samsrf_doublegamma(Model.TR, [fP(1:2) 1 1 fP(3) 0 32])' * fP(4);
    else
        disp(' Using Subject-specific HRF provided');
    end
end
new_line; 

%% Generate prediction matrix
if isempty(Model.Seed_Fine_Fit) % Only if running coarse fit
    if ~exist([pwd filesep SearchspaceFile], 'file') 
        disp('Generating predictions...');
        [X,S] = prf_generate_searchspace(Model.Prf_Function, ApFrm, ApXY, Model.Param1, Model.Param2, Model.Param3, Model.Param4, Model.Param5, ...
                                                                          Model.Param6, Model.Param7, Model.Param8, Model.Param9, Model.Param10, Model.Polar_Search_Space);    
        save(SearchspaceFile, 'X', 'S', '-v7.3');
        t1 = toc(t0); 
        disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
    else
        disp('Loading pre-defined predictions...');
        load([pwd filesep SearchspaceFile]);
        disp([' Loading ' SearchspaceFile]);
        % Does number of grid points match model?
        if size(S,2) ~= length(Model.Param1) * length(Model.Param2) * length(Model.Param3) * length(Model.Param4) * length(Model.Param5) ...
                      * length(Model.Param6) * length(Model.Param7) * length(Model.Param8) * length(Model.Param9) * length(Model.Param10)
            error('Mismatch between saved search space and model definition!');
        end
        % Does length of search space match apertures?
        if size(ApFrm,2) ~= size(X,1)
            error('Mismatch between length of saved search space and apertures!');
        end
        % Does length of search parameter & prediction matrix match?
        if size(S,2) ~= size(X,2)
            % Shouldn't ever happen unless someone screwed with the search space file
            error('Search space is corrupt! Mismatch between number of parameters & predictions!');
        end
    end
    disp(['Using search space with ' num2str(size(S,2)) ' grid points.']);
    new_line; 

    %% Convolution with HRF
    disp('Convolving predictions with HRF...');
    cX = NaN(size(Tc,1),size(X,2)); % Convolved X (has lower number of volumes than X, if downsampling) 
    for p = 1:size(X,2)        
        cX(:,p) = prf_convolve_hrf(X(:,p), Model.Hrf, Model.Downsample_Predictions); % Convolve each prediction & downsample if desired
    end
    X = cX; % Replace X with convolution
    new_line; 
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
  disp(['Loading ' Model.Seed_Fine_Fit ' to seed fine fit...']);
  SeedMap = load(EnsurePath(Model.Seed_Fine_Fit));
  SeedMap.Srf = samsrf_expand_srf(SeedMap.Srf);
  Pimg = SeedMap.Srf.Data(2:length(Model.Param_Names)+1,:); % Fitted parameter maps
  Rimg = SeedMap.Srf.Data(1,:); % R^2 map
  % In case saved as singles
  Pimg = double(Pimg);
  Rimg = double(Rimg);
else
  % Coarse fitting procedure
  disp('Coarse fitting...');
  Srf.X = zeros(size(Tc,1), size(Tc,2));  % Matrix with predictions
  Pimg = zeros(length(Model.Param_Names), size(Srf.Vertices,1)); % Fitted parameter maps
  Rimg = zeros(1, size(Srf.Vertices,1)); % R^2 map
  % If only running coarse fit
  if Model.Coarse_Fit_Only
      Bimg = zeros(2,size(Srf.Vertices,1)); % Beta map
  end
  
  % Which percentile of correlations to include in parameter estimate
  if Model.Coarse_Fit_Percentile == 100
    disp(' Using only the maximal coarse fit R^2 as parameter estimate');
  else
    disp([' Including top ' num2str(100-Model.Coarse_Fit_Percentile) '% of coarse fit R^2s in parameter estimates']);
  end
  
  % Loop through mask vertices (in blocks if Matlab R2012a or higher)
  disp([' Block size: ' num2str(cfvb) ' vertices']);
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
  disp(['Coarse fitting completed in ' num2str(t2/60) ' minutes.']);
end
new_line;

%% Run fine-fit or break now?
if Model.Coarse_Fit_Only 
    %% Store coarse fit parameters
    disp('Only running coarse fit!');
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
    disp('Fine fitting...');
    % Restrict matrices to mask
    Tc = Srf.Y(:,mver); % Always use unsmoothed data for fine fit
    Rimg = Rimg(1,mver);
    Pimg = Pimg(:,mver);
    
    % Are we modelling compressive nonlinearity?
    if Model.Compressive_Nonlinearity
        disp(' Modelling neural response with compressive spatial summation');
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
    disp(['Fine fitting completed in ' num2str(t3/60/60) ' hours.']);
    new_line;
    
    disp('Fitting betas & storing parameter estimates...');   
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

disp('Tidying up final results structure...');
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
new_line;

% Compress to save space
Srf = samsrf_compress_srf(Srf, mver);
% Save fine fit map files
disp('Saving pRF fitting results...');
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60/60);
new_line; disp(['Whole analysis completed in ' EndTime ' hours.']);
disp('******************************************************************');
new_line; new_line;
