function OutFile = samsrf_fit_prf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_fit_prf(Model, SrfFiles, [Roi = ''])
%
% Fits a pRF model using the slow fine-fitting procedure (fminsearch in Optimization toolbox).
%
%   Model:          Contains the parameters defining the eccentricity/scaling factor of the space.
%   SrfFiles:       Cell array of file names without extension (e.g. {'lh_Bars1' 'lh_Bars2'})
%                       Files will be concatenated in this order.
%   Roi:            ROI label to restrict the analysis (default = '') 
%                       Optional, but without this the analysis can take forever.
%
% You must be in the folder containing the surface data files as well as
% the various parameter files (i.e. apertures, searchspace, and HRF).
%
% Returns the name of the map file it saved.
%
% 21/06/2018 - SamSrf 6 version (DSS)
% 09/11/2018 - Added option to use coarse fit parameters for bad fits (DSS)
% 28/11/2018 - Added support for noise ceiling, if calculated (DSS)
% 02/12/2018 - Added back an option for smoothed coarse-fit as in versions < 6 
%              Added option for only running the coarse-fit (DSS)
% 10/12/2018 - Fixed bug when coarse-fit only flag was undefined (DSS)
% 17/02/2020 - Added option for thresholding for what goes into fine fit (DSS)
% 18/02/2020 - Added option to use a seed map for fine fit (DSS)
% 20/02/2020 - Turned off searchspace generation when using seed map (DSS)
% 21/02/2020 - If only coarse fit is run the file name is suffixed with '_CrsFit' (DSS)
% 03/04/2020 - Removed all dependencies on spm_hrf (DSS)
%

%% Defaults & constants
% Check if waitbar is to be used
wb = samsrf_waitbarstatus;
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
if ~isfield(Model, 'Seed_Fine_Fit')
    Model.Seed_Fine_Fit = ''; % Use no seed map, so run coarse fit instead
end 
if ~isfield(Model, 'Replace_Bad_Fits')
    Model.Replace_Bad_Fits = false; % Don't replace bad fine fits with coarse fit
end
if ~isfield(Model, 'Smoothed_Coarse_Fit')
    Model.Smoothed_Coarse_Fit = 0; % No smoothing on coarse fit
end
if ~isfield(Model, 'Coarse_Fit_Only')
    Model.Coarse_Fit_Only = false; % Only run coarse fit & then save
end
if ~isfield(Model, 'Fine_Fit_Threshold')
    Model.Fine_Fit_Threshold = 0.01; % Include coarse fits with R^2>0.01 in fine fit
end

%% If coarse fit only suffix filename
if Model.Coarse_Fit_Only 
    if ~isempty(Model.Seed_Fine_Fit) % In case stupid choices were made
        error('No point running only coarse fit when seeding the fine fit!');
    end
end

%% MatLab R2012a or higher can do fast coarse-fit
if verLessThan('matlab','7.13')
    cfvb = 1;
else
    % Coarse-fit vertex block size
    cfvb = 10000;
end

%% Start time of analysis
t0 = tic; new_line;  
disp('*** SamSrf pRF model fitting ***');
disp([' pRF model: ' PrfFcnName]);
new_line;

%% Load apertures
disp('Load stimulus apertures...');
load(Model.Aperture_File);  % Loads a variable called ApFrm
disp([' Loading ' Model.Aperture_File '.']);
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
    disp([' Loading ' SrfFiles{f} ': ' num2str(size(Srf.Vertices,1)) ' vertices']);
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

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
% Store raw time courses
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = [];  % Clear data field

% Do aperture & data length match?
if size(Tc,1) ~= size(ApFrm,3)
    error('Mismatch between length of apertures and data!');
end

% Smooth for coarse fit?
if Model.Smoothed_Coarse_Fit > 0
    Srf.Data = Tc; % Store full time course in Srf
    Srf = samsrf_smooth_sphere(Srf, Model.Smoothed_Coarse_Fit, Roi, 0); % Smooth within ROI
    Srf = rmfield(Srf, 'Raw_Data'); % Remove raw data field
    Tc = Srf.Data; % Put smoothed data back into time course variable
    Srf.Data = []; % Clear data field
end

%% Load or generate HRF
disp('Haemodynamic response function...')
if isempty(Model.Hrf)
    disp(' Using canonical HRF');
    Model.Hrf = samsrf_hrf(Model.TR);
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
        [X,S] = prf_generate_searchspace(Model.Prf_Function, ApFrm, Model.Param1, Model.Param2, Model.Param3, Model.Param4, Model.Param5, wb);    
        save(SearchspaceFile, 'X', 'S', '-v7.3');
        t1 = toc(t0); 
        disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
    else
        disp('Loading pre-defined predictions...');
        load([pwd filesep SearchspaceFile]);
        disp([' Loading ' SearchspaceFile]);
        % Does length of search space match apertures?
        if size(ApFrm,3) ~= size(X,1)
            error('Mismatch between length of search space and apertures!');
        end
    end
    new_line; 

    %% Convolution with HRF
    disp('Convolving predictions with HRF...');
    for p = 1:size(X,2)
        cX = conv(X(:,p), Model.Hrf);
        X(:,p) = cX(1:size(Tc,1)); % Truncate back to original length
    end
    new_line; 
end

%% Coarse fit / Load seed map
if ~isempty(Model.Seed_Fine_Fit)
  % Load a previous map as seeds for fine fit
  disp(['Loading ' Model.Seed_Fine_Fit ' to seed fine fit...']);
  SeedMap = load(Model.Seed_Fine_Fit);
  SeedMap.Srf = samsrf_expand_srf(SeedMap.Srf);
  Pimg = SeedMap.Srf.Data(2:length(Model.Scaled_Param)+1,:); % Fitted parameter maps
  Rimg = SeedMap.Srf.Data(1,:); % R^2 map
  % Renormalise the data
  for i = 1:size(Pimg,1)
      if Model.Scaled_Param(i)
          Pimg(i,:) = Pimg(i,:) / Model.Scaling_Factor;  % Renormalise this parameter
      end
  end  
else
  % Coarse fitting procedure
  disp('Coarse fitting...');
  if wb h = waitbar(0, ['Coarse fitting... ' strrep(OutFile,'_','-')], 'Units', 'pixels', 'Position', [100 100 360 70]); end
  Srf.X = zeros(size(Tc));  % Matrix with predictions
  Pimg = zeros(length(Model.Param_Names), size(Srf.Vertices,1)); % Fitted parameter maps
  Rimg = zeros(1, size(Srf.Vertices,1)); % R^2 map
  % If only running coarse fit
  if Model.Coarse_Fit_Only
      Bimg = zeros(2,size(Srf.Vertices,1)); % Beta map
  end

  % Loop through mask vertices (in blocks if Matlab R2012a or higher)
  for vs = 1:cfvb:length(mver)
      % Starting index of current vertex block
      ve = vs + cfvb - 1;
      if ve > length(mver)
          ve = length(mver);
      end
      vx = mver(vs:ve);
      % Find best prediction
      Y = Tc(:,vx);  % Time course of current vertex
      R = corr(Y,X).^2; % Best correlating prediction (squared to allow for negative betas!)
      mR = max(R,[],2); % Find best fit
      for v = 1:length(vx)
          rx = find(R(v,:) == mR(v)); % Matrix position of best prediction
          if ~isempty(rx)
              rx = rx(1); % Only first instance 
              % Store prediction
              Srf.X(:,vx(v)) = X(:,rx);  % Best fitting prediction
              % Store parameters
              for p = 1:length(Model.Param_Names)
                  Pimg(p,vx(v)) = S(p,rx); % Add the pth fitted parameter
              end
              Rimg(1,vx(v)) = mR(v);  % Variance explained
              % If running coarse fit only determine betas now
              if Model.Coarse_Fit_Only 
                  % Fit betas for amplitude & intercept
                  B = [ones(length(Y(:,v)),1) X(:,rx)] \ Y(:,v); % GLM fit 
                  Bimg(1,vx(v)) = B(2); % Amplitude
                  Bimg(2,vx(v)) = B(1); % Intercept
              end            
          end
         if wb waitbar((vs+v-1)/length(mver),h); end
      end
  end
  if wb close(h); end
  t2 = toc(t0); 
  disp(['Coarse fitting completed in ' num2str(t2/60) ' minutes.']);
end
new_line;

%% Run fine-fit or break now?
if Model.Coarse_Fit_Only 
    %% Store coarse fit parameters
    disp('Only running coarse fit!');
    Srf.X = zeros(size(Tc)); % Matrix with unconvolved predictions
    fPimg = Pimg; % Coarse fitted parameter maps
    fRimg = Rimg; % R^2 map
    fBimg = Bimg; % Beta map
    OutFile = [OutFile '_CrsFit']; % Suffix to indicate coarse fit only 
else
    %% Run fine fit for each vertex
    warning off % Because of rank deficiency warnings in GLM for artifactual vertices
    disp('Slow fine fitting...');
    if wb h = waitbar(0, ['Slow fine fitting... ' strrep(OutFile,'_','-')], 'Units', 'pixels', 'Position', [100 100 360 70]); end
    Tc = Srf.Y; % Always use unsmoothed data for fine fit
    Srf.X = zeros(size(Tc)); % Matrix with unconvolved predictions
    fPimg = zeros(length(Model.Param_Names), size(Srf.Vertices,1)); % Fine fitted parameter maps
    fRimg = zeros(1,size(Srf.Vertices,1)); % R^2 map
    fBimg = zeros(2,size(Srf.Vertices,1)); % Beta map
    % Keep track of redundancies
    Fitted = zeros(1,size(Srf.Vertices,1)); % Toggle if vertex was already analysed
    % Only loop through mask vertices
    for v = 1:length(mver)
        % Index of current vertex
        vx = mver(v);

        if Fitted(vx) == 0 && Rimg(vx) >= Model.Fine_Fit_Threshold % Only non-redundant & reasonable coarse fits
            % Find redundant vertices
            rd = samsrf_find_redundancy(Tc,vx);  
            % Mark all redundant vertices
            Fitted(rd) = 1; 

            % Find best prediction
            Y = Tc(:,vx);  % Time course of current vertex
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y), ...
                        Pimg(:,vx)', optimset('TolX',1e-2,'TolFun',1e-2,'Display','off'));  % Lower than default tolerance             

            % Loop thru fitted parameters
            IsGoodFit = true;
            for p = 1:length(fP)
                % If scaled parameter is out of bounds
                if Model.Scaled_Param(p) && abs(fP(p)) > 3
                    IsGoodFit = false;
                end
                % If a parameter is negative but shouldn't be
                if Model.Only_Positive(p) && fP(p) < 0
                    IsGoodFit = false;
                end
            end

            % Only keep good parameters
            if IsGoodFit
                fR = 1 - fR; % 1 minus unexplained variance
                % Generate predicted time course
                Rfp = Model.Prf_Function(fP, size(ApFrm,1)*2);
                fX = prf_predict_timecourse(Rfp, ApFrm, false);
                % Store prediction without HRF convolution
                Srf.X(:,rd) = repmat(fX, 1, length(rd));  % Best fitting prediction
                % Loop thru fitted parameters
                for p = 1:length(fP)
                    fPimg(p,rd) = fP(p); % Add the pth fitted parameter
                end
                % Convolve predicted time course with HRF
                fX = conv(fX, Model.Hrf);
                fX = fX(1:length(Y)); % Truncate back to original length
                % Fit betas for amplitude & intercept
                fB = [ones(length(Y),1) fX] \ Y; % GLM fit 
                fBimg(1,rd) = fB(2); % Amplitude
                fBimg(2,rd) = fB(1); % Intercept
                fRimg(rd) = fR;  % Variance explained
            else
                % Replace bad fits with coarse fit?
                if Model.Replace_Bad_Fits
                    fR = Rimg(vx); % Goodness of coarse fit
                    fP = Pimg(:,vx); % Parameter estimates from coarse fit
                    % Generate predicted time course
                    Rfp = Model.Prf_Function(fP, size(ApFrm,1)*2);
                    fX = prf_predict_timecourse(Rfp, ApFrm, false);
                    % Store prediction without HRF convolution
                    Srf.X(:,rd) = repmat(fX, 1, length(rd));  % Best fitting prediction
                    % Loop thru fitted parameters
                    for p = 1:length(fP)
                        fPimg(p,rd) = fP(p); % Add the pth fitted parameter
                    end
                    % Convolve predicted time course with HRF
                    fX = conv(fX, Model.Hrf);
                    fX = fX(1:length(Y)); % Truncate back to original length
                    % Fit betas for amplitude & intercept
                    fB = [ones(length(Y),1) fX] \ Y; % GLM fit 
                    fBimg(1,rd) = fB(2); % Amplitude
                    fBimg(2,rd) = fB(1); % Intercept
                    fRimg(rd) = fR;  % Variance explained
                end
            end
        end
        if wb waitbar(v/length(mver),h); end
    end
    if wb close(h); end
    warning on
    new_line;
end

% Rescale the data
for i = 1:size(fPimg,1)
    if Model.Scaled_Param(i)
        fPimg(i,:) = fPimg(i,:) * Model.Scaling_Factor;  % Scale this parameter
    end
end

% Prepare surface structure
Srf.Functional = PrfFcnName; % pRF function name
Srf.Data = [fRimg; fPimg; fBimg]; % Parameter maps
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
end
Srf.Values = Srf.Values'; % So it is the same as Srf.Data

% Compress to save space
Srf = samsrf_compress_srf(Srf, mver);
% Save fine fit map files
disp('Saving pRF fitting results...');
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t3 = toc(t0); 
EndTime = num2str(t3/60/60);
new_line; disp(['Whole analysis completed in ' EndTime ' hours.']);
