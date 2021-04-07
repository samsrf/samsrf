function OutFile = samsrf_fit_prf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_fit_prf(Model, SrfFiles, [Roi = ''])
%
% Fits a pRF model using the slow fine-fitting procedure (fminsearch in Optimization toolbox).
%
% IMPORTANT NOTE: In SamSrf 7 the precision of the model fitting algorithm was improved!
%                 Model fits using this version are -NOT- compatible with previous versions!
%                 Under some conditions previous versions sometimes reverted to the search grid.
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
% 20/07/2020 - SamSrf 7 version (DSS) 
% 23/07/2020 - Cosmetic changes to command window outputs (DSS)
% 03/02/2021 - Fixed crashing bug when replacing bad fine fits (DSS)
% 05/02/2021 - Fixed another smaller bug with time series when replacing bad fine fits (DSS)  
% 07/04/2021 - Added parameter option for only allowing positive coarse fits to pass (DSS)   
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
if ~isfield(Model, 'Noise_Ceiling_Threshold')
    Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
end
if ~isfield(Model, 'Polar_Search_Space')
    Model.Polar_Search_Space = false; % Search space is Cartesian
end
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
if ~isfield(Model, 'Only_Positive_Coarse_Fits')
    Model.Only_Positive_Coarse_Fits = false; % Coarse fit can either be negative or positive correlation to pass 
end
if ~isfield(Model, 'Coarse_Fit_Block_Size')
    Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
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

%% Load apertures
disp('Load stimulus apertures...');
load(Model.Aperture_File);  % Loads a variable called ApFrm
disp([' Loading ' Model.Aperture_File ': ' num2str(size(ApFrm,3)) ' volumes']);
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
        disp(['Limiting analysis to ' num2str(size(mver,1)) ' vertices above noise ceiling ' num2str(Model.Noise_Ceiling_Threshold)]);
        new_line;
    end
end

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
        [X,S] = prf_generate_searchspace(Model.Prf_Function, ApFrm, Model.Param1, Model.Param2, Model.Param3, Model.Param4, Model.Param5, Model.Polar_Search_Space);    
        save(SearchspaceFile, 'X', 'S', '-v7.3');
        t1 = toc(t0); 
        disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
    else
        disp('Loading pre-defined predictions...');
        load([pwd filesep SearchspaceFile]);
        disp([' Loading ' SearchspaceFile]);
        % Does number of grid points match model?
        if size(S,2) ~= length(Model.Param1) * length(Model.Param2) * length(Model.Param3) * length(Model.Param4) * length(Model.Param5)
            error('Mismatch between saved search space and model definition!');
        end
        % Does length of search space match apertures?
        if size(ApFrm,3) ~= size(X,1)
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
    for p = 1:size(X,2)
        X(:,p) = prf_convolve_hrf(X(:,p), Model.Hrf);
    end
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
  Srf.X = zeros(size(Tc));  % Matrix with predictions
  Pimg = zeros(length(Model.Param_Names), size(Srf.Vertices,1)); % Fitted parameter maps
  Rimg = zeros(1, size(Srf.Vertices,1)); % R^2 map
  % If only running coarse fit
  if Model.Coarse_Fit_Only
      Bimg = zeros(2,size(Srf.Vertices,1)); % Beta map
  end
  
  % Loop through mask vertices (in blocks if Matlab R2012a or higher)
  disp([' Block size: ' num2str(cfvb) ' vertices']);
  for vs = 1:cfvb:length(mver)
      % Starting index of current vertex block
      ve = vs + cfvb - 1;
      if ve > length(mver)
          ve = length(mver);
      end
      vx = mver(vs:ve);
      % Find best prediction
      Y = Tc(:,vx);  % Time course of current vertex
      if Parameters.Only_Positive_Coarse_Fits
      	 R = corr(Y,X); % Best correlating prediction 
      	 mR = max(R,[],2); % Find best fit
      	 R = R.^2; % Now turn into R^2
      else
      	 R = corr(Y,X).^2; % Best correlating prediction (squared to allow for negative betas!)
      	 mR = max(R,[],2); % Find best fit
      end
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
                warning off % In case of rank deficient GLM
                B = [ones(length(Y(:,v)),1) X(:,rx)] \ Y(:,v); % GLM fit 
                warning on
                Bimg(1,vx(v)) = B(2); % Amplitude
                Bimg(2,vx(v)) = B(1); % Intercept
              end            
          end
          if cfvb > 1 && v == length(vx)
              disp([' ' num2str(round(ve/length(mver)*100)) '% completed']);
          end
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
    Srf.X = zeros(size(Tc)); % Matrix with unconvolved predictions
    fPimg = Pimg(:,mver); % Coarse fitted parameter maps
    fRimg = Rimg(1,mver); % R^2 map
    fBimg = Bimg(:,mver); % Beta map
    OutFile = [OutFile '_CrsFit']; % Suffix to indicate coarse fit only 
else
    %% Run fine fit for each vertex
    disp('Fine fitting...');
    % Restrict matrices to mask
    Tc = Srf.Y(:,mver); % Always use unsmoothed data for fine fit
    Rimg = Rimg(1,mver);
    Pimg = Pimg(:,mver);
    
    % Run fine fit (must be separate function for parallel computing toolbox
    [fPimg, fRimg] = samsrf_fminsearch_loop(Model, Tc, ApFrm, Rimg, Pimg);
    t3 = toc(t0);
    disp(['Fine fitting completed in ' num2str(t3/60/60) ' hours.']);
    new_line;
    
    disp('Fitting beta parameters & storing fitted models...');   
    % Additional data fields
    fBimg = zeros(2,length(mver)); % Beta maps
    Srf.X = zeros(size(Tc,1), size(Srf.Vertices,1)); % Matrix with unconvolved predictions
    
    % Process & fit betas for mask vertices
    for v = 1:length(mver)
        % Loop thru fitted parameters
        IsGoodFit = true;
        for p = 1:size(fPimg,1)
            % If scaled parameter is out of bounds
            if Model.Scaled_Param(p) && abs(fPimg(p,v)) > 2 
                IsGoodFit = false;
            end
            % If a parameter is negative but shouldn't be
            if Model.Only_Positive(p) && fPimg(p,v) < 0
                IsGoodFit = false;
            end
        end

        % Store parameters & predicted time course 
        if IsGoodFit 
            % Only keep good parameters
            Rfp = Model.Prf_Function(fPimg(:,v)', size(ApFrm,1)*2); % Generate predicted pRF
            fX = prf_predict_timecourse(Rfp, ApFrm); % Generate predicted time course
            % Store prediction without HRF convolution
            Srf.X(:,mver(v)) = fX;  % Best fitting prediction
        else
            % Replace bad fits with coarse fit?
            if Model.Replace_Bad_Fits
                % Generate predicted time course
                Rfp = Model.Prf_Function(Pimg(:,v)', size(ApFrm,1)*2);
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
            % Convolve predicted time course with HRF
            fX = prf_convolve_hrf(fX, Model.Hrf);
            % Fit betas for amplitude & intercept
            warning off % In case of rank deficient GLM
            fB = [ones(length(Y),1) fX] \ Y; % GLM fit
            warning on
            fBimg(:,v) = fB([2 1]); % Amplitude & intercept
        end
    end
end

disp('Tidying up final results structure...');
% Rescale the data
for i = 1:size(fPimg,1)
    if Model.Scaled_Param(i)
        fPimg(i,:) = fPimg(i,:) * Model.Scaling_Factor;  % Scale this parameter
    end
end

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
