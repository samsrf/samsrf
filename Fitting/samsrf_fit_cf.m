function OutFile = samsrf_fit_cf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_fit_cf(Model, SrfFiles, [Roi = ''])
%
% Fits a CF model using a fast grid search fitting procedure.
% For each vertex in the region of interest Roi, determines the best CF parameters
% for which the time series in the seed region Model.SeedRoi can predict its response.
%
%   Model:          Contains the parameters defining the connective field model.
%   SrfFiles:       Cell array of file names without extension (e.g. {'lh_Bars1' 'lh_Bars2'})
%                       Files will be concatenated in this order.
%   Roi:            ROI label to restrict the analysis (default = '') 
%                       Optional, but without this the analysis can take forever.
%
% Returns the name of the map file it saved.
%
% 22/05/2021 - Writen (DSS) 
% 24/05/2021 - Displays asterisks & new lines when analysis is complete (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end
% Define search space name
SearchspaceFile = ['src_' Model.Name '.mat'];  

%% Default model parameters
if ~isfield(Model, 'Noise_Ceiling_Threshold')
    Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
end
if ~isfield(Model, 'Smoothing')
    Model.Smoothing = 0; % No smoothing on coarse fit
end
if ~isfield(Model, 'Only_Positive_Coarse_Fits')
    Model.Only_Positive_Coarse_Fits = false; % Coarse fit can either be negative or positive correlation to pass 
end
if ~isfield(Model, 'Coarse_Fit_Block_Size')
    Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
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
disp('*** SamSrf CF model fitting ***');
[vn, vd] = samsrf_version;
disp([' Version ' num2str(vn) ' - ' vd]);
new_line;
disp('Current working directory:');
disp([' ' pwd]);
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
Srf.Data = Tc; % Store full time course in Srf

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

%% Load seed ROI 
disp(['Loading seed ROI: ' Model.SeedRoi]);
svx = samsrf_loadlabel(Model.SeedRoi);
Srf.SeedVx = svx; % Store for posterity

%% Load template map
disp(['Loading template map: ' Model.Template]);
Temp = load(Model.Template);
Temp.Srf = samsrf_expand_srf(Temp.Srf);
new_line;

%% Add version number
Srf.Version = samsrf_version;

%% Generate prediction matrix
if ~exist([pwd filesep SearchspaceFile], 'file') 
    disp('Generating predictions...');
    [X,S] = cf_generate_searchspace(Srf, svx, Model.Sizes);    
    save(SearchspaceFile, 'X', 'S', '-v7.3');
    t1 = toc(t0); 
    disp([' Search space generated in ' num2str(t1/60) ' minutes.']);
else
    disp('Loading pre-defined predictions...');
    load([pwd filesep SearchspaceFile]);
    disp([' Loading ' SearchspaceFile]);
    % Does number of grid points match model?
    if size(S,2) ~= length(svx) * length(Model.Sizes)
        error('Mismatch between saved search space and model definition!');
    end
    % Does length of search parameter & prediction matrix match?
    if size(S,2) ~= size(X,2)
        % Shouldn't ever happen unless someone screwed with the search space file
        error('Search space is corrupt! Mismatch between number of parameters & predictions!');
    end
end
disp(['Using search space with ' num2str(size(S,2)) ' grid points.']);
new_line; 

%% Smooth for coarse fit?
if Model.Smoothing > 0
    Srf = samsrf_smooth_sphere(Srf, Model.Smoothing, Roi, 0); % Smooth within ROI
    Srf = rmfield(Srf, 'Raw_Data'); % Remove raw data field
    Tc = Srf.Data; % Put smoothed data back into time course variable
end

%% Preprocess data
% Store raw time courses
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = [];  % Clear data field

%% Coarse fit 
disp('Coarse fitting...');
Srf.X = zeros(size(Tc));  % Matrix with predictions
Vimg = zeros(1, size(Srf.Vertices,1)); % Fitted vertex number map
Simg = zeros(1, size(Srf.Vertices,1)); % Fitted size parameter map
Pimg = zeros(2, size(Srf.Vertices,1)); % Template parameter maps
Rimg = zeros(1, size(Srf.Vertices,1)); % R^2 map
Bimg = zeros(2,size(Srf.Vertices,1)); % Beta map
  
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
  if Model.Only_Positive_Coarse_Fits
     R = corr(Y,X); % Mean corrected correlation 
     mR = max(R,[],2); % Find best fit
     R = R.^2; % Now turn into R^2
  else
      R = corr(Y,X).^2; % Mean corrected correlation (squared to allow for negative betas!)
      mR = max(R,[],2); % Find best fit
  end
  for v = 1:length(vx)
      rx = find(R(v,:) == mR(v)); % Matrix position of best prediction
      if ~isempty(rx)
          rx = rx(1); % Only first instance 
          % Store prediction
          Srf.X(:,vx(v)) = X(:,rx);  % Best fitting prediction
          % Store parameters
          Vimg(1,vx(v)) = S(1,rx); % Add vertex number
          Simg(1,vx(v)) = S(2,rx) / 3; % Add size parameter (rescaled to sigma)
          Rimg(1,vx(v)) = mR(v);  % Variance explained
          Pimg(:,vx(v)) = Temp.Srf.Data(2:3, S(1,rx)); % Retrieve template coordinates 
          % Fit betas for amplitude & intercept
          warning off % In case of rank deficient GLM
          B = [ones(length(Y(:,v)),1) X(:,rx)] \ Y(:,v); % GLM fit 
          warning on
          Bimg(1,vx(v)) = B(2); % Amplitude
          Bimg(2,vx(v)) = B(1); % Intercept
      end
      if cfvb > 1 && v == length(vx)
          disp([' ' num2str(round(ve/length(mver)*100)) '% completed']);
      end
  end
end
t2 = toc(t0); 
disp(['Coarse fitting completed in ' num2str(t2/60) ' minutes.']);
new_line;

%% Store coarse fit parameters
disp('Only running coarse fit!');
Srf.X = zeros(size(Tc)); % Matrix with unconvolved predictions
fPimg = Pimg(:,mver); % Template parameter maps
fRimg = Rimg(1,mver); % R^2 map
fVimg = Vimg(1,mver); % Vertex number map
fSimg = Simg(1,mver); % CF size map
fBimg = Bimg(:,mver); % Beta map

% Prepare surface structure
Srf.Functional = 'Connective Field Fast-Fit'; % pRF function name
Data = [fRimg; fPimg; fSimg; fBimg; fVimg];
Srf.Data = zeros(size(Data,1), size(Srf.Vertices,1));
Srf.Data(:,mver) = Data; % Add parameter maps into full data matrix
% Add parameter names
Srf.Values = {'R^2'; 'x0'; 'y0'; 'Sigma'; 'Beta'; 'Baseline'; 'Vx'}; 
% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
    Srf = samsrf_normr2(Srf); % Calculate normalised R^2
end
new_line;

% Compress to save space
Srf = samsrf_compress_srf(Srf, mver);
% Save fine fit map files
disp('Saving CF fitting results...');
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60/60);
new_line; disp(['Whole analysis completed in ' EndTime ' hours.']);
disp('******************************************************************');
new_line; new_line;

