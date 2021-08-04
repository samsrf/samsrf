function OutFile = samsrf_revcor_prf(Model, SrfFiles, Roi)
%
% OutFile = samsrf_revcor_prf(Model, SrfFiles, [Roi])
%
% Calculates the pRF profiles using a reverse correlation technique:
% It first creates a regressor for each pixel in the stimulus apertures
% convolved with the HRF and then correlates these with the observed time
% series. This produces a map of correlation coefficients that is then
% thresholded, downsampled, and saved. A standard retinotopic map with
% X0 and Y0 coordinates (peak correlation position) is also saved.
%
% Note that this procedure does not estimate pRF size. In order to do this,
% you will need to run a post-processing function to fit pRF models to the
% reverse correlation profiles.
%
%   Model:          Contains the parameters defining the pRF model,
%                       the search space for the coarse fitting,
%                       and the eccentricity/scaling factor of the space.
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
% 18/07/2020 - SamSrf 7 version (DSS)
% 27/07/2020 - More info in command window (DSS)
% 24/07/2020 - Added option to limit data by noise ceiling (DSS)
%              Reorganised analysis loop but parallel processing isn't working yet (DSS)  
% 29/03/2021 - Fixed horrendous bug when fitting bilateral surface meshes (DSS)    
% 24/05/2021 - Displays asterisks & new lines when analysis is complete (DSS)
% 30/06/2021 - Added new-fangled old-school command-line progress-bars (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end
if ~isfield(Model, 'Noise_Ceiling_Threshold')
    Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
end

%% Start time of analysis
t0 = tic; new_line;  
disp('*** SamSrf reverse-correlation pRF analysis ***');
[vn, vd] = samsrf_version;
disp([' Version ' num2str(vn) ' - ' vd]);
new_line;
disp('Current working directory:');
disp([' ' pwd]);
new_line;

%% Load apertures
disp('Load stimulus apertures...');
load(Model.Aperture_File);  % Loads a variable called ApFrm
disp([' Loading '  Model.Aperture_File ': ' num2str(size(ApFrm,3)) ' volumes']);
new_line; 
% Regressor file name
[~,ApsName] = fileparts(Model.Aperture_File);
RegressorFile = ['reg_' ApsName '.mat'];  

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

%% Generate regressors
if ~exist([pwd filesep RegressorFile], 'file') 
    disp('Generating regressors...');
    % Reshape apertures to matrix
    X = reshape(ApFrm, size(ApFrm,1)*size(ApFrm,2), size(ApFrm,3))';
    % Convolve with HRF
    for i = 1:size(X,2) 
        X(:,i) = prf_convolve_hrf(X(:,i),Model.Hrf); 
    end
    % X & Y coordinates
    [xc,yc] = meshgrid(linspace(-Model.Scaling_Factor, Model.Scaling_Factor, size(ApFrm,2)), ...
                       linspace(-Model.Scaling_Factor, Model.Scaling_Factor, size(ApFrm,1)));
    yc = flipud(yc);  % Flip vertically to be in Cartesian space
    % Save regressors
    save(RegressorFile, 'X', 'xc', 'yc', '-v7.3');
    t1 = toc(t0); 
    disp([' Regressors generated in ' num2str(t1) ' seconds.']);
else
    disp('Loading regressors...');
    load([pwd filesep RegressorFile]);
    disp([' Loading ' RegressorFile]);
    % Does length of regressors match apertures?
    if size(ApFrm,3) ~= size(X,1)
        error('Mismatch between length of search space and apertures!');
    end
end
new_line; 

%% Calculate reverse correlation map 
disp('Calculating pRF profiles...');
Srf.Regs = X;  % Design matrix (regressors per pixel)
% R-map for each vertex 
Rmaps = zeros(Model.Rdim^2, length(mver));  
% Parameters for pRFs
fXimg = zeros(1,length(mver)); % X-coordinate map
fYimg = zeros(1,length(mver)); % Y-coordinate map
fSimg = zeros(1,length(mver)); % Sigma map
fBimg = zeros(1,length(mver)); % Beta map
fRimg = zeros(1,length(mver)); % R^2 map
% Loop through mask vertices 
samsrf_progbar(0);
for v = 1:length(mver)
    % Calculate r-map
    Y = Tc(:,mver(v));  % Time course of current vertex
    warning off
    cM = [Y ones(size(Y,1),1)] \ X; % Linear regression
    warning on
    cM = cM(1,:); % Remove intercept beta
    mM = max(cM); % Find peak activation in each map
    mM = mM(1); % Ensure only one value
    gp = cM > mM/2; % Full area at half maximum
    m = find(cM==mM,1); % Find peak coordinate
    mR = corr(Y, X(:,m)); % Correlation coefficient at peak
    if ~isempty(mR) && ~isnan(mR) && mR >= 0.1
        cM = reshape(cM,size(ApFrm,1),size(ApFrm,2)); % Reshape into a map
        cM = imresize(cM,[Model.Rdim Model.Rdim]); % Down-sample r-map
        cM = cM(:); % Vectorise again
        % Store pRF profile
        Rmaps(:,v) = cM; % Activation map as vector 
        fXimg(v) = xc(m);  % X-coordinate
        fYimg(v) = yc(m);  % Y-coordinate
        fSimg(v) = sqrt(mean(gp) * (Model.Scaling_Factor*2)^2); % Full width at half maximum
        fBimg(v) = mM;  % Activation peak
        fRimg(v) = mR^2;  % Variance explained
    end
    % Progress report
    samsrf_progbar(v/length(mver));
end
t2 = toc(t0); 
disp(['Analysis completed in ' num2str(t2/60) ' minutes.']);

% Coordinates for contour/surf plots
Xc = imresize(xc,[Model.Rdim Model.Rdim]); 
Yc = imresize(yc,[Model.Rdim Model.Rdim]); 
Srf.X_coords = Xc;
Srf.Y_coords = Yc;

% Save as surface structure
Srf.Functional = 'Reverse correlation';
Srf.Data = zeros(5, size(Srf.Vertices,1));
Srf.Data(:,mver) = [fRimg; fXimg; fYimg; fSimg; fBimg];
Srf.Values = {'R^2'; 'x0'; 'y0'; 'Fwhm'; 'Beta'};
Srf.Rmaps = zeros(Model.Rdim^2, size(Srf.Vertices,1));
Srf.Rmaps(:,mver) = Rmaps; % Add activation maps 
% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
end

% Save map files
disp('Saving pRF results...');
Srf = samsrf_compress_srf(Srf, mver);
OutFile = [OutFile '_Rcp'];
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t3 = toc(t0); 
EndTime = num2str(t3/60);
new_line; disp(['Whole analysis completed in ' EndTime ' minutes.']);
disp('******************************************************************');
new_line; new_line;

