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
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
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
disp([' Loading ' Model.Aperture_File ': ' num2str(size(ApFrm,3)) ' volumes']);
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
Srf.Rmaps = zeros(Model.Rdim^2, size(Srf.Vertices,1));  
% Parameters for pRFs
fXimg = zeros(1,size(Srf.Vertices,1)); % X-coordinate map
fYimg = zeros(1,size(Srf.Vertices,1)); % Y-coordinate map
fSimg = zeros(1,size(Srf.Vertices,1)); % Sigma map
fBimg = zeros(1,size(Srf.Vertices,1)); % Beta map
fRimg = zeros(1,size(Srf.Vertices,1)); % R^2 map
% Keep track of redundancies
Fitted = zeros(1,size(Srf.Vertices,1));   % Toggle if vertex was already analysed
% Loop through mask vertices (in blocks if Matlab R2012a or higher)
for v = 1:length(mver)
    vx = mver(v); % Current vertex
    % Calculate r-map
    Y = Tc(:,vx);  % Time course of current vertex
    if Fitted(vx) == 0 % Only non-redundant vertices
        warning off
        cM = [Y ones(size(Y,1),1)] \ X; % Linear regression
        warning on
        cM = cM(1,:); % Remove intercept beta
        mM = max(cM); % Find peak activation in each map
        mM = mM(1); % Ensure only one value
        gp = cM > mM/2; % Full area at half maximum
        m = find(cM==mM,1); % Find peak coordinate
        mR = corr(Y, X(:,m)); % Correlation coefficient at peak
        rd = samsrf_find_redundancy(Tc,vx); % Find redundant vertices
        if ~isempty(mR) && ~isnan(mR) && mR >= 0.1
            cM = reshape(cM,size(ApFrm,1),size(ApFrm,2)); % Reshape into a map
            cM = imresize(cM,[Model.Rdim Model.Rdim]); % Down-sample r-map
            cM = cM(:); % Vectorise again
            Fitted(rd) = 1; % Mark all redundant vertices
            % Store pRF profile
            Srf.Rmaps(:,rd) = repmat(cM, 1, length(rd)); % Activation map as vector 
            fXimg(rd) = xc(m);  % X-coordinate
            fYimg(rd) = yc(m);  % Y-coordinate
            fSimg(rd) = sqrt(mean(gp) * (Model.Scaling_Factor*2)^2); % Full width at half maximum
            fBimg(rd) = mM;  % Activation peak
            fRimg(rd) = mR^2;  % Variance explained
        else
            Fitted(rd) = 1; % Mark all redundant vertices
        end
    end
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
Srf.Data = [fRimg; fXimg; fYimg; fSimg; fBimg];
Srf.Values = {'R^2'; 'x0'; 'y0'; 'Fwhm'; 'Beta'};
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
