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
% Note that by default this procedure does not estimate pRF size, only a 
% guestimate based on the area the pRF profile subtends in visual space. 
% In order to estimate a pRF size, you will need to fit a 2D pRF model to 
% the reverse correlation profiles (Model.Prf_Function).
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
% IMPORTANT: Results are not completely identical to those from SamSrf 7!
%            While parameter estimates are mostly very similar, there are 
%            some small differences (& also larger ones especially for 
%            peripheral pRFs). Goodness-of-fit (R^2_{pRF}) is subtly but
%            (extremely) significantly greater; however, some voxels that 
%            surpassed thresholding in SamSrf 7 are rejected in SamSrf 8. 
%            
%            The reason for this is that SamSrf 7 used the downscaled
%            correlation profiles for 2D model fitting while SamSrf 8 
%            recomputes the full profiles before fitting. Since this seems 
%            objectively more appropriate, we decided not to add full 
%            compatibility as an option in SamSrf 8. 
%            (Obviously, this issue only relates to the 2D fits but the 
%             reverse correlation results should be perfectly identical)
%
% 14/03/2022 - Added option to remove reverse correlation profiles (DSS)
%              Now also saves the noise ceiling if it exists in raw data file (DSS)
% 13/04/2022 - Now checks that vectors defining parameters are all same length (DSS)
%              Added Hooke-Jeeves algorithm & customisable Nelder-Mead tolerance (DSS)
% 14/04/2022 - pRF model fitting now uses parallel computing (DSS)
% 15/04/2022 - Warns if both Hooke-Jeeves steps & Nelder-Mead tolerance are defined (DSS)
%              Outsourced check for default parameters so no longer needs to check these (DSS)
%              Final duration now reported in hours (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
%              Added explanation about discrepancies with SamSrf 7 fits (DSS)
% 23/06/2022 - Now warns if aperture matrix contains negative values (DSS)
% 10/08/2022 - Removed initial thresholding of reverse correlation profiles (DSS)
%              Now supports convex hull algorithm for fast estimation of pRF parameters (DSS)
% 31/08/2023 - Added option for negative peaks (DSS)
%              Removed unnecessary line (DSS)
% 01/09/2023 - Beta now reflects negative peaks (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end

%% Default model parameters
Model = samsrf_model_defaults('samsrf_revcor_prf', Model);

%% Start time of analysis
t0 = tic; new_line;  
disp('*** SamSrf reverse-correlation pRF analysis ***');
[vn, vd] = samsrf_version;
disp([' Version ' num2str(vn) ' - ' vd]);
new_line;
disp('Current working directory:');
disp([' ' pwd]);
new_line;
% Are we also fitting pRF model?
if isfield(Model, 'Prf_Function')
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
    new_line;
end

%% Load apertures
disp('Load stimulus apertures...');
load(EnsurePath(Model.Aperture_File));  % Loads a variable called ApFrm
disp([' Loading '  Model.Aperture_File ': ' num2str(size(ApFrm,3)) ' volumes']);
if sum(ApFrm(:)<0) > 0
    disp(' Warning: Apertures contain negative values!');
end
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
        X(:,i) = prf_convolve_hrf(X(:,i), Model.Hrf, 1); % Added 3rd input to ensure it works but dowwnsampling has not yet been implemented!
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
    % If peak can be negative
    PeakSign = 1;
    if Model.Allow_Negative_Peaks
        % Invert sign if peak is negative
        if abs(max(cM)) < abs(min(cM))
            cM = -cM;
            PeakSign = -1;
        end
    end
    % Determine peak
    mM = max(cM); % Find positive peak in each map
    gp = cM > mM/2; % Full area at half maximum
    m = find(cM==mM,1); % Find peak coordinate
    mR = corr(Y, X(:,m)); % Peak correlation with stimulus design
    % Create profile & store parameters
    if ~isempty(mR) && ~isnan(mR) 
        cM = reshape(cM,size(ApFrm,1),size(ApFrm,2)); % Reshape into a map
        cM = imresize(cM,[Model.Rdim Model.Rdim]); % Down-sample r-map
        cM = cM(:); % Vectorise again
        % Store pRF parameters
        Rmaps(:,v) = cM; % Activation map as vector 
        fXimg(v) = xc(m);  % X-coordinate
        fYimg(v) = yc(m);  % Y-coordinate
        fSimg(v) = sqrt(mean(gp) * (Model.Scaling_Factor*2)^2); % Full width at half maximum
        fBimg(v) = mM * PeakSign;  % Activation peak & inverted sign if necessary
        fRimg(v) = mR^2;  % Variance explained
    end
    % Progress report
    samsrf_progbar(v/length(mver));
end
t2 = toc(t0); 
disp(['pRF profiles completed in ' num2str(t2/60) ' minutes.']);
new_line;

% Coordinates for contour/surf plots
Srf.X_coords = imresize(xc,[Model.Rdim Model.Rdim]);
Srf.Y_coords = imresize(yc,[Model.Rdim Model.Rdim]);

% Save as surface structure
Srf.Functional = 'Reverse correlation';
Srf.Data = zeros(5, size(Srf.Vertices,1));
Srf.Data(:,mver) = [fRimg; fXimg; fYimg; fSimg; fBimg];
Srf.Values = {'R^2'; 'x0'; 'y0'; 'Fwhm'; 'Beta'};
Srf.Rmaps = NaN; % If fitting pRF models, need to calculate profiles anew

%% Fit 2D pRF models?
if isfield(Model, 'Prf_Function')
    disp('Fitting pRF models to reverse correlation profiles...');
    Srf.Raw_Data = Srf.Data; % Store reverse correlations in raw data
    Srf.Raw_Values = Srf.Values; % Also store value names cause they'll change
    Srf.Values = {}; % Clear value field
    Srf.Values{1} = 'R^2'; % Goodness of model fit
    Srf.Values(2:length(Model.Param_Names)+1) = Model.Param_Names; % pRF parameters
    Srf.Values{end+1} = 'Beta'; % Amplitude 
    Srf.Values{end+1} = 'Baseline'; % Baseline
    Srf.Values = Srf.Values'; % Ensure not row vector
    % Initialise data field
    Srf.Data = NaN(length(Srf.Values), size(Srf.Raw_Data,2));
    
    % Which reverse correlation data to include?
    disp([' Using reverse correlations profiles with R^2 > ' num2str(Model.R2_Threshold)]);
    GoF = find(Srf.Raw_Data(1,:) > Model.R2_Threshold); % Vertices with good reverse correlation profiles 
   
    disp('Seeding parameters using:');
    disp(Model.SeedPar_Function);
    
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
    
    % Fit 2D models
    Srf = samsrf_revcorprf_loop(Srf, GoF, Model, ApFrm, AlgorithmParam);
    t3 = toc(t0); 
    disp(['pRF parameter fitting completed in ' num2str(t3/60/60) ' hours.']);
else
    % Using convex hull algorithm to determine pRF 
    disp('Using convex hull algorithm to estimate parameters...');
    % Loop through mask vertices 
    samsrf_progbar(0);
    for v = 1:length(mver)
        % Convex hull estimation
        vx = mver(v);
        Rmap = prf_contour(Srf, vx);
        txy = [xc(:) yc(:)]; % Matrix with pixel coordinates in visual space
        txy = txy(Rmap(:) > max(Rmap(:))/2,:); % Limit pixel matrix to suprathreshold pixels
        Dt = delaunayTriangulation(txy); % Delaunay triangulation of CF profile
        try
            % Determine pRF parameters
            [chp, cha] = convexHull(Dt); % Points & area of convex hull
            warning off
            ps = polyshape(Dt.Points(chp,:)); % Polygon of convex hull
            warning on
            [chx, chy] = ps.centroid; % Centroid of convex hull
            Srf.Data(2,vx) = chx; % X-coordinate
            Srf.Data(3,vx) = chy; % Y-coordinate
            Srf.Data(4,vx) = sqrt(cha) / (2*sqrt(2*log(2))); % Visual FWHM converted to Sigma
        catch
            % Failed determining parameters
            Srf.Data(:,vx) = 0;
        end        
        % Progress report
        samsrf_progbar(v/length(mver));
    end    
    t3 = toc(t0); 
    disp(['pRF parameter estimation completed in ' num2str(t3/60/60) ' hours.']);
end
new_line;

% Add noise ceiling if it has been calculated
if isfield(Srf, 'Noise_Ceiling')
    Srf.Data = [Srf.Data; Srf.Noise_Ceiling];
    Srf.Values{end+1} = 'Noise Ceiling'; 
    Srf = rmfield(Srf, 'Noise_Ceiling'); % Remove separate field
end

% Are we saving correlation profiles?
if Model.Save_Rmaps
    disp('Saving pRF profiles in data file.');
    Srf.Rmaps = zeros(Model.Rdim^2, size(Srf.Vertices,1));
    Srf.Rmaps(:,mver) = Rmaps; % Add activation profiles
else
    disp('Not saving pRF profiles...');
end

%% Save map files
disp('Saving pRF results...');
Srf = samsrf_compress_srf(Srf, mver);
OutFile = [OutFile '_Rcp'];
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60/60);
new_line; disp(['Whole analysis completed in ' EndTime ' hours.']);
disp('******************************************************************');
new_line; new_line;

