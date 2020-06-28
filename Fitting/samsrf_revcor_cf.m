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
% You must be in the folder containing the surface data files as well as
% the various parameter files (i.e. apertures, searchspace, and HRF).
%
% Returns the name of the map file it saved.
%
% 09/08/2018 - SamSrf 6 version (DSS)
% 05/05/2020 - Added option to smooth connectivity profiles (DSS)
% 27/05/2020 - Streamlined how waitbar is handled (DSS)
%

%% Defaults & constants
% If no ROI defined, analyse whole brain...
if nargin < 3
    Roi = ''; 
end

%% Start time of analysis
t0 = tic; new_line;  
disp('*** Connective field analysis ***');

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

%% Load seed ROI 
svx = samsrf_loadlabel(Model.SeedRoi);
Srf.SeedVx = svx; % Store for posterity
X = Tc(:,svx); % Time courses in seed ROI

%% Calculate cortical distances within seed ROI
if Model.Smoothing > 0 % No point when not smoothing
    Ds = samsrf_geomatrix(Srf.Vertices, Srf.Faces, svx); % Cortical distance matrix
    Ws = exp(-(Ds.^2)/(2*Model.Smoothing.^2)); % Smoothing weight matrix
end

%% Load template map
Temp = load(Model.Template);
Temp.Srf = samsrf_expand_srf(Temp.Srf);

%% Add version number
Srf.Version = samsrf_version;

%% Preprocess data
% Store raw time courses
Srf.Y = Tc; % Raw time coarse stored away
Srf.Data = NaN(length(svx),size(Srf.Vertices,1)); % Connective field profile for each vertex 

%% Calculate reverse correlation map 
disp('Calculating CF profiles...');
% Keep track of redundancies
Fitted = zeros(1,size(Srf.Vertices,1));   % Toggle if vertex was already analysed
% Loop through mask vertices (in blocks if Matlab R2012a or higher)
for v = 1:length(mver)
    % Index of current vertex
    vx = mver(v);

    if Fitted(vx) == 0 % Only non-redundant vertices    
        % Calculate r-map
        Y = Tc(:,vx);  % Time course of current vertex
        R = corr(Y,X); % Correlation of time courses with regressors
        % Smooth profile if needed
        if Model.Smoothing > 0
            sR = repmat(R,length(svx),1) .* Ws; % Smoothed profiles for each seed vertex 
            R = sum(sR,2) ./ sum(Ws,2); % Smoothed seed ROI    
        end
        % Determine parameters
        mR = max(R); % Find peak activation in each map
        mR = mR(1); % Ensure only one value
        rd = samsrf_find_redundancy(Tc,vx); % Find redundant vertices
        if ~isnan(mR) && mR > 0
            Fitted(rd) = 1; % Mark all redundant vertices
            % Store CF profile
            Srf.Data(:,rd) = repmat(R(:), 1, length(rd)); % Activation map as vector 
        else
            Fitted(rd) = 1; % Mark all redundant vertices
        end
    end
end
t2 = toc(t0); 
disp(['Correlation analysis completed in ' num2str(t2/60) ' minutes.']);

%% Determine CF parameters
% Keep track of redundancies
Fitted(:) = 0;  % Toggle if vertex was already analysed
% CF parameters
fVimg = zeros(1,size(Srf.Vertices,1)); % Peak vertex index
fXimg = zeros(1,size(Srf.Vertices,1)); % Left-Right coordinate
fZimg = zeros(1,size(Srf.Vertices,1)); % Inferior-Superior coordinate
fWimg = zeros(1,size(Srf.Vertices,1)); % Full width half maximum
fRimg = zeros(1,size(Srf.Vertices,1)); % R^2 map
% % Smooth connective field profiles if desired
if Model.Smoothing > 0   
    % CURRENTLY NOT YET IMPLEMENTED!
end
Srf.ConFlds = Srf.Data; % Smoothed correlation profile
Srf.Data = [];
disp('Estimating CF parameters...');
% Keep track of redundancies
Fitted = zeros(1,size(Srf.Vertices,1));   % Toggle if vertex was already analysed
% Loop through mask vertices (in blocks if Matlab R2012a or higher)
for v = 1:length(mver)
    % Index of current vertex
    vx = mver(v);

    % Retrieve r-map
    R = Srf.ConFlds(:,vx); % Correlation of time courses with regressors
    % Determine parameters
    mR = max(R); % Find peak activation in each map
    mR = mR(1); % Ensure only one value
    fwhm = sqrt(sum(Srf.Area(R > mR/2))); % Square root of area above half maximum
    m = find(R==mR,1); % Find peak coordinate
    rd = samsrf_find_redundancy(Tc,vx); % Find redundant vertices 
    if ~isnan(mR) && mR > 0
        Fitted(rd) = 1; % Mark all redundant vertices
        fVimg(rd) = svx(m);  % Peak vertex in seed ROI
        fXimg(rd) = Temp.Srf.Data(2,svx(m)); % Left-Right coordinate
        fZimg(rd) = Temp.Srf.Data(3,svx(m)); % Inferior-Superior coordinate
        fWimg(rd) = fwhm; % Full width half maximum guestimate
        fRimg(rd) = mR^2;  % Peak correlation squared 
    else
        Fitted(rd) = 1; % Mark all redundant vertices
    end
end
t3 = toc(t0); 
disp(['Parameter estimates completed in ' num2str(t3/60) ' minutes.']);

% Save as surface structure
Srf.Functional = 'Connective field';
Srf.Data = [fRimg; fXimg; fZimg; fWimg; fVimg];
Srf.Values = {'R^2'; 'x0'; 'y0'; 'Fwhm'; 'Vx'};

% Save map files
disp('Saving CF results...');
Srf = samsrf_compress_srf(Srf, mver);
save(OutFile, 'Model', 'Srf', '-v7.3');
disp(['Saved ' OutFile '.mat']); 

% End time
t4 = toc(t0); 
EndTime = num2str(t4/60);
new_line; disp(['Whole analysis completed in ' EndTime ' minutes.']);
