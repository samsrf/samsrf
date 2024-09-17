function [fVimg, fXimg, fYimg, fWimg, fRimg, fSimg, fBimg, Rmaps] = samsrf_cfparam_loop(Tc, X, Area, SeedVx, Temp, Weights, FitPrf, SaveRmaps) 
%
% [fVimg, fXimg, fZimg, fWimg, fRimg, fSimg, fBimg, Rmaps] = samsrf_cfparam_loop(Tc, X, Area, SeedVx, Temp, Weights, FitPrf, SaveRmaps) 
%
% Internal function looping through solving the CF parameter fit for each mask vertex.
% Separated from main function for clarity & to allow progress bar.
%
%   Tc:         Time courses for all mask vertices
%   X:          Time courses for all seed vertices
%   Area:       Srf.Area from data structure 
%   SeedVx:     Seed ROI vertex indices
%   Temp:       Srf.Data from template structure
%   Weights:    Matrix of geodesic distance matrices for smoothing ([] means no smoothing)
%   FitPrf:     Toggles wehther to fit a 2D Gaussian pRF to the CF correlation profile
%               []:         Use convex hull to estimate parameters (default)
%               NaN:        Use Nelder-Mead algorithm with default tolerance
%               [NaN Tol]:  Use Nelder-Mead algorithm with parameter tolerance Tol
%               Vector:     Use Hooke-Jeeves algorithm with these initial step sizes
%               Inf:        Use simple summary statistics of significant template vertices
%   SaveRmaps:  Toggles whether Rmaps returns correlation profiles (default = false)    
%
% Returns the parameter maps & if desired, the reverse correlation profiles in Rmaps.
%
% 07/04/2022 - Fitting pRF now thresholds correlations by half-maximum (DSS)
% 08/04/2022 - Now uses iterative search to home in on pRF size (DSS)
% 12/04/2022 - Removed inconsequential erroneous comment (DSS)
% 14/04/2022 - Added Hooke-Jeeves alogrithm & adjustable Nelder-Mead tolerance (DSS)
%              Removed non-parallel computing option (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 06/05/2022 - Implemented parameter estimation based on convex hull (DSS)
% 16/05/2022 - Incorporated reverse correlation step for computational efficiency (DSS)
%              Now has option of estimating parameters using summary statistics (DSS)
%              Corrected error in help section (DSS)
% 10/08/2022 - Fixed error with computing summary statistics of CF profile (DSS)
%              Implemented region growing algorithm for better CF size estimation (DSS)
%              Added convex hull estimation of inhibitory surround (DSS)
% 07/09/2022 - Fixed bug when saving correlation profiles (DSS)
% 02/12/2023 - Bugfix when using single (32bit) data (DSS)
%

if nargin < 7
    FitPrf = NaN;
end
if nargin < 8
    SaveRmaps = false;
end

% Dimensions
nver = size(Tc,2); % Number of mask vertices
nsvx = size(X,2); % Number of seed vertices

% Output matrices
fRimg = zeros(1,nver); % R^2 map
fVimg = zeros(1,nver); % Seed vertex map
fXimg = zeros(1,nver); % X map
fYimg = zeros(1,nver); % Y map
fWimg = zeros(1,nver); % Width map
fSimg = zeros(1,nver); % Sigma map for convex hull
if isempty(FitPrf)
    fIimg = zeros(1,nver); % Sigma map
end
fBimg = zeros(2,nver); % Beta estimates 
if SaveRmaps
    Rmaps = NaN(nsvx, nver); % Connective field profile for each vertex
else
    Rmaps = NaN; % Connective field profiles are not stored
end

% Parallel processing?
try 
  gcp;
  samsrf_disp(' Parallel computing!');
catch
  error(' No parallel computing!');  
end

% Prepare progress report
try
    Q = parallel.pool.DataQueue; % Queue variable
    afterEach(Q, @percentanalysed); % Reporting function
    ProgReps = true; % Progress will be reported
catch
    Q = NaN; % No queue variable exists
    ProgReps = false; % No progress can be reported when using parallel computing
    samsrf_disp(' No progress reports possible :(');
end

%% Loop thru mask vertices 
vc = 1;
% Run parallel loop
if ProgReps
    samsrf_progbar(0);
end
parfor v = 1:nver
    % Calculate correlation profile
    Y = Tc(:,v);  % Time course of current vertex
    R = corr(Y,X); % Correlation of time courses with regressors
    % Smooth profile if needed
    if ~isempty(Weights) 
        sR = repmat(R,nsvx,1) .* Weights; % Smoothed profiles for each seed vertex 
        R = nansum(sR,2) ./ sum(Weights,2); % Smoothed seed ROI    
    end    
    % Determine parameters
    mR = max(R); % Find peak activation in each map
    mR = mR(1); % Ensure only one value
    fwhm = sqrt(sum(Area(R > mR/2))); % Square root of area above half maximum
    m = find(R==mR,1); % Find peak coordinate
    % Store correlation profile?
    if SaveRmaps
        Rmaps(:,v) = R; % Correlation profile as vector 
    end        
    % If good correlation
    if ~isnan(mR) && mR > 0
        % Parameter estimation
        fVimg(v) = SeedVx(m);  % Peak vertex in seed ROI
        fWimg(v) = fwhm; % Cortical full width half maximum guestimate
        fRimg(v) = mR^2;  % Peak correlation squared            
        % Fitting pRF parameters?
        if isempty(FitPrf)
            % Convex hull estimation
            txy = double(Temp(2:3,SeedVx))'; % Template X & Y coordinates
            [txy,u] = unique(txy, 'rows'); % Unique coordinates in seed map 
            R = R(u); % Remove duplicate correlation coefficients
            mR = max(R); % New maximum after redundant vertices removed

            % Region growing for positive subfield
            warning off
            tri = delaunay(txy); % Delaunay triangles
            pos = samsrf_clusterroi(find(R==mR,1), R, mR/2, tri); % Points in positive subfield
            
            % Quantify positive subfield
            ps = polyshape(txy(pos,:)); % Points in positive subfield
            ps = ps.convhull; % Convex hull of positive subfield
            cha = ps.area; % Area of convex hull
            [chx, chy] = ps.centroid; % Centroid of convex hull
            fXimg(v) = chx; % X-coordinate
            fYimg(v) = chy; % Y-coordinate
            fSimg(v) = sqrt(cha) / (2*sqrt(2*log(2))); % Positive FWHM converted to Sigma
            
            % Quantify negative subfield
            nR = min(R); % Minimum correlation
            neg = R < nR/2; % Points in negative subfield
            ps = polyshape(txy(neg,:)); % Points in negative subfield
            ps = ps.convhull; % Convex hull of negative subfield
            cha = ps.area; % Area of negative convex hull
            fIimg(v) = sqrt(cha) / (2*sqrt(2*log(2))); % Negative FWHM converted to Sigma
            fBimg(:,v) = [mR; nR];
            warning on
            
        elseif isinf(FitPrf)
            % Summary statistics only
            txy = Temp(2:3,SeedVx)'; % Template X & Y coordinates
            txy = unique(txy(R > mR/2,:), 'rows'); % Unique coordinates in CF profile
            fXimg(v) = nanmedian(txy(:,1)); % Median X-coordinate
            fYimg(v) = nanmedian(txy(:,2)); % Median Y-coordinate
            fSimg(v) = sqrt(mad(txy(:,1),1)^2 + mad(txy(:,2),1)^2); % Euclidean of median absolute deviations
            
        else
            % Explicit model fitting
            fP = [];
            fR = [];
            
            % Threshold correlations
            R = R - mR/2; % Threshold correlations by half maximum
            R(R<0) = 0; % Subthreshold pRF coordinates set to zero
            
            % Attempt multiple searches for Sigma
            for s = 10.^(-2:2)
                [cP,cR] = samsrf_fit2dprf(R, @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), Temp(2:3,SeedVx)'), [Temp(2:3,SeedVx(m))' s], [1 0 0 0], [], FitPrf); % Fit 2D model to pRF coordinates
                fP = [fP; cP]; % Parameters per iteration
                fR = [fR; cR]; % Model fit per iteration
            end
            ms = find(fR==nanmax(fR),1); % Search iteration with best correlation
            fP = fP(ms,:); % Best parameters
            fR = fR(ms); % Peak correlation
            
            % Store fit parameters
            fRimg(v) = fR; % Replace peak correlation with goodness of fit
            fXimg(v) = fP(1); % X-coordinate
            fYimg(v) = fP(2); % Y-coordinate
            fSimg(v) = fP(3); % Sigma parameter
            fBimg(:,v) = fP(4:5)'; % Beta parameters
        end
    end
    % Report back?
    if ProgReps
        send(Q,v); 
    end
end

% If convex hull fit combine sigmas
fSimg = [fSimg; fIimg];

    %% Nested progress report function
    function percentanalysed(~)
        samsrf_progbar(vc/nver);
        vc = vc + 1;
    end
end
