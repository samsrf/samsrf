function [fVimg, fXimg, fYimg, fWimg, fRimg, fSimg, fBimg] = samsrf_cfparam_loop(Area, ConFlds, SeedVx, Temp, FitPrf) 
%
% [fVimg, fXimg, fZimg, fWimg, fRimg, fSimg, fBimg] = samsrf_cfparam_loop(Area, ConFlds, SeedVx, Temp, FitPrf) 
%
% Internal function looping through solving the CF parameter fit for each mask vertex.
% Separated from main function for clarity & to allow progress bar.
%
%   Area:     Srf.Area from data structure 
%   ConFlds:  Matrix with CF correlation profiles, restricted to mask
%   SeedVx:   Seed ROI vertex indices
%   Temp:     Srf.Data from template structure
%   FitPrf:   Toggles wehther to fit a 2D Gaussian pRF to the CF correlation profile
%
% 05/11/2021 - Written (DSS)
% 15/02/2022 - Added option to fit to pRF coordinates (DSS)
% 07/04/2022 - Fitting pRF now thresholds correlations by half-maximum (DSS)
% 08/04/2022 - Now uses iterative search to home in on pRF size (DSS)
% 12/04/2022 - Removed inconsequential erroneous comment (DSS)
%

% Number of vertices
nver = size(ConFlds,2);

% Output matrices
fVimg = zeros(1,nver); % Seed vertex map
fXimg = zeros(1,nver); % X map
fYimg = zeros(1,nver); % Y map
fWimg = zeros(1,nver); % Width map
fRimg = zeros(1,nver); % R^2 map
if nargin < 5
    FitPrf = true;
end
if FitPrf
    % Fitting pRF parameters
    fSimg = zeros(1,nver); % pRF parameter estimates 
    fBimg = zeros(2,nver); % pRF parameter estimates 
else
    % No pRF fitting 
    fSimg = [];
    fBimg = [];
end

% Parallel processing?
try 
  gcp;
  IsParallel = true;
  disp(' Parallel computing!');
catch
  IsParallel = false;
  disp(' No parallel computing!');  
end

% Prepare progress report
try
    Q = parallel.pool.DataQueue; % Queue variable
    afterEach(Q, @percentanalysed); % Reporting function
    ProgReps = true; % Progress will be reported
catch
    Q = NaN; % No queue variable exists
    ProgReps = false; % No progress can be reported when using parallel computing
    if IsParallel
        disp(' No progress reports possible :(');
    end
end

%% Loop thru mask vertices 
vc = 1;
if IsParallel
    % Run parallel loop
    if ProgReps
        samsrf_progbar(0);
    end
    parfor v = 1:nver
        % Retrieve r-map
        R = ConFlds(:,v); % Correlation of time courses with regressors
        % Determine parameters
        mR = max(R); % Find peak activation in each map
        mR = mR(1); % Ensure only one value
        fwhm = sqrt(sum(Area(R > mR/2))); % Square root of area above half maximum
        m = find(R==mR,1); % Find peak coordinate
        % Threshold correlations
        R = R - mR/2; % Threshold correlations by half maximum
        R(R<0) = 0; % Subthreshold pRF coordinates set to zero
        % If good correlation
        if ~isnan(mR) && mR > 0
            fVimg(v) = SeedVx(m);  % Peak vertex in seed ROI
            fXimg(v) = Temp(2,SeedVx(m)); % Template X coordinate
            fYimg(v) = Temp(3,SeedVx(m)); % Template Y coordinate
            fWimg(v) = fwhm; % Full width half maximum guestimate
            fRimg(v) = mR^2;  % Peak correlation squared            
            % Fitting pRF parameters?
            if FitPrf
                % Attempt mulitple searches for Sigma
                fP = [];
                fR = [];
                for s = 10.^(-2:2)
                    [cP,cR] = samsrf_fit2dprf(R, @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), Temp(2:3,SeedVx)'), [fXimg(v) fYimg(v) s], [1 0 0 0]); % Fit 2D model to pRF coordinates
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
else
    % Run normal loop
    samsrf_progbar(0);
    for v = 1:nver
        % Retrieve r-map
        R = ConFlds(:,v); % Correlation of time courses with regressors
        % Determine parameters
        mR = max(R); % Find peak activation in each map
        mR = mR(1); % Ensure only one value
        fwhm = sqrt(sum(Area(R > mR/2))); % Square root of area above half maximum
        m = find(R==mR,1); % Find peak coordinate
        % Threshold correlations
        R = R - mR/2; % Threshold correlations by half maximum
        R(R<0) = 0; % Subthreshold pRF coordinates set to zero
        % If good correlation
        if ~isnan(mR) && mR > 0
            fVimg(v) = SeedVx(m);  % Peak vertex in seed ROI
            fXimg(v) = Temp(2,SeedVx(m)); % Left-Right coordinate
            fYimg(v) = Temp(3,SeedVx(m)); % Inferior-Superior coordinate
            fWimg(v) = fwhm; % Full width half maximum guestimate
            fRimg(v) = mR^2;  % Peak correlation squared 
            % Fitting pRF parameters?
            if FitPrf
                % Attempt mulitple searches for Sigma
                fP = [];
                fR = [];
                for s = 10.^(-2:2)
                    [cP,cR] = samsrf_fit2dprf(R, @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), Temp(2:3,SeedVx)'), [fXimg(v) fYimg(v) s], [1 0 0 0]); % Fit 2D model to pRF coordinates
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
        % Reports back 
        samsrf_progbar(vc/nver);
        vc = vc + 1;
    end
end

    %% Nested progress report function
    function percentanalysed(~)
        samsrf_progbar(vc/nver);
        vc = vc + 1;
    end
end
