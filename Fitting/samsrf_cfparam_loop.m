function [fVimg, fXimg, fYimg, fWimg, fRimg] = samsrf_cfparam_loop(Area, ConFlds, SeedVx, Temp) 
%
% [fVimg, fXimg, fZimg, fWimg, fRimg] = samsrf_cfparam_loop(Area, ConFlds, SeedVx, Temp) 
%
% Internal function looping through solving the CF parameter fit for each mask vertex.
%
% Loops thru fitting parameters to reverse correlation CFs. Separated from
%  main function for clarity & to allow progress bar.
%
%   Area:     Srf.Area from data structure 
%   ConFlds:  Matrix with CF correlation profiles, restricted to mask
%   SeedVx:   Seed ROI vertex indices
%   Temp:     Srf.Data from template structure
%
% 05/11/2021 - Written (DSS)
%

% Number of vertices
nver = size(ConFlds,2);

% Output matrices
fVimg = zeros(1,nver); % Seed vertex map
fXimg = zeros(1,nver); % X map
fYimg = zeros(1,nver); % Y map
fWimg = zeros(1,nver); % Width map
fRimg = zeros(1,nver); % R^2 map

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
        if ~isnan(mR) && mR > 0
            fVimg(v) = SeedVx(m);  % Peak vertex in seed ROI
            fXimg(v) = Temp(2,SeedVx(m)); % Left-Right coordinate
            fYimg(v) = Temp(3,SeedVx(m)); % Inferior-Superior coordinate
            fWimg(v) = fwhm; % Full width half maximum guestimate
            fRimg(v) = mR^2;  % Peak correlation squared 
        end
        % Report back?
        if ProgReps
            send(Q,v); % Only reports every 1000 vertices
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
        if ~isnan(mR) && mR > 0
            fVimg(v) = SeedVx(m);  % Peak vertex in seed ROI
            fXimg(v) = Temp(2,SeedVx(m)); % Left-Right coordinate
            fYimg(v) = Temp(3,SeedVx(m)); % Inferior-Superior coordinate
            fWimg(v) = fwhm; % Full width half maximum guestimate
            fRimg(v) = mR^2;  % Peak correlation squared 
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


%% Output function to prevent loop from getting stuck  
function stop = samsrf_fminsearch_outfun(x, optimValues, state) 
    if isfinite(optimValues.fval)
        stop = false;
    else
        stop = true;
    end
end
