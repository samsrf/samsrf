function Rmaps = samsrf_cfcorrel_loop(Tc, X, Ws) 
%
% Rmaps = samsrf_cfcorrel_loop(Tc, X, Ws) 
%
% Internal function looping through calculating CF correlation profiles for each mask vertex.
% Separated from main function for clarity & to allow progress bar.
%
%   Tc: Time course data (limited to mask vertices)
%   X:  Time courses for seed region 
%   Ws: Weight matrix for smoothing (leave empty if no smoothing)
%
% 22/02/2022 - Written (DSS)
%

% Dimensions
nsvx = size(X,2); % Number of seed vertices
nmvx = size(Tc,2); % Number of mask vertices

% Output matrix
Rmaps = NaN(nsvx, nmvx); % Connective field profile for each vertex

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
    parfor v = 1:nmvx
        % Calculate r-map
        Y = Tc(:,v);  % Time course of current vertex
        R = corr(Y,X); % Correlation of time courses with regressors
        % Smooth profile if needed
        if ~isempty(Ws) 
            sR = repmat(R,nsvx,1) .* Ws; % Smoothed profiles for each seed vertex 
            R = nansum(sR,2) ./ sum(Ws,2); % Smoothed seed ROI    
        end
        % Determine parameters
        mR = max(R); % Find peak activation in each map
        mR = mR(1); % Ensure only one value
        if ~isnan(mR) && mR > 0
            Rmaps(:,v) = R(:); % Activation map as vector 
        end
        % Report back?
        if ProgReps
            send(Q,v); % Only reports every 1000 vertices
        end
    end
else
    % Run normal loop
    samsrf_progbar(0);
    for v = 1:nmvx
        % Calculate r-map
        Y = Tc(:,v);  % Time course of current vertex
        R = corr(Y,X); % Correlation of time courses with regressors
        % Smooth profile if needed
        if ~isempty(Ws) 
            sR = repmat(R,nsvx,1) .* Ws; % Smoothed profiles for each seed vertex 
            R = nansum(sR,2) ./ sum(Ws,2); % Smoothed seed ROI    
        end
        % Determine parameters
        mR = max(R); % Find peak activation in each map
        mR = mR(1); % Ensure only one value
        if ~isnan(mR) && mR > 0
            Rmaps(:,v) = R(:); % Activation map as vector 
        end
        % Reports back 
        samsrf_progbar(vc/nmvx);
        vc = vc + 1;
    end
end

    %% Nested progress report function
    function percentanalysed(~)
        samsrf_progbar(vc/nmvx);
        vc = vc + 1;
    end
end
