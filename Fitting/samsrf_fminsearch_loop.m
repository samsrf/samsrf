function [fPimg, fRimg] = samsrf_fminsearch_loop(Model, Y, ApFrm, Rimg, Pimg) 
%
% [fPimg, fRimg] = samsrf_fminsearch_loop(Model, Y, ApFrm, Rimg, Pimg) 
%
% Internal function looping through solving the fine fit for each mask vertex.
%
% Loops thru solving fine fit. 'Tis necessary to ensure the parallel computing 
% toolbox accepts it. Go to hell, Matlab!
%
%   Model:  Model structure for the pRF fit
%   Y:      Matrix with time series for every vertex, restricted to mask
%   ApFrm:  Matrix with stimulus apertures
%   Rimg:   Coarse fit/seed map R^2 data for thresholding, restricted to mask
%   Pimg:   Coarse fit/seed map parameters, restricted to mask
%
% 20/07/2020 - SamSrf 7 version (DSS)
% 03/08/2020 - Fixed bug where loop could get stuck on bad fits (IA & DSS)
% 21/12/2020 - No progress reports if parallel computing but using older Matlab versions (DSS)
% 22/12/2020 - Bugfix for when progress reports are turned off (DSS)
% 30/06/2021 - Added new-fangled old-school command-line progress-bars (DSS)
% 22/09/2021 - Now supports downsampling of predictions if TR mismatches stimulus timing (DSS)
%

% Number of vertices
nver = size(Y,2);

% Output matrices
fPimg = zeros(length(Model.Param_Names), nver); % Parameter maps
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

% Display off & default tolerances (=slow)
OptimOpts = optimset('Display', 'off'); 
% Add output function to prevent loop from getting stuck
OptimOpts.OutputFcn = @samsrf_fminsearch_outfun;

%% Loop thru mask vertices 
vc = 1;
if IsParallel
    % Run parallel loop
    if ProgReps
        samsrf_progbar(0);
    end
    parfor v = 1:nver
        if Rimg(v) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
            % Find best prediction
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y(:,v), Model.Downsample_Predictions), Pimg(:,v)', OptimOpts);  
            fPimg(:,v) = fP;
            fRimg(1,v) = 1 - fR;
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
        if Rimg(v) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
            % Find best prediction
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y(:,v), Model.Downsample_Predictions), Pimg(:,v)', OptimOpts);  
            fPimg(:,v) = fP;
            fRimg(1,v) = 1 - fR;
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
