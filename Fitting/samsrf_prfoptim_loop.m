function [fPimg, fRimg] = samsrf_prfoptim_loop(Model, Y, ApFrm, ApXY, Rimg, Pimg) 
%
% [fPimg, fRimg] = samsrf_prfoptim_loop(Model, Y, ApFrm, ApXY, Rimg, Pimg) 
%
% Internal function looping through solving the fine fit for each mask vertex.
%
% Loops thru solving fine fit. 'Tis necessary to ensure the parallel computing 
% toolbox accepts it. Go to hell, Matlab!
%
%   Model:  Model structure for the pRF fit
%   Y:      Matrix with time series for every vertex, restricted to mask
%   ApFrm:  Matrix with vectorised stimulus apertures
%   ApXY:   Pixel coordinates for vectorised apertures
%   Rimg:   Coarse fit/seed map R^2 data for thresholding, restricted to mask
%   Pimg:   Coarse fit/seed map parameters, restricted to mask
%
% 06/07/2022 - Rewritten for vectorised apertures (DSS)
%

% Number of vertices
nver = size(Y,2);

% Output matrices
fPimg = zeros(length(Model.Param_Names), nver); % Parameter maps
fRimg = zeros(1,nver); % R^2 map

% Which algorithm?
if isfield(Model, 'Hooke_Jeeves_Steps')
    UseHookeJeeves = true; % Use HookeJeeves pattern search
else
    UseHookeJeeves = false; % Use Nelder-Mead (fminsearch)
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

% Display off & default tolerances (=slow)
OptimOpts = optimset('Display', 'off'); 
% Add output function to prevent loop from getting stuck
OptimOpts.OutputFcn = @samsrf_fminsearch_outfun;
% Was parameter estimation tolerance defined?
if isfield(Model, 'Nelder_Mead_Tolerance')
    OptimOpts.TolX = Model.Nelder_Mead_Tolerance;
end

%% Loop thru mask vertices 
vc = 1;
% Run parallel loop
if ProgReps
    samsrf_progbar(0);
end
parfor v = 1:nver
    if Rimg(v) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
        % Find best prediction
        if UseHookeJeeves
            % Use Hooke-Jeeves
            [fP,fR] = samsrf_hookejeeves(@(P) prf_errfun(Model.Prf_Function, ApFrm, ApXY, Model.Hrf, P, Y(:,v), Model.Downsample_Predictions, Model.Compressive_Nonlinearity), ...
                        Pimg(:,v)', Model.Hooke_Jeeves_Steps, Model.Only_Positive, 15, 3);  
        else
            % Use Nelder-Mead
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, ApXY, Model.Hrf, P, Y(:,v), Model.Downsample_Predictions, Model.Compressive_Nonlinearity), Pimg(:,v)', OptimOpts);  
        end
        fPimg(:,v) = fP;
        fRimg(1,v) = 1 - fR;
    end
    % Report back?
    if ProgReps
        send(Q,v); 
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
