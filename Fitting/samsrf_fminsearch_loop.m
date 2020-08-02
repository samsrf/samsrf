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
%

% Number of vertices
nver = size(Y,2);

% Output matrices
fPimg = zeros(length(Model.Param_Names), nver); % Parameter maps
fRimg = zeros(1,nver); % R^2 map

% Prepare progress report
try
    Q = parallel.pool.DataQueue;
    afterEach(Q, @percentanalysed);
    IsParallel = true;
    disp(' Parallel computing - progress update every 5000 vertices');
catch
    IsParallel = false;
    disp(' No parallel computing - progress update every 100 vertices');
end

% Display off & default tolerances (=slow)
OptimOpts = optimset('Display', 'off'); 
% Add output function to prevent loop from getting stuck
OptimOpts.OutputFcn = @samsrf_fminsearch_outfun;

%% Loop thru mask vertices 
vc = 1;
if IsParallel
    % Run parallel loop
    parfor v = 1:nver
        if Rimg(v) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
            % Find best prediction
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y(:,v)), Pimg(:,v)', OptimOpts);  
            fPimg(:,v) = fP;
            fRimg(1,v) = 1 - fR;
        end
        % Report back?
        send(Q,v); % Only reports every 1000 vertices
    end
else
    % Run normal loop
    for v = 1:nver
        if Rimg(v) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
            % Find best prediction
            [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y(:,v)), Pimg(:,v)', OptimOpts);  
            fPimg(:,v) = fP;
            fRimg(1,v) = 1 - fR;
        end
        % Reports back ?
        if mod(vc,100) == 0 % Every 100 vertices
            disp([' ' num2str(round(vc/nver*100)) '% completed']);
        end
        vc = vc + 1;
    end
end

    %% Nested progress report function
    function percentanalysed(~)
    % Reports back every 5000 vertices
        if mod(vc,5000) == 0
            disp([' ' num2str(round(vc/nver*100)) '% completed']);
        end
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
