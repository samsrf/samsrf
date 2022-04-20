function Srf = samsrf_revcorprf_loop(Srf, GoF, Model, ApFrm, FitParam)
%
% Srf = samsrf_revcorprf_loop(Srf, GoF, Model, ApFrm, FitParam)
%
% Internal function looping through solving the pRF parameter fit for reverse-correlation profiles.
% Separated from main function for clarity & to allow progress bar.
%
%   Srf:        Input surface data structure 
%   GoF:        Vector with vertex indices meeting the goodness-of-fit threshold
%   Model:      Model structure with parameters
%   ApFrm:      Aperture frames
%   FitParam:   Parameters for optimisation algorithm
%
% Returns the surface data structure with the parameter fits.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Parallel processing?
try 
  gcp;
  disp(' Parallel computing!');
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
    disp(' No progress reports possible :(');
end

%% Loop thru good vertices
vc = 1;
Data = NaN(length(Model.Scaled_Param)+3, length(GoF)); % 3 additional rows for R^2 & betas
% Run parallel loop
if ProgReps
    samsrf_progbar(0);
end
parfor v = 1:length(GoF)
    vx = GoF(v); % Current vertex
    [fP,fR] = samsrf_fit2dprf(prf_contour(Srf,vx), Model.Prf_Function, Model.SeedPar_Function(Srf.Raw_Data(:,vx)), [Model.Scaling_Factor Model.Scaled_Param], ApFrm, FitParam); % Fit 2D model
    % Is good fit?
    KeepFit = true;  
    for p = 1:length(Model.Param_Names)
        if Model.Scaled_Param(p)
            if abs(fP(p)) > Model.Scaling_Factor*1.5
                KeepFit = false;
            end
        end
    end
    % Only store good fits
    if KeepFit
        Data(:,v) = [fR fP]';
    end
    
    % Report back?
    if ProgReps
        send(Q,v); 
    end
end
% Store fits in Srf
Srf.Data(:,GoF) = Data;

    %% Nested progress report function
    function percentanalysed(~)
        samsrf_progbar(vc/length(GoF));
        vc = vc + 1;
    end
end



