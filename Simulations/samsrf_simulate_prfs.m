function Srf = samsrf_simulate_prfs(GtPars, PrfFcn, ApFrm, Model)
%
% Srf = samsrf_simulate_prfs(GtPars, PrfFcn, ApFrm, [Model])
%
% Simulates time courses for ground truth pRF parameters based on the 
%  pRF function PrfFcn and the stimulus apertures defined by ApFrm.
%  Responses are modelled in terms of percentage of pRF activation. 
%
% Note: When generating noisy time courses it is advisable to z-normalise 
%  them but you must do this *after* adding noise.
%
% GtPars can be a matrix where each column is a "vertex" and each row 
%  is a parameter necessary for the pRF function (excluding betas).
%  In this case, parameters should be defined in aperture space, so that
%  the left edge of the mapping stimulus is x = -1 and the right at x = 1.
%
% Alternatively, it can be a Srf structure containing a ground truth map.
%  In this case you must also define Model so that the function knows your
%  eccentricity / scaling factor & which parameters need to be scaled.
%  Moreover, you should prefilter your data by setting any vertices in
%  Srf.Data(1,:) to 1 or 0, if they should be simulated or not, respectively.
%
% Model.Hrf can define the HRF to use just like in a normal model fit. 
%  If Model is undefined, this defaults to the canonical HRF and TR = 1s.
%  Model can also define downsampling of predictions if TR mismatches stimulus timing.
%
% Returns a Srf structure Srf which will contain the simulated timeseries
%  in Srf.Data. Srf.Ground_Truth contains the ground truth parameters. 
%  These timeseries are noise-free. You can add noise to them as desired.
%  You may then also z-score them but this must be after adding the noise.
%
%  If GtPars was a Srf structure the output will retain all the same fields
%  as the original (except for Srf.Data and Srf.Ground_Truth obviously) so 
%  you may wish to remove fields you don't care about. 
%
%  If the input was just a matrix with parameters then this output Srf is
%  not fully functional (contains no brain structure etc).
%
% 02/06/2020 - SamSrf 7 version (DSS) 
% 12/07/2021 - Added stand-by message since parallel progress reports are a pain (DSS)
% 22/09/2021 - Now also supports downsampling of simulated timeseries (DSS)
%

%% What type of input?
if isstruct(GtPars)
    % Input is a map file
    Srf = samsrf_expand_srf(GtPars);
    Srf.Ground_Truth = Srf.Data(2:end,:); % Store ground truth 
    Srf.Data = NaN(size(ApFrm,3), size(Srf.Ground_Truth,2)); % Time course data
    % Rescale parameters
    if nargin < 4
        error('Model must be defined!');
    end
    for p = 1:length(Model.Scaled_Param)
        if Model.Scaled_Param(p)
            Srf.Ground_Truth(p,:) = Srf.Ground_Truth / Model.Scaling_Factor;
        end
    end
else
    % Input is a parameter matrix
    Srf = struct;
    Srf.Version = samsrf_version;
    Srf.Hemisphere = 'sim'; % Dummy
    Srf.Vertices = NaN(size(GtPars,2),3); % Dummy
    Srf.Ground_Truth = GtPars; % Store ground truth
    Srf.Data = NaN(size(ApFrm,3), size(Srf.Ground_Truth,2)); % Time course data
end

%% Load or generate HRF
if nargin < 4
    % Use canonical HRF if undefined
    Model = struct;
    Model.TR = 1;
    Model.Hrf = [];
    Model.Downsample_Predictions = 1; 
end
% Downsampling timeseries?
if ~isfield(Model, 'Downsample_Predictions')
    Model.Downsample_Predictions = 1; % Downsampling factor by which Model.TR mismatches the true TR
end

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

%% Simulate pRF timecourses
disp('Simulating data...');
disp(' Please stand by...');
Gt = Srf.Ground_Truth; % Ground truth data
Data = zeros(size(Srf.Data)); % Simulated data
parfor v = 1:size(Srf.Data,2)
    Rfp = PrfFcn(Gt(:,v)', size(ApFrm,1)*2); % pRF profile
    Y = prf_predict_timecourse(Rfp, ApFrm); % Prediction without z-normalisation!
    Y = prf_convolve_hrf(Y, Model.Hrf, Model.Downsample_Predictions); % Convolve with HRF & downsample if desired
    Data(:,v) = Y; % Store simulation
end
Srf.Data = Data; % Save simulation in Srf