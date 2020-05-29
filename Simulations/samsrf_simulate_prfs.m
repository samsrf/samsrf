function Srf = samsrf_simulate_prfs(GtPars, PrfFcn, ApFrm, Model)
%
% Srf = samsrf_simulate_prfs(GtPars, PrfFcn, ApFrm, [Model])
%
% Simulates time courses for ground truth pRF parameters based on the 
%  pRF function PrfFcn and the stimulus apertures defined by ApFrm.
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
%
% Returns a Srf structure Srf which will contain the simulated timeseries
%  in Srf.Data. Srf.Ground_Truth contains the ground truth parameters. 
%  These timeseries are noise-free. You can add noise to them as desired.
%
%  If GtPars was a Srf structure the output will retain all the same fields
%  as the original (except for Srf.Data and Srf.Ground_Truth obviously) so 
%  you may wish to remove fields you don't care about. 
%
%  If the input was just a matrix with parameters then this output Srf is
%  not fully functional (contains no brain structure etc).
%
% 27/05/2020 - Written (DSS)
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
h = samsrf_waitbar('Simulating pRF timecourses...');
for v = 1:size(Srf.Data,2)
    Rfp = PrfFcn(Srf.Ground_Truth(:,v)', size(ApFrm,1)*2); % pRF profile
    Y = prf_predict_timecourse(Rfp, ApFrm, false); % Prediction without z-normalisation!
    Y = conv(Y, Model.Hrf); % Convolve with HRF
    Y = Y(1:size(ApFrm,3)); % Truncate tail end
    Srf.Data(:,v) = Y; % Store simulation
    samsrf_waitbar(v/size(Srf.Data,2), h);
end
samsrf_waitbar('', h);