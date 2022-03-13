function Simulate_Noisefree_Prfs(ModAps, SimAps)
%
% Simulate_Noisefree_Prfs(ModAps, [SimAps])
%
% Fits a standard 2D Gaussian pRF model to simulated noise-free pRF data.
%
% ModAps defines the aperture file to model (omitting 'aps_' prefix).
%
% SimAps is optional & defines the apertures to use for -simulating- data.
%   Use this if you want to use different underlying data than you use for
%   fitting the pRF model (e.g. a scotoma model or different design).
%   Defaults to being equal to ModAps.
%
% This is an example script of the simulation analyses you can run.
%

if nargin < 2
    SimAps = ModAps; % Use same apertures for simulation & model fit
end

%% Standard 2D Gaussian pRF
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Name = ['pRF_Sim_' SimAps '_Mod_' ModAps]; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1]; % Which parameters must be positive?
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (e.g. eccentricity) - not changed here
Model.TR = 1; % Repetition time (TR) of pulse sequence - standard in our experiments
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = ['aps_' ModAps]; % Box standard sweeping bars design we typically use

%% Search grid for coarse fit
Model.Polar_Search_Space = true; % (Optional) If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 10 : 350; % Polar angle search grid
Model.Param2 = 2 .^ (-5 : 0.2 : 0.6); % Eccentricity  search grid
Model.Param3 = 2 .^ (-5.6 : 0.2 : 1); % Sigma search grid
Model.Param4 = 0; % Unused
Model.Param5 = 0; % Unused

%% Simulate pRFs 
[Theta, Rho, Sigma] = ndgrid(0:15:345, 2.^[-4.5 -3.5:.5:.5], .05*2.^(0:5)); % A simulated visual field with a range of pRF sizes
[X,Y] = pol2cart(Theta/180*pi, Rho); % Convert polar to Cartesian coordinates
Ground_Truth = [X(:) Y(:) Sigma(:)]'; % Matrix of ground truth parameters
load(EnsurePath(['aps_' SimAps])); % Load apertures to simulate responses for
Srf = samsrf_simulate_prfs(Ground_Truth, @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth), ApFrm, Model); % Simulate time courses
save(['sim_' SimAps], 'Srf'); % Save simulated data 
        
%% Fit pRF model
samsrf_fit_prf(Model, ['sim_' SimAps]);

