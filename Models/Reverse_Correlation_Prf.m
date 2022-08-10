function Reverse_Correlation_Prf(DataPath, SrfFiles, Roi)
%
% Runs a reverse correlation pRF analysis
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis (see note below)
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%
% NOTE: This analysis uses a lot of memory. You will need a ROI as you will 
%       probably run out of memory otherwise. If you run a whole-brain analysis
%       it may be necessary to split that up into several ROIs.
%

%% Mandatory parameters 
Model.Name = 'pRF'; % File name to indicate type of pRF model
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pRF'; % Aperture file

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [.1 .1 .1]; % Use Hooke-Jeeves algorithm with these initial step sizes (in visual space) 
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in visual space)

%% For fitting 2D pRF models: Commment this whole section to use convex hull estimation!
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.SeedPar_Function = @(V) [V(2:3); V(4)/(2*sqrt(2*log(2)))]; % Seed parameter function
Model.R2_Threshold = 0.1; % Reverse correlation threshold

%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_revcor_prf(Model, SrfFiles, Roi);

%% Return home
cd(HomePath); 