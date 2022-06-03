function Gaussian_Tuning_Curve_Vertical(DataPath, SrfFiles, Roi)
%
% Fits a 1D Gaussian tuning curve model with (horizontal) bars moving vertically
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Mandatory parameters 
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(0, P(1), P(2), ApWidth); % Which pRF model function? 
Model.Name = 'pTC_ver'; % File name to indicate type of pRF model
Model.Param_Names = {'Mu'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 1]; % Which parameters must be positive? (refer to ModelHelp for issues with Nelder-Mead algorithm)
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pTC_ver'; % Aperture file

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [.01 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space)
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in aperture space)

% Search grid for coarse fit
Model.Polar_Search_Space = false; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = -1.05 : 0.05 : 1.05; % Mu search grid
Model.Param2 = 2 .^ (-5.6 : 0.1 : 1); % Sigma search grid
        
%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Return home
cd(HomePath); 