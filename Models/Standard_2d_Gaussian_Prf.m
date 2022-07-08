function Standard_2d_Gaussian_Prf(DataPath, SrfFiles, Roi)
%
% Fits a standard 2D Gaussian pRF model
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Mandatory parameters 
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Name = 'pRF_Gaussian'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1]; % Which parameters must be positive? (refer to ModelHelp for issues with Nelder-Mead algorithm)
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pRF'; % Aperture file

%% Search grid for coarse fit 
% Some parameters are multiplied with scaling factor as this must be in stimulus space!
Model.Polar_Search_Space = true; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 10 : 350; % Polar search grid
Model.Param2 = 2 .^ (-5 : 0.2 : 0.6) * Model.Scaling_Factor; % Eccentricity  search grid
Model.Param3 = 2 .^ (-5.6 : 0.2 : 1) * Model.Scaling_Factor; % Sigma search grid

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [.01 .01 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space)
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in aperture space)
        
%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Return home
cd(HomePath); 