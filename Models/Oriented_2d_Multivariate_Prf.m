function Oriented_2d_Multivariate_Prf(DataPath, SrfFiles, Roi)
%
% Fits an oriented multivariate 2D pRF model
%
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%
% IMPORTANT NOTE: The search space for the coarse fit in this model has not
% yet been tested! You may also instead wish to use a seed map based on a
% standard 2D Gaussian instead of a coarse fit.
%

%% Mandatory parameters 
Model.Name = 'pRF_Multivariate'; % File name to indicate type of pRF model
Model.Prf_Function = @(P,ApWidth) prf_multivariate_rf(P(1), P(2), P(3), P(4), atan2(P(2),P(1))/pi*180 + P(5), ApWidth); % Which pRF model function? 
Model.Param_Names = {'x0' 'y0' 'Sigma1' 'Sigma2' 'Theta'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1 1 0]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1 1 0]; % Which parameters must be positive? (refer to ModelHelp for issues with Nelder-Mead algorithm)
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = 0; % HRF file or vector to use (0 = SPM canonical, [] = de Haas canonical)
Model.Aperture_File = ''; % Aperture file must be defined!

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [.01 .01 .01 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space!)
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in visual space!)

%% Search grid for coarse fit
% Some parameters are multiplied with scaling factor as this must be in stimulus space!
Model.Polar_Search_Space = true; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 10 : 350; % Polar search grid
Model.Param2 = 2 .^ (-5 : 0.2 : 0.6) * Model.Scaling_Factor; % Eccentricity  search grid
Model.Param3 = 2 .^ (-5:1) * Model.Scaling_Factor; % Horizontal sigma search grid
Model.Param4 = 2 .^ (-5:1) * Model.Scaling_Factor; % Vertical sigma search grid
Model.Param5 = -90 : 15 : 90; % Orientation offset (relative to radial axis) search grid

%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Post-processing
load(MapFile); % Load map we just analysed
Srf.Data(6,:) = mod(Srf.Data(6,:), 180); % Ensure that 0 <= Theta < 180 
Srf.Data(6,Srf.Data(:,6) > 90) = Srf.Data(6,Srf.Data(:,6) > 90) - 180; % Ensure -90 < Theta <= 90
% Add aspect ratio
Srf.Values{end+1} = 'Aspect Ratio'; 
Srf.Data = [Srf.Data; log2(Srf.Data(4,:) ./ Srf.Data(5,:))]; % Logarithm of Radial/Tangential 
% Save again
save(MapFile, 'Srf', 'Model', '-v7.3'); 

%% Return home
cd(HomePath); 