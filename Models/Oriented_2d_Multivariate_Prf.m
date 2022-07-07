function Oriented_2d_Multivariate_Prf(DataPath, SrfFiles, Roi)
%
% Fits an oriented multivariate 2D pRF model
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%
% IMPORTANT NOTE: The search space for the coarse fit in this model has not
% yet been tested! You may also instead wish to use a seed map based on a
% standard 2D Gaussian instead of a coarse fit.
%

%% Mandatory parameters 
Model.Prf_Function = @(P,ApWidth) prf_multivariate_rf(P(1), P(2), P(3), P(4), P(5), ApWidth); % Which pRF model function? 
Model.Name = 'pRF_Multivariate'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma1'; 'Sigma2'; 'Phi'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1 1 0]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1 1 0]; % Which parameters must be positive? (refer to ModelHelp for issues with Nelder-Mead algorithm)
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pRF'; % Aperture file

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [.01 .01 .01 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space)
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in aperture space)

% Search grid for coarse fit
Model.Polar_Search_Space = true; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 10 : 350; % Polar search grid
Model.Param2 = 2 .^ (-5 : 0.2 : 0.6) * Model.Scaling_Factor; % Eccentricity  search grid
Model.Param3 = 2 .^ (-5:1) * Model.Scaling_Factor; % Horizontal sigma search grid
Model.Param4 = 2 .^ (-5:1) * Model.Scaling_Factor; % Vertical sigma search grid
Model.Param5 = 0 : 15 : 345; % Orientation search grid

%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Post-processing
load(MapFile); % Load map we just analysed
% Phi stays between -180 & +180
Srf.Data(6,:) = mod(Srf.Data(6,:), 360); % Ensure nothing above 360 or below 0
x = Srf.Data(6,:) > 180; % Index phi > 180 degrees
Srf.Data(6,x) = Srf.Data(6,x) - 360; % Greater than 180 is now negative
% Flip phi if sigma1 < sigma2
x = Srf.Data(4,:) < Srf.Data(5,:); % Is sigma1 < sigma2?
Srf.Data(6,x) = -Srf.Data(6,x); % Flip sign if sigma1 < sigma2
Srf.Data(6,:) = mod(Srf.Data(6,:), 180); % Now ensure it's within 0-180
% Add aspect ratio
Srf.Values{9} = 'Aspect Ratio'; 
Srf.Data(9,:) = abs(log(abs(Srf.Data(4,:)) ./ abs(Srf.Data(5,:)))); % Logarithm of Radial/Tangential (can't be negative)
% Save again
save(MapFile, 'Srf', 'Model', '-v7.3'); 

%% Return home
cd(HomePath); 