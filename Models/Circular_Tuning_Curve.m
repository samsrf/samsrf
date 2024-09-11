function Circular_Tuning_Curve(DataPath, SrfFiles, Roi)
%
% Fits a 1D circular Gaussian tuning curve model with polar wedge apertures
%
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Mandatory parameters 
Model.Name = 'pTC_cir'; % File name to indicate type of pRF model
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(cosd(P(1))/2, sind(P(1))/2, P(2), ApWidth); % Which pRF model function? 
Model.Param_Names = {'Mu' 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [0 0]; % Which of these parameters are scaled 
Model.Only_Positive = [0 1]; % Which parameters must be positive? (refer to ModelHelp for issues with Nelder-Mead algorithm)
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (stays fixed for this model!)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = ''; % Aperture file must be defined!

%% Optional fine-fitting parameters
% Model.Hooke_Jeeves_Steps = [5 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space!)
% Model.Nelder_Mead_Tolerance = 0.01; % Define parameter tolerance for Nelder-Mead algorithm (in visual space!)

%% Search grid for coarse fit
% Some parameters are multiplied with scaling factor as this must be in stimulus space!
Model.Polar_Search_Space = false; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = -179 : +180; % Mu search grid
Model.Param2 = 2 .^ (-5.6 : 0.1 : 0) * Model.Scaling_Factor; % Sigma search grid
        
%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Post-processing
load(MapFile); % Load map we just analysed
% Mu stays between -180 & +180
Srf.Data(2,:) = mod(Srf.Data(2,:), 360); % Ensure nothing above 360 or below 0
x = Srf.Data(2,:) > 180; % Index mu > 180 degrees
Srf.Data(2,x) = Srf.Data(2,x) - 360; % Greater than 180 is now negative
save(MapFile, 'Srf', 'Model', '-v7.3'); % Save again

%% Return home
cd(HomePath); 