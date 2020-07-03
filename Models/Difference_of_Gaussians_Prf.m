function Difference_of_Gaussians_Prf(SrfFiles, Roi)
%
% Fits a Difference-of-Gaussians 2D pRF model
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%
% IMPORTANT NOTE: The search space for the coarse fit in this model has not
% yet been tested! You may also instead wish to use a seed map based on a
% standard 2D Gaussian instead of a coarse fit.
%

%% Difference of Gaussians pRF
Model.Prf_Function = @(P,ApWidth) prf_dog_rf(P(1), P(2), P(3), P(4), P(5), ApWidth); % Which pRF model function? 
Model.Name = 'pRF_DoG'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Centre'; 'Surround'; 'Delta'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1 1 0]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1 1 1]; % Which parameters must be positive?
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pRF'; % Aperture file
Model.Seed_Fine_Fit = ''; % (Optional) Define a Srf file to use as seed map
Model.Replace_Bad_Fits = false; % (Optional) If true, uses coarse fit for bad slow fits
Model.Smoothed_Coarse_Fit = 0; % (Optional) If > 0, smoothes data for coarse fit
Model.Coarse_Fit_Only = false; % (Optional) If true, only runs the coarse fit
Model.Fine_Fit_Threshold = 0.01; % (Optional) Define threshold for what to include in fine fit

% Search grid for coarse fit
Model.Polar_Search_Space = true; % (Optional) If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 15 : 345; % Polar search grid
Model.Param2 = 2 .^ (-5 : 0.4 : 0.6); % Eccentricity  search grid
Model.Param3 = 2 .^ (-5:1); % Centre sigma search grid
Model.Param4 = 2 .^ (-3:2); % Surround sigma search grid
Model.Param5 = 0 : 0.5 : 1; % Delta (surround/centre amplitude) search grid

%% Open dialogs if needed
HomePath = pwd;
% Choose data files
if nargin == 0
    [SrfFiles, PathName] = uigetfile('*h_*.mat', 'Choose SamSrf files', 'MultiSelect', 'on');
    if SrfFiles ~= 0
        cd(PathName);
    else
        error('No data files selected!');
    end
end
% Choose ROI label
if nargin <= 1
    [Roi, RoiPath] = uigetfile('*.label', 'Choose ROI label');
    if Roi ~= 0 
        Roi = [RoiPath Roi(1:end-6)];
    else
        Roi = '';
    end    
end

%% Fit pRF model
MapFile = samsrf_fit_prf(Model, SrfFiles, Roi);

%% Post-processing
load(MapFile);

% Calculate FWHM
Srf = samsrf_dog_fwhm(Srf);

% Field sign, CMF & smooth
R2_Threshold = 0.05; % R^2 threshold for surface calculations
Eccentricity_Range = [1 Model.Scaling_Factor]; % Eccentricity range for surface calculations
Smoothing_Kernels = [10 3]; % First kernel for field sign & second kernel for everything else
Srf = samsrf_surfcalcs(Srf, Roi, R2_Threshold, Eccentricity_Range, 'S', Smoothing_Kernels, false);

%% Save again
save(MapFile, 'Srf', 'Model', '-v7.3'); 

%% Return home
cd(HomePath); 