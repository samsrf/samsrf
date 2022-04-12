function Gaussian_Tuning_Curve_Vertical(SrfFiles, Roi)
%
% Fits a 1D Gaussian tuning curve model with bars moving vertically
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Standard 2D Gaussian pRF
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(0, P(1), P(2), ApWidth); % Which pRF model function? 
Model.Name = 'pTC_ver'; % File name to indicate type of pRF model
Model.Param_Names = {'Mu'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 1]; % Which parameters must be positive?
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pTC_ver'; % Aperture file

% Optional parameters
Model.Noise_Ceiling_Threshold = 0; % Limit data to above certain noise ceiling?
Model.Replace_Bad_Fits = false; % If true, uses coarse fit for bad slow fits
Model.Smoothed_Coarse_Fit = 0; % If > 0, smoothes data for coarse fit
Model.Coarse_Fit_Only = false; % If true, only runs the coarse fit
Model.Seed_Fine_Fit = ''; % Define a Srf file to use as seed map
Model.Fine_Fit_Threshold = 0.01; % Define threshold for what to include in fine fit
Model.Coarse_Fit_Block_Size = 10000; % Defines block size for coarse fit (reduce if using large search space)
% Model.Downsample_Predictions = 10; % Use for microtime resolution if stimulus timing is faster than TR
% Model.Hooke_Jeeves_Steps = [.01 .01]; % Use Hooke-Jeeves algorithm with these initial step sizes (in aperture space)
% Model.Nelder_Mead_Tolerance = 0.01; % When using Nelder-Mead algorithm, use this parameter tolerance (in aperture space)

% Search grid for coarse fit
Model.Polar_Search_Space = false; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = -1.05 : 0.05 : 1.05; % Mu search grid
Model.Param2 = 2 .^ (-5.6 : 0.1 : 1); % Sigma search grid
Model.Param3 = 0; % Unused
Model.Param4 = 0; % Unused
Model.Param5 = 0; % Unused
        
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

%% Return home
cd(HomePath); 