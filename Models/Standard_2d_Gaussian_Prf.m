function Standard_2d_Gaussian_Prf(SrfFiles, Roi)
%
% Fits a standard 2D Gaussian pRF model
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Standard 2D Gaussian pRF
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Name = 'pRF_Gaussian'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1]; % Which parameters must be positive?
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
Model.Param1 = 0 : 10 : 350; % Polar search grid
Model.Param2 = 2 .^ (-5 : 0.2 : 0.6); % Eccentricity  search grid
Model.Param3 = 2 .^ (-5.6 : 0.2 : 1); % Sigma search grid
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

%% Post-processing
load(MapFile); % Load map we just analysed
Srf = samsrf_normr2(Srf); % Convert raw R^2 to normalised R^2
save(MapFile, 'Srf', 'Model', '-v7.3'); % Save again

%% Return home
cd(HomePath); 