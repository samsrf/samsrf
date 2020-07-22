function Circular_Tuning_Curve(SrfFiles, Roi)
%
% Fits a 1D circular Gaussian tuning curve model with polar wedge apertures
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Circular tuning function 
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(cosd(P(1))/2, sind(P(1))/2, P(2), ApWidth); % Which pRF model function? 
Model.Name = 'pTC_cir'; % File name to indicate type of pRF model
Model.Param_Names = {'Phi'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [0 0]; % Which of these parameters are scaled 
Model.Only_Positive = [0 1]; % Which parameters must be positive?
Model.Scaling_Factor = 1; % Scaling factor of the stimulus space (stays fixed for this model!)
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pTC_cir'; % Aperture file

% Optional parameters
Model.Noise_Ceiling_Threshold = 0; % Limit data to above certain noise ceiling?
Model.Replace_Bad_Fits = false; % If true, uses coarse fit for bad slow fits
Model.Smoothed_Coarse_Fit = 0; % If > 0, smoothes data for coarse fit
Model.Coarse_Fit_Only = false; % If true, only runs the coarse fit
Model.Seed_Fine_Fit = ''; % Define a Srf file to use as seed map
Model.Fine_Fit_Threshold = 0.01; % Define threshold for what to include in fine fit
Model.Coarse_Fit_Block_Size = 10000; % Defines block size for coarse fit (reduce if using large search space)
Model.Polar_Search_Space = false; % If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates

% Search grid for coarse fit
Model.Param1 = -179 : +180; % Mu search grid
Model.Param2 = 2 .^ (-5.6 : 0.1 : 0); % Sigma search grid
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

%% Post-processing
load(MapFile); % Load map we just analysed
% Phi stays between -180 & +180
Srf.Data(2,:) = mod(Srf.Data(2,:), 360); % Ensure nothing above 360 or below 0
x = Srf.Data(2,:) > 180; % Index phi > 180 degrees
Srf.Data(2,x) = Srf.Data(2,x) - 360; % Greater than 180 is now negative
save(MapFile, 'Srf', 'Model', '-v7.3'); % Save again

%% Return home
cd(HomePath); 