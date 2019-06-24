function Radial_2d_Multivariate_Prf(SrfFiles, Roi)
%
% Fits a radially oriented multivariate 2D pRF model
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% This is an example model. Move this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Radial 2D multivariate pRF
Model.Prf_Function = @(P,ApWidth) prf_multivariate_rf(P(1), P(2), P(3), P(4), atan2(P(2),P(1))/pi*180, ApWidth); % Which pRF model function? 
Model.Name = 'pRF_Radial'; % File name to indicate type of pRF model
Model.Param_Names = {'x0'; 'y0'; 'Sigma1'; 'Sigma2'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1 1]; % Which parameters must be positive?
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_Bars'; % Aperture file
Model.Replace_Bad_Fits = false; % (Optional) If true, uses coarse fit for bad slow fits
Model.Smoothed_Coarse_Fit = 0; % (Optional) If > 0, smoothes data for coarse fit
Model.Coarse_Fit_Only = false; % (Optional) If true, only runs the coarse fit

% Search grid for coarse fit
Model.Param1 = -1 : 0.2 : 1; % X0 search grid
Model.Param2 = -1 : 0.2 : 1; % Y0 search grid
Model.Param3 = 2 .^ (-5:1); % Horizontal sigma search grid
Model.Param4 = 2 .^ (-5:1); % Vertical sigma search grid
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
load(MapFile);

% Add aspect ratio
Srf.Values{8} = 'Aspect Ratio'; 
Srf.Data(8,:) = log(abs(Srf.Data(4,:)) ./ abs(Srf.Data(5,:))); % Logarithm of Radial/Tangential 

% Field sign, CMF & smooth
R2_Threshold = 0.05; % R^2 threshold for surface calculations
Eccentricity_Range = [1 Model.Scaling_Factor]; % Eccentricity range for surface calculations
Smoothing_Kernels = [10 3]; % First kernel for field sign & second kernel for everything else
Srf = samsrf_surfcalcs(Srf, Roi, R2_Threshold, Eccentricity_Range, 'S', Smoothing_Kernels, false);

%% Save again
save(MapFile, 'Srf', 'Model', '-v7.3'); 

%% Return home
cd(HomePath); 