function Reverse_Correlation_Prf(SrfFiles, Roi)
%
% Runs a reverse correlation pRF analysis
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%

%% Parameters for reverse correlation analysis
Model.Name = 'pRF'; % File name to indicate type of pRF model
Model.Scaling_Factor = 10; % Scaling factor of the stimulus space (e.g. eccentricity)
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = []; % HRF file or vector to use (empty = canonical)
Model.Aperture_File = 'aps_pRF'; % Aperture file
Model.Rdim = 50; % Side length of correlation map
Model.Noise_Ceiling_Threshold = 0.2; % Limit data to above certain noise ceiling?
Model.Save_Rmaps = true; % Saving pRF correlation profiles in data file (optional)

%% Parameters when fitting 2D pRF models (comment this section if not desired)
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Param_Names = {'x0'; 'y0'; 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.R2_Threshold = 0.1; % Reverse correlation threshold
Model.SeedPar_Function = @(V) [V(2:3); V(4)/(2*sqrt(2*log(2)))]; % Seed parameter function

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
MapFile = samsrf_revcor_prf(Model, SrfFiles, Roi);

%% Return home
cd(HomePath); 