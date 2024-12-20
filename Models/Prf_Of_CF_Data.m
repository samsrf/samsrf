function Prf_Of_CF_Data(DataPath, SrfFiles, Roi, SeedRoi, TempMap)
% 
% Fits a standard 2D Gaussian pRF model but using apertures created by
% backprojecting activity from a seed region (see BackprojAps.m). 
% So effectively, this is a CF model.
%
%	DataPath:	Path where the mapping data are
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
%   SeedRoi:    Seed ROI label 
%   TempMap:    Template map file
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

%% Standard 2D Gaussian pRF
Model.Name = 'pRF-CF'; % File name to indicate type of pRF model
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); % Which pRF model function? 
Model.Param_Names = {'x0' 'y0' 'Sigma'}; % Names of parameters to be fitted
Model.Scaled_Param = [1 1 1]; % Which of these parameters are scaled 
Model.Only_Positive = [0 0 1]; % Which parameters must be positive?
Model.TR = 1; % Repetition time (TR) of pulse sequence
Model.Hrf = 1; % No HRF!
Model.SeedRoi = SeedRoi; % ROI label for seed region
Model.Template = TempMap; % Atlas retinotopic map of seed region
Model.Aperture_File = ''; % Aperture file is automatically set in this model!
Model.Scaling_Factor = NaN; % Scaling factor is automatically set in this model!

%% Search grid for coarse fit 
% Some parameters are multiplied with scaling factor as this must be in stimulus space!
Model.Polar_Search_Space = true; % (Optional) If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Param1 = 0 : 15 : 345; % Polar search grid 
Model.Param2 = 2 .^ (-5 : 0.1 : 0);  % Eccentricity search grid will be automatically scaled in this model!
Model.Param3 = 2 .^ (-5.6 : 0.1 : -0.5);  % Sigma search grid will be automatically scaled in this model!
            
%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit pRF model
OutFile = samsrf_fit_prf(Model, SrfFiles, Roi);
    
%% Return home
cd(HomePath); 