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
Model.Alpha = 10^-6; % Significant p-value for correlation maps
Model.Rdim = 50; % Side length of correlation map

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