function Connective_Fields(SrfFiles, Roi, SeedRoi, TempMap)
%
% Runs a connective field reverse correlation analysis
%   SrfFiles:   Cell array with SamSrf data files (without extension)
%   Roi:        ROI label to restrict analysis 
%   SeedRoi:    Seed ROI label 
%   TempMap:    Template map file
% Both inputs are optional. If undefined, a dialog is opened for user selection.
%
% Includes code examples for using either forward-model fast fit or reverse correlation analysis.
%
% This is an example model. Copy this file to your parent data folder and
% adapt the model parameters to suit your personal needs and desires.
%

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
% Choose seed ROI label
if nargin <= 2
    [SeedRoi, SeedRoiPath] = uigetfile('*.label', 'Choose seed ROI label');
    if SeedRoi ~= 0 
        SeedRoi = [SeedRoiPath SeedRoi(1:end-6)];
    else
        error('Seed region undefined!');
    end    
end
% Choose template map
if nargin <= 3
    [TempMap, TempMapPath] = uigetfile('*h_*.mat', 'Choose template map');
    if TempMap ~= 0 
        TempMap = [TempMapPath TempMap(1:end-4)];
    else
        error('Template map undefined!');
    end    
end

%% Mandatory parameters for CF analysis
Model.Name = 'CF'; % File name for output map
Model.SeedRoi = SeedRoi; % Seed ROI for analysis
Model.Template = TempMap; % Vertex number of seed region origin

%% Optional in forward-model CF - You can use only one in this order of priority!
% Model.Polar = 0:15:315; % Use polar wedges as CFs
% Model.Eccentricity = [0 1 2 4 8 16 32 64 90]; % Use eccentricity bands as CFs
% Model.Sizes = 5:5:20; % Search space for CF sizes (in geodesic steps)

%% Fit pRF model
MapFile = samsrf_revcor_cf(Model, SrfFiles, Roi); % Reverse correlation analysis
% MapFile = samsrf_fit_cf(Model, SrfFiles, Roi); % Forward-model fast fit (computationally expensive & not fully tested)

%% Return home
cd(HomePath); 