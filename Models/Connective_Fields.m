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

%% Parameters for analysis
Model.Name = 'CF'; % File name for output map
Model.SeedRoi = SeedRoi; % Seed ROI for analysis
Model.Template = TempMap; % Vertex number of seed region origin
Model.Smoothing = 0; % Smoothing kernel (works differently for forward-model & reverse-correlation)
Model.Sizes = 5:5:20; % Range of sizes for CFs in geodesic steps (forward-model only)

%% Fit pRF model
MapFile = samsrf_fit_cf(Model, SrfFiles, Roi); % Forward-model fast fit
% MapFile = samsrf_revcor_cf(Model, SrfFiles, Roi); % Reverse correlation analysis

%% Return home
cd(HomePath); 