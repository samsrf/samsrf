function Connective_Fields(DataPath, SrfFiles, Roi, SeedRoi, TempMap)
%
% Runs a connective field reverse correlation analysis
%	DataPath:	Path where the mapping data are
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

%% Mandatory parameters 
Model.Name = 'CF'; % File name for output map
Model.SeedRoi = SeedRoi; % Seed ROI for analysis
Model.Template = TempMap; % Vertex number of seed region origin

%% Optional in forward-model CF - You can use only one in this order of priority!
% Model.Polar = 0:15:315; % Use polar wedges as CFs
% Model.Eccentricity = [0 1 2 4 8 16 32 64 90]; % Use eccentricity bands as CFs
% Model.Sizes = 5:5:20; % Search space for CF sizes (in geodesic steps)

%% Go to data 
HomePath = pwd;
cd(DataPath);

%% Fit CF model (either reverse-correlation or forward-model)
MapFile = samsrf_revcor_cf(Model, SrfFiles, Roi); % Reverse correlation analysis
% MapFile = samsrf_fit_cf(Model, SrfFiles, Roi); % Forward-model fast fit (computationally expensive & not fully tested)

%% Return home
cd(HomePath); 