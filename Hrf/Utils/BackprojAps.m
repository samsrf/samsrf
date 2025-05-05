function [ApName, ApMax] = BackprojAps(SeedMap, SeedRoi, Data)
%
% [ApName, ApMax] = BackprojAps(SeedMap, SeedRoi, Data)
%
% This is an internal function used for pRF from CF analysis. 
% You should not usually have to worry about using this directly.
% An example Model function is included in SamSrf.
%
% Uses the seed region & map to project brain activity into visual space &
% then saves this as an aperture file named in the same folder as Data. 
% That file is called: aps_[Data]_[SeedROi].mat
%
% These apertures are then used in pRF analysis to estimate the CFs from 
% that seed region. Importantly, this uses the maximum eccentricity for 
% pixels, so Model.Scaling_Factor is also automatically set to this value.
%
% The function returns the name of the apertures in ApName & 
% the maximum of ApXY in ApMax.
%
% 07/03/2024 - Written (DSS)
% 13/09/2024 - Fixed help section & updated for SamSrf X usage (DSS)
% 17/09/2024 - Now ensures correct path for seed ROI label (DSS)
% 18/10/2024 - Fixed bug when only single functional run in Srf input (DSS)
%              Fixed bug with folder names (DSS)
%

% Load seed map
load(EnsurePath(SeedMap));
Srf = samsrf_expand_srf(Srf);
% Load seed region
roi = samsrf_loadlabel(EnsurePath(SeedRoi));
% pRF coordinates in seed region
ApXY = double(Srf.Data(2:3,roi)');
% Maximum value needed for scaling
ApMax = max(abs(ApXY(:))); 

% Data provided or load file?
if isstruct(Data)
    Srf = Data;
    if ischar(Srf.Functional)
        Srf.Functional = {Srf.Functional};
    end
    Data = Srf.Functional{1}; % Update name to first 
else
    load(EnsurePath(Data));
end
Srf = samsrf_expand_srf(Srf);
% Brain activity per pRF per volume
ApFrm = double(Srf.Data(:,roi)');

% Save apertures
[DataPath,Data] = fileparts(Data);
if isempty(DataPath)
    DataPath = '.';
end
[~,SeedRoi] = fileparts(SeedRoi);
ApName = [DataPath filesep 'aps_' Data '_' SeedRoi];
save(ApName, 'ApXY', 'ApFrm');