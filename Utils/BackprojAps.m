function BackprojAps(SeedMap, SeedRoi, Data)
%
% BackprojAps
%
% Uses the seed region & map to project brain activity into visual space &
% then saves this as an aperture file. You can then use this as input to
% pRF analysis to estimate the CFs from that seed region. Importantly, you
% must use the maximum eccentricity in the seed ROI as Model.Scaling_Factor!
% An example Model function is included in SamSrf.
%
% IMPORTANT: It probably goes without saying that the meshes for the seed map 
%            & the brain activity surface must be identical!
%
% 07/03/2024 - Written (DSS)
%

% Load seed map
load(SeedMap);
Srf = samsrf_expand_srf(Srf);
% Load seed region
roi = samsrf_loadlabel(SeedRoi);
% pRF coordinates in seed region
ApXY = double(Srf.Data(2:3,roi)');

% Load brain activity
load(Data);
Srf = samsrf_expand_srf(Srf);
% Brain activity per pRF per volume
ApFrm = double(Srf.Data(:,roi)');

% Save apertures
save(['aps_' Data], 'ApXY', 'ApFrm');