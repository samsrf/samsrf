function [Ptc, S] = cf_generate_searchspace(Srf, SeedVx, Sizes)
%
% [Ptc, S] = cf_generate_searchspace(Srf, SeedVx, Sizes)
%
% Search grid Ptc for the coarse fit of a Gaussian CF model.
% Each column of Ptc contains a predicted time course for a point in the search grid.
% Each column in S contains the parameters for the predictions in Ptc.
%
% Srf is surface data file with the response time series in Srf.Data.
% SeedVx defines the vertices of the seed ROI label.
% Sizes defines the different CF sizes to fit for each ROI vertex (in geodesic steps).
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Response time series
Y = Srf.Data;

% Prediction matrix
Ptc = NaN(size(Y,1), length(SeedVx)*length(Sizes));

% Search grid 
[S1, S2] = ndgrid(SeedVx, Sizes);

% Parameter matrix
S = [S1(:) S2(:)]'; 

% Generating predictions
samsrf_disp(' Please stand by...');
parfor n = 1:numel(S1)
    Cir = samsrf_georoi(S1(n), S2(n), Srf.Vertices, Srf.Faces); % Circular CF patch
    svx = ismember(Cir,SeedVx); % Vertices in patch overlapping seed ROI
    Cir = Cir(svx); % Patch inside seed ROI
    Yc = Y(:,Cir); % Time series in patch
    Ptc(:,n) = roi_1steigvar(Yc); % 1st eigenvariate of patch
end
