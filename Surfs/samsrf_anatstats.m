function [Area, Thickness, Curvature] = samsrf_anatstats(Srf, Ecc, Roi)
%
% [Area, Thickness, Curvature] = samsrf_anatstats(Srf, [Ecc=NaN, Roi=''])
%
% Returns anatomical cortical statistics in an eccentricity range defined 
%  by the two-element vector Ecc. The measurement can be restricted to the 
%  label Roi. If Ecc is NaN the function calculates these statistics for 
%  the whole ROI or the whole hemisphere (this is the default).
%
% There is deliberately no R^2 threshold in this function. If restricting 
%  the analysis to an eccentricity range, you want to use very smooth maps 
%  with good cortical coverage for this. If you care about R^2 you should 
%  use samsrf_plot instead.
%
% Returns the surface area, mean cortical thickness, and mean curvature.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 2
    Ecc = NaN;
end
if nargin < 3
    Roi = '';
end

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Load region of interest
if ~isempty(Roi)
    Vs = samsrf_loadlabel(Roi);
else
    Vs = 1:size(Srf.Vertices,1);  
end

% Anatomical data
A = Srf.Area; % Surface area
T = Srf.Thickness; % Thickness
C = Srf.Curvature; % Curvature

% Calculate anatomical stats
if isnan(Ecc)
    Area = sum(A(Vs)); % Surface area
    Thickness = mean(T(Vs)); % Thickness
    Curvature = mean(C(Vs)); % Curvature
else
    A = A(Vs); % Surface area for whole ROI
    T = T(Vs); % Thickness for whole ROI
    C = C(Vs); % Curvature for whole ROI
    D = Srf.Data(:,Vs);  % pRF data for ROI
    E = sqrt(D(2,:).^2 + D(3,:).^2);  % Eccentricity for vertices
    VsE = E > Ecc(1) & E <= Ecc(2); % Eccentricities within range
    Area = nansum(A(VsE)); % Sum of area in eccentricity range
    Thickness = nanmean(T(VsE)); % Mean thickness in eccentricity range
    Curvature = nanmean(C(VsE)); % Mean curvature in eccentricity range
end