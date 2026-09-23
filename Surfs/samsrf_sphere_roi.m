function [Vs, Ds] = samsrf_sphere_roi(v, d, S)
%
% [Vs, Ds] = samsrf_sphere_roi(v, d, S)
%
% Finds the vertices in the neighbourhood of v 
% within Euclidean distance d in the sphere mesh S.
%
% 03/09/2026 - Written (DSS)
%

Ds = sqrt(sum((S - S(v,:)).^2,2)); % Euclidean distance of v to all sphere vertices
Vs = find(Ds <= d); % ROI vertices
Ds = Ds(Vs); % Distances of ROI vertices
