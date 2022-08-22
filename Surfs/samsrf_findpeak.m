function PeakVx = samsrf_findpeak(Map, Roi)
%
% PeakVx = samsrf_clusterroi(Map, [Roi=''])
%
% Finds the peak statistic in Map (a row vector for each vertex in a Srf).
% This can be restricted to a ROI label. Returns the vertex index of the peak.
%
% 22/08/2022 - Thus 'tis writ (DSS)
%

% No ROI defined
if nargin < 2 || isempty(Roi)
    RoiVs = 1:size(Map,2); % All vertex indeces
else
    RoiVs = samsrf_loadlabel(Roi); % Load ROI vertex indeces
end

% Find peak statistic
Map = Map(RoiVs); % Restrict to ROI
v = find(Map == nanmax(Map),1); % Index of maximum in ROI
PeakVx = RoiVs(v); % Peak vertex in complete Srf map


