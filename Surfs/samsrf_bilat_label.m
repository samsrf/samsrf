function samsrf_bilat_label(Srf, Roi)
%
% samsrf_bilat_label(Srf, Roi)
%
% Saves a ROI label combining the left & right hemisphere. This requires a
%  Srf structure that combines the hemispheres as input.
%
%   Srf:  Surface data structure combining hemispheres
%   Roi:  ROI name without any prefix (e.g. 'V1').
%
% The function automatically searches for the labels prefixed as  'lh_' and
%  'rh_' to combine them into a single label file. They do not both need to 
%  exist, but if one doesn't the function returns a warning.
% The function returns an error if Srf.Nvert_Lhem doesn't exist, as this
%  indicates that the Srf is not for combined hemispheres.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if ~isfield(Srf, 'Nvert_Lhem')
    error('Not a combined hemispheres Srf!'); 
end

% Path to label
[p,f] = fileparts(Roi);
if isempty(p)
    p = '.';
end

% Load data
Lvx = samsrf_loadlabel([p filesep 'lh_' f]);
Rvx = samsrf_loadlabel([p filesep 'rh_' f]);

% Do labels exist?
if isnan(Lvx)
    warning(['lh_' Roi '.label does not exist']);
    Lvx = [];
end
if isnan(Rvx)
    warning(['rh_' Roi '.label does not exist']);
    Rvx = [];
end

% Adjust right hemisphere indices
Rvx = Rvx + Srf.Nvert_Lhem;

% Combined vertices
Vx = [Lvx; Rvx];

% Save new label
samsrf_srf2label(Srf, Roi, 1, Vx);
