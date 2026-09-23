function samsrf_hemi_label(Srf, Roi)
%
% samsrf_hemi_label(Srf, Roi)
%
% Saves two ROI labels splitting the bilateral label Roi into the left & right hemisphere. 
%  This requires a Srf structure of the combined hemispheres as input.
%
%   Srf:  Surface data structure combining hemispheres
%   Roi:  ROI name without any prefix (e.g. 'V1').
%
% The function saves the labels prefixed as  'lh_' and 'rh_'.  
% If the label is unilateral it the function returns a warning.
% The function returns an error if Srf.Nvert_Lhem doesn't exist, as this
%  indicates that the Srf is not for combined hemispheres.
%
% 02/09/2026 - Written (DSS)
%

if ~isfield(Srf, 'Nvert_Lhem')
    samsrf_error('Not a combined hemispheres Srf!'); 
end

% Path to label
[p,f] = fileparts(Roi);
if isempty(p)
    p = '.';
end

% Load data
Vx = samsrf_loadlabel([p filesep f]);
Lvx = Vx(Vx <= Srf.Nvert_Lhem); % Vertices in left hemisphere
Rvx = Vx(Vx > Srf.Nvert_Lhem); % Vertices in right hemisphere

% Do labels exist?
if isempty(Lvx)
    samsrf_disp(['WARNING: ' Roi '.label does not exist in left hemisphere']);
    Lvx = [];
end
if isempty(Rvx)
    samsrf_disp(['WARNING: ' Roi '.label does not exist in right hemisphere']);
    Rvx = [];
end

% Adjust right hemisphere indices
Rvx = Rvx - Srf.Nvert_Lhem; 

% Split hemisheres
[Lsrf, Rsrf] = samsrf_hemi_srfs(Srf);

% Save new labels
samsrf_srf2label(Lsrf, ['lh_' Roi], 1, Lvx);
samsrf_srf2label(Rsrf, ['rh_' Roi], 1, Rvx);
