function V = samsrf_loadlabel(Roi)
%
% V = samsrf_loadlabel(Roi)
%
% Returns a vector with the vertex indices of the label Roi.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if exist('OCTAVE_VERSION', 'builtin') <= 4
    % Matlab
    V = Read_FreeSurfer([Roi '.label']);
else
    % Octave
    V = textread([Roi '.label']); % Reads in a column vector
    V = V(~isnan(V)); % Remove NaNs (should only be at start)
    V = V(2:5:end); % Only extract the vertex indeces
end
V = V(:,1)+1;
