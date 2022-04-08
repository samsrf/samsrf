function V = samsrf_loadlabel(Roi)
%
% V = samsrf_loadlabel(Roi)
%
% Returns a vector with the vertex indices of the label Roi.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 09/04/2022 - Removed Octave 4 support since this didn't work anyway (DSS)
%

V = Read_FreeSurfer([Roi '.label']);
V = V(:,1)+1;
