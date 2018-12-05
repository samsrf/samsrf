function V = samsrf_loadlabel(Roi)
%
% V = samsrf_loadlabel(Roi)
%
% Returns a vector with the vertex indices of the label Roi.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

V = Read_FreeSurfer([Roi '.label']);
V = V(:,1)+1;
