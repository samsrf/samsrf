function Vs = samsrf_borderpath(Srf, LabelVs)
%
% Vs = samsrf_borderpath(Srf, LabelVs)
%
% Creates a path surrounding the border of all the vertices in Srf that
% are in vector LabelVs. Returns a vector of vertex indeces in Vs.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

%% Calculate path
Vs = samsrf_neighbours(LabelVs, Srf.Faces); % Neighbours to labelled vertices
