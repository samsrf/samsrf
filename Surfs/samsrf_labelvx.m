function Vs = samsrf_labelvx(Srf, RoiVx)
%
% Vs = samsrf_labelvx(Srf, RoiVx)
%
% Returns a boolean vector labelling the vertices in Srf inside a ROI.
% Use this to transform a label list imported with samsrf_loadlabel into a mask.
%
% 05/11/2021 - Written (DSS)
% 13/01/2022 - Fixed mistake in help section (DSS)
%

Vs = false(1, size(Srf.Vertices,1));
Vs(RoiVx) = true;