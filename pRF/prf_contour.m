function R = prf_contour(Srf, v)
%
% R = prf_contour(Srf, v)
%
% Returns the pRF profile for vertex v in Srf as a contour matrix. 
%
% Note: This function requires a Srf analysed with reverse correlation.
%
% 24/02/2020 - Ensured axes are square now (DSS)
% 07/04/2020 - Removed plotting to avoid redundancy with samsrf_showprf (DSS)
%

% Reverse correlation profile vector
Rmap = Srf.Rmaps(:,v);
% Square side width
Width = sqrt(size(Rmap,1));
% Reshape vector into a map
R = reshape(Rmap, Width, Width);
