function [Vs, Rho, Ws] = samsrf_visualroi(xyr, X, Y, Sigma)
%
% [Vs, Rho, Ws] = samsrf_visualroi(xyr, X, Y, Sigma)
%
% Finds the vertices in the neighbourhood of xyr based on the pRF data 
% defined by X, Y and Sigma. xyr is a 3-element vector that defines the
% centre position and the radius of the ROI.
% 
% Returns in Vs the vertices inside the ROI, in Rho the distance, 
% and in Ws the weight of each vertex based on the distance and pRF size.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 05/07/2021 - Fixed output dimensions (DSS)
%

% Individual ROI parameters
x = xyr(1);
y = xyr(2);
r = xyr(3);

% Distance from ROI centre
Rho = sqrt((X-x).^2 + (Y-y).^2);

% Distance from ROI edge
Rho = Rho - r;

% All vertices within ROI have zero distance
Rho(Rho < 0) = 0;

% Vertices inside the ROI
Vs = find(Rho == 0);

% Calculate weights
Ws = exp(-(Rho.^2 ./ (2*Sigma.^2)));

% Ensure column outputs
if size(Vs,1) < size(Vs,2)
    Vs = Vs';
    Rho = Rho';
    Ws = Ws';
end