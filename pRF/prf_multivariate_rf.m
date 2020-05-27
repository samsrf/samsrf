function Rfp = prf_multivariate_rf(x0, y0, s1, s2, phi, screen_width)
% Rfp = prf_multivariate_rf(x0, y0, s1, s2, phi, [screen_width])
%
% Returns an multivariate oriented Gaussian receptive field profile 
% with coordinates (x0,y0), standard deviations s1 and s2 for the principal
% and orthogonal axes and orientation (phi) in degrees, along the unit circle. 
%
% The optional input screen_width determines the size of the stimulus model
% screen in pixels. This should be twice the width of the apertures.
%
% Note that if you use this pRF model you should add an unwrapping step to
% the postprocessing of your model data (see e.g. Oriented_2d_Multivariate_Prf
% or Circular_Tuning_Curve in SamSrf/Models).
%
% 20/08/2018 - SamSrf 6 version (date added - no changes from v5) (DSS)
% 27/05/2020 - Added unwrapping of angles
%

if nargin < 6
    screen_width = 200;
end
eccentricity = screen_width / 4;

% Convert into pixels
x0 = x0 * eccentricity + screen_width / 2;
y0 = y0 * eccentricity + screen_width / 2;
s1 = s1 * eccentricity;
s2 = s2 * eccentricity;

% Unwrap orientation
phi = mod(phi, 360);

% Orientation from degrees to radians
phi = phi / 180 * pi;

% As we'll flip the Y dimension (see below), flip the orientation index as well.
phi = (2 * pi) - phi;

% Search grid
[X, Y] = meshgrid(1:screen_width, fliplr(1:screen_width));

% Calculate the gaussian
a = (cos(phi)^2 / (2 * s1^2)) + (sin(phi)^2 / (2 * s2^2));
b = -(sin(2 * phi) / (2 * s1^2)) + (sin(2 * phi) / (2 * s2^2));
c = (sin(phi)^2 / (2 * s1^2)) + (cos(phi)^2 / (2 * s2^2));
Rfp = exp(-a*(X - x0).^2 - b*(X - x0).*(Y - y0) - c*(Y - y0).^2);
      
