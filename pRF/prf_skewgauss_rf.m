function Rfp = prf_skewgauss_rf(x0, y0, sigma, alpha_x, alpha_y, screen_width)
%
% Rfp = prf_skewgauss_rf(x0, y0, sigma, alpha_x, alpha_y, [screen_width])
%
% Returns a skewed Gaussian receptive field profile with coordinates (x0,y0),
%   a standard deviation of sigma, and skew factors along of alpha_x/y along 
%   the x and y dimension, respectively. Units are in proportion of screen space 
%   with the origin at the centre of gaze.
%
%   The optional final input screen_width determines the size of the stimulus 
%   model screen in pixels. This should be twice the width of the apertures.
%
% 08/10/2021 - SamSrf 7 version (DSS)
%

if nargin < 6
    screen_width = 200;
end
eccentricity = screen_width/4;

% Convert into pixels
x0 = x0 * eccentricity + screen_width/2;
y0 = y0 * eccentricity + screen_width/2;
sigma = sigma * eccentricity;

% Calculate the skewed Gaussian
[X Y] = meshgrid(1:screen_width, fliplr(1:screen_width));
Rfp = 2 * exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma.^2)) .* normcdf((X-x0)*alpha_x) .* normcdf((Y-y0)*alpha_y);

