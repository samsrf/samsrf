function Rfp = prf_gaussian_rf(x0, y0, sigma, screen_width)
%
% Rfp = prf_gaussian_rf(x0, y0, sigma, [screen_width])
%
% Returns a Gaussian receptive field profile with coordinates (x0,y0)
%   and a standard deviation of sigma. Units are in proportion of 
%   screen space with the origin at the centre of gaze.
%
%   The optional fourth input screen_width determines the size of the 
%   stimulus model screen in pixels. This should be twice the width
%   of the apertures.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 4
    screen_width = 200;
end
eccentricity = screen_width/4;

% Convert into pixels
x0 = x0 * eccentricity + screen_width/2;
y0 = y0 * eccentricity + screen_width/2;
sigma = sigma * eccentricity;

% Calculate the Gaussian
[X Y] = meshgrid(1:screen_width, fliplr(1:screen_width));
Rfp = exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma.^2));

