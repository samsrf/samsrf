function Rfp = prf_dog_rf(x0, y0, sigma1, sigma2, delta, screen_width)
%
% Rfp = prf_dog_rf(x0, y0, sigma1, sigma2, delta, [screen_width])
%
% Returns a Difference-of-Gaussian receptive field profile with 
%   coordinates (x0,y0) and standard deviations sigma1 and sigma2.
%   The ratio between the two amplitudes is given by delta.
%   Units are in proportion of screen space with the origin 
%   at the centre of gaze.
%
%   The optional sixth input screen_width determines the size of the 
%   stimulus model screen in pixels. This should be twice the width
%   of the apertures.
%
% 20/08/2018 - SamSrf 6 version (data added - no changes from v6) (DSS)
%

if nargin < 6
    screen_width = 200;
end
eccentricity = screen_width/4;

% Convert into pixels
x0 = x0 * eccentricity + screen_width/2;
y0 = y0 * eccentricity + screen_width/2;
sigma1 = sigma1 * eccentricity;
sigma2 = sigma2 * eccentricity;

% Calculate the Gaussian
[X Y] = meshgrid(1:screen_width, fliplr(1:screen_width));
Rfp = exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma1.^2)) - delta*exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma2.^2));

