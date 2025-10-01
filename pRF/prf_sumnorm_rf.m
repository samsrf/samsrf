function Rfp = prf_sumnorm_rf(x0, y0, sigma, kappa, screen_width)
%
% Rfp = prf_sumnorm_rf(x0, y0, sigma, kappa, [screen_width])
%
% Returns a Gaussian receptive field profile plus the same Gaussian
% multiplied with its logarithm. Any values > 1 are clipped to 1.
% This allows fitting a complex shape varying from centre-surround
% inhibition to nonlinear summation by only fitting one extra parameter.
% 
% 	The pRF is described by its position (x0,y0), the standard deviation 
% 	of sigma, and the modulating factor kappa. Units are in proportion of 
%   screen space with the origin at the centre of gaze.
%
%   The optional fifth input screen_width determines the size of the 
%   stimulus model screen in pixels. This should be twice the width
%   of the apertures.
%
%   Alternatively, screen_width can contain a n*2 matrix where each column
%   contains X and Y coordinates of pRFs. These must be in stimulus space, 
%   not apertures space! This is used internally for fitting pRFs & CFs 
%   but you may find other uses for this feature.
%
% 30/09/2025 - Written (DSS)
% 01/10/2025 - Fixed issue with NaN outputs (DSS) 
%              Renamed modulating parameter to kappa (DSS)
%

if nargin < 5
    screen_width = 200;
end

% Prepare pRF profile
if size(screen_width,1) == 1
    eccentricity = screen_width/4;

    % Convert into pixels
    x0 = x0 * eccentricity + screen_width/2;
    y0 = y0 * eccentricity + screen_width/2;
    sigma = sigma * eccentricity;

    [X Y] = meshgrid(1:screen_width, fliplr(1:screen_width));
else
    % pRF coordinates provided as input
    X = screen_width(:,1);
    Y = screen_width(:,2);
end

% Calculate the Gaussian
Rfp = exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma.^2));
% Modulate by itself
Rfp = Rfp + kappa * log(Rfp) .* Rfp;
% Clip to 1
Rfp(Rfp > 1) = 1;
% Remove bad points
Rfp(isnan(Rfp)) = 0;