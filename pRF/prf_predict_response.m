function [y, overlaid, stim_mask] = prf_predict_response(Stim, Rfp)
%
% [y, overlaid, stim_mask] = prf_predict_response(Stim, Rfp)
%
% Returns the predicted response by overlaying the stimulus mask Stim
% with the receptive field profile Rfp and averaging across the frame.
% This is expressed as a percentage (in previous versions this was 
% calculated as the pixel sum instead - this shouldn't matter but it is
% still somewhat awkward whereas a percentage of the active frame makes
% more intuitive sense).
%
% Also returns the stimulus aperture and the overlay of pRF and aperture.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    DispTc = false;
end

% Constant paramaters
eccentricity = size(Stim,2)/2;
screen_width = size(Rfp,2);
midpoint = screen_width/2;

% Centre stimulus mask in screen
if size(Rfp,1) ~= size(Rfp,2)
    stim_mask = Stim;
else
    stim_mask = zeros(screen_width, screen_width);
    stim_mask(midpoint-eccentricity+1:midpoint+eccentricity, midpoint-eccentricity+1:midpoint+eccentricity) = Stim;
end

% Overlay the two images
overlaid = Rfp .* stim_mask;

% Predicted response
y = mean(overlaid(:)) / mean(Rfp(:)) * 100; % 100% = whole pRF is stimulated
