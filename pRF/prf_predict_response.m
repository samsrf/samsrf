function [y overlaid stim_mask] = prf_predict_response(Stim, Rfp)
%
% y = prf_predict_response(Stim, Rfp)
%
% Returns the predicted response by overlaying the stimulus mask Stim
% with the receptive field profile Rfp and summing.
%
% 20/08/2018 - SamSrf 6 version (date added - no changes from v6) (DSS)
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
y = sum(overlaid(:));