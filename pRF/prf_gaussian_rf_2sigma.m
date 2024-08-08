function Rfp = prf_gaussian_rf_2sigma(x0, y0, sigma1, sigma2, XY)
%
% Rfp = prf_gaussian_rf_2sigma(x0, y0, sigma1, sigma2, XY)
%
% Returns a Gaussian receptive field profile with coordinates (x0,y0) standard deviation of sigma1 & sigma2. 
% This is an example of a multi-condition model where pRF size can vary between 2 stimulus conditions. 
%
%	Units are in stimulus space relative to fixation.
%
%   The seventh input XY contains a n*2 matrix where each column contains X and Y coordinates of pRFs. 
%	These must be in stimulus space, not apertures space! This is used internally for fitting pRFs.
%	There is no 2D version for this pRF model but you can create that by reshaping the output vectors.
%
% 01/08/2024 - Created, based on standard 2D Gaussian model (DSS)
% 02/08/2024 - Changed to two-condition model (DSS)
% 04/08/2024 - Renamed to be more descriptive & changed help section (DSS)
%

% Prepare pRF profile
if size(XY,1) == 1
    error('This function only works with coordinate input!');
else
    % pRF coordinates provided as input
    X = XY(:,1);
    Y = XY(:,2);
	
	Rfp = zeros(size(XY,1),2);
end

% Calculate the Gaussian
sigma = [sigma1 sigma2]; % Current sigma
for s = 1:2
	Rfp(:,s) = exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma(s).^2)); % pRF profile for current sigma
end
