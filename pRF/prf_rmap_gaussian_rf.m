function Rfp = prf_rmap_gaussian_rf(X,Y,x0,y0,sigma,beta)
%
% Rfp = prf_rmap_gaussian_rf(X,Y,x0,y0,sigma,beta)
%
% Returns a Gaussian receptive field profile with coordinates (x0,y0)
%   a standard deviation of sigma, and an amplitude of beta. 
%   X and Y are mesh grids in the eccentricity range of the stimuli.
%
% This function will be important for fitting a pRF model to r-maps.
%   Currently this is not used by anything!!! 
%
% 20/08/2018 - SamSrf 6 version (data added - no changes from v6) (DSS)
%

% Calculate the Gaussian
Rfp = beta * exp(-((X-x0).^2+(Y-y0).^2) ./ (2*sigma.^2));
