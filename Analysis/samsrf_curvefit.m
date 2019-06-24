function [mP,mR,oF] = samsrf_curvefit(xyci, ft)
%
% [mP,mR,oF] = samsrf_curvefit(xyci, ft)
%
% Fits a curve to pRF data from an individual subject (Sigma, CMF, etc).
% This requires Matlab's Curve Fitting toolbox.
%
%   xyci:   Matrix output from samsrf_plot (e.g. eccentricity, value, CIs, N). 
%           If CIs are included, points are weighted by their reciprocal.
%
%   ft:     What function to fit. Refer to the List of Library Models for 
%           Curve and Surface Fitting in Matlab Help. In addition you can 
%           use 'cumulgauss' which applies: a*(1+erf((x-b)/(sqrt(2)*c))).  
%           Or you also define your own function as a string, e.g.:
%               'a * sin(x*b+c) + d'
%
% Returns the fitted parameters mP, the goodness-of-fit mR structure, 
% and the fitted function oF (use as oF(x)). 
%
% 03/08/2018 - SamSrf 6 version (DSS)
% 

% Split variables
x = xyci(:,1); % Eccentricity
y = xyci(:,2); % Dependent variable
if size(xyci,2) > 2
    % CIs are defined so weight it
    w = 1./(xyci(:,4)-xyci(:,3)); 
else
    % No CIs are defined so no weighting
    w = ones(size(x));
end

% Remove bad data
x(isnan(y)) = [];
w(isnan(y)) = [];
y(isnan(y)) = [];
x(isnan(w)|isinf(w)) = [];
y(isnan(w)|isinf(w)) = [];
w(isnan(w)|isinf(w)) = [];

% Are we fitting cumulative Gaussian?
if strcmpi(ft, 'cumulgauss')
    ft = 'a*(1+erf((x-b)/(sqrt(2)*c)))';
    % Fit function robustly
    [oF mR] = fit(x, y, ft, 'weight', w, 'startpoint', [max(m)/2 mean(x) range(m)]);
else
    % Fit function robustly
    [oF mR] = fit(x, y, ft, 'weight', w);
end

% Return fitted parameters
mP = coeffvalues(oF); 
