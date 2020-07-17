function Alpha = circmean(X)
% Alpha = circmean(X)
%
% Calculates the circular mean of the angles (in degrees) in X.

X = X / 180 * pi; % Convert to radians
n = length(X); % Number of observations
Alpha = atan2(sum(sin(X)) / n, sum(cos(X)) / n); % Circular mean
Alpha = Alpha / pi * 180; % Convert back into degrees