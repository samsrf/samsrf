function Alpha = circmean(X)
% Alpha = circmean(X)
%
% Calculates the circular mean of the angles (in degrees) in X.

X = X / 180 * pi; % Convert to radians
n = sum(~isnan(X)); % Number of observations
Alpha = atan2(nansum(sin(X)) ./ n, nansum(cos(X)) ./ n); % Circular mean
Alpha = Alpha / pi * 180; % Convert back into degrees