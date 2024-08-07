function scatter_size(h)
%
% scatter_size(h)
%
% Adjusts the size of symbols in the scatter plot in handle h to ensure
% that Sigma corresponds to the radius of the circles in the plot.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

U = get(gca, 'Units'); % What are current units of axes?
set(gca, 'Units', 'Points'); % Set units to points
Ax = get(gca, 'Position'); % Axes dimensions in points
set(gca, 'Units', U); % Set back to original units
W = 2 / diff(xlim) * Ax(3); % Marker width in points
S = get(h, 'SizeData'); % Extract size data
set(h, 'SizeData', S * W^2); % Set size data in scatter plot

