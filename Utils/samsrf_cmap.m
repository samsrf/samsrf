function rgb = samsrf_cmap(MapName, Nrows)
%
% rgb = samsrf_cmap(MapName, Nrows)
%
% Reads in the colour look-up table in the file [MapName].csv & converts it
% if necessary to have Nrows rows. This allows you to define colour maps
% without having Matlab functions & is essential for standalone app.
%
% 19/09/2024 - Written (DSS)
%

% global SamSrfXPath

% If sign undefined
if MapName(1) ~= '+' && MapName(1) ~= '-'
    MapName = ['+' MapName];
end

% Load colour table
rgb = csvread([MapName(2:end) '.csv']);

% Inverting colour map?
if MapName(1) == '+'
    rgb = flipud(rgb);
end

% Rescale?
if Nrows ~= size(rgb,1)
    rgb = imresize(rgb, [Nrows 3]);
end