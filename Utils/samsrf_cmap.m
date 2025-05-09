function rgb = samsrf_cmap(MapName, Nrows)
%
% rgb = samsrf_cmap(MapName, [Nrows=256])
%
% Reads in the colour look-up table in the file [MapName].csv & converts it
% if necessary to have Nrows rows. This allows you to define colour maps
% without having Matlab functions & is essential for standalone app.
%
% 19/09/2024 - Written (DSS)
%

if nargin < 2
    Nrows = 256;
end

% If sign undefined
if MapName(1) ~= '+' && MapName(1) ~= '-'
    MapName = ['+' MapName];
end

% Load colour table
if isdeployed 
    % Running in deployed SamSrfX app
    global SamSrfXPath
    rgb = readmatrix([SamSrfXPath filesep 'Colours' filesep MapName(2:end) '.csv']);
else
    % Running in Matlab command window
    rgb = readmatrix([MapName(2:end) '.csv']);
end

% Inverting colour map?
if MapName(1) == '-'
    rgb = flipud(rgb);
end

% Rescale?
if Nrows ~= size(rgb,1)
    rgb = imresize(rgb, [Nrows 3]);
end