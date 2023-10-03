function samsrf_drawpaths(PatchHdl, Paths, Colour, Srf)
%
% samsrf_drawpaths(PatchHdl, Paths, Colour, [Srf])
%
% Draws a list of paths onto a surface map created with samsrf_surf. 
% You can use this if you want to plot paths of different colours on the same map 
% (It would be more complicated adding this to samsrf_surf & I have better things to do...)
% Automatically loads the default path thickness from SamSrf_defaults.
%
%   PatchHdl:   The patch handle to the map returned by samsrf_surf
%   Paths:      String or cell array with the paths (.mat or .label file, or a vector)
%   Colour:     1x3 vector with colour RGB value or Inf to make transparent path
%   Srf:        Srf with surface data - only required if loading .label files
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 03/10/2023 - Removed overly verbose defaults message (DSS)
%

%% Load default parameters?
load('SamSrf_defaults.mat');

%% Ensure cell array
if ~iscell(Paths)
    Paths = {Paths};
end

%% Vector of path vertices
Vs_paths = [];
% Loop thru paths
for i = 1:length(Paths)
    if ~ischar(Paths{i})             
        % If paths contain vertices
        Vs_paths = [Vs_paths; Paths{i}];
    else
        % Is a string
        if strfind(Paths{i}, '.label')
            % Load label
            Label = samsrf_loadlabel(Paths{i}(1:end-6));
            if nargin < 4
                error('Loading labels requires a Srf as input!');
            else
                Vs_paths = [Vs_paths; samsrf_borderpath(Srf, Label)];
            end
        else
            % Load paths
            Vs_paths = [Vs_paths; samsrf_loadpaths(Paths{i})];
        end
    end
end

%% Default path thickness defined?
if exist('def_pathwidth', 'var')
    for i = 1:def_pathwidth-1
        Vs_paths = [Vs_paths; samsrf_neighbours(Vs_paths, Srf.Faces)];
    end
end

%% Draw new paths
rgb = get(PatchHdl, 'FaceVertexCData');
if isinf(Colour)
    % Transparent paths
    rgb(Vs_paths,:) = 1 - rgb(Vs_paths,:);
else
    % Colour defined by RGB
    rgb(Vs_paths,:) = repmat(Colour, length(Vs_paths), 1);
end
set(PatchHdl, 'FaceVertexCData', rgb);
