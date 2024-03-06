function Mesh = samsrf_backproj_del(Response, pRF_Data, Threshold)
%
% Mesh = samsrf_backproj_del(Response, pRF_Data, [Threshold=[0, 0, Inf]])
%
% Projects the activity values in Response back into visual space as a 
% triangle mesh based on Delaunay tesselation of the pRFs locations in pRF_Data.
% pRF_Data contains the first 4 rows of a pRF map (see backprojection functions).
%
% The optional input Threshold defines the minimal R^2 of the pRFs to be projected 
% and the second and third value define the inner and outer eccentricity to be plotted. 
% Defaults to the most inclusive criteria, i.e. all R^2 and eccentricities.
%
% Returns the mesh data for plotting in a structure where X and Y are column vectors 
% with pRF coordinates, Tri is the output of the Delaunay triangulation, and Resp is 
% a matrix with the responses per vertex in rows and conditions in columns. 
% You can use samsrf_heatmap_del for plotting.
%           
% 20/04/2022 - SamSrf 8 version (DSS)
% 05/03/2024 - Fixed bug when using 32-bit (single) data 
%

if nargin < 3
    Threshold = [];
end
if isempty(Threshold)
    Threshold = [0 0 Inf]; 
elseif length(Threshold) == 1
    Threshold = [Threshold 0 Inf];
elseif length(Threshold) == 2
    Threshold = [Threshold Inf];
end

% Convert into doubles
pRF_Data = double(pRF_Data);

% pRF map data
gof = pRF_Data(1,:); % Goodness of fit
ecc = sqrt(pRF_Data(2,:).^2 + pRF_Data(3,:).^2); % Eccentricity
sigma = pRF_Data(4,:); % pRF size

% Filter vertices
nanresp = isnan(Response(1,:)); % Determine NaNs 
Response = Response(:, gof > Threshold(1) & ecc > Threshold(2) & ecc < Threshold(3) & sigma > 0 & nanresp == 0); % Throw out NaNs 
pRF_Data = pRF_Data(:, gof > Threshold(1) & ecc > Threshold(2) & ecc < Threshold(3) & sigma > 0 & nanresp == 0); % Throw out NaNs 

% Output structure  
Mesh = struct;
Mesh.X = pRF_Data(2,:)';
Mesh.Y = pRF_Data(3,:)';
Mesh.Tri = delaunay(Mesh.X, Mesh.Y); % Delaunay triangulation
Mesh.Resp = Response';

