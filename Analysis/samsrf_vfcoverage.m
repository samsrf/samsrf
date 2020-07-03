function [Img, Mov] = samsrf_vfcoverage(Srf, Ecc, Roi, R2Thrsh, Clipping, Raw)
%
% [Img, Mov] = samsrf_vfcoverage(Srf, Ecc, [Roi='', R2Thrsh=0.05, Clipping=1, Raw])
%
% Creates an image showing the visual field coverage of a particular ROI. 
% Also plots the image & overlays the pRF centres as black dots.
%
% ROI is a string to a ROI label. If it is undefined or empty, all vertices are displayed.
%   If ROI is a vector it instead directly refers to the ROI vertices.
%
% Ecc defines the eccentricity of the mapping stimulus. 
%   If this is a vector, the second & third element define the 
%   eccentricity range (min, max) to include in the coverage plot.
%   This range defaults to [0 2*Ecc(1)].
%
% R2Thrsh defines the R^2 threshold (defaults = 0.05). 
%
% Clipping defines the value above which all pixels are coloured the same 
%   (default = Inf, i.e. no clipping).
%
% Raw toggles whether raw data are used.
%
% The optional output Mov contains all the pRF profiles as individual frames.
%
% 01/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 3
    Roi = '';
    R2Thrsh = 0.05;
    Clipping = Inf;
    Raw = false;
elseif nargin < 4
    R2Thrsh = 0.05;
    Clipping = Inf;
    Raw = false;
elseif nargin < 5
    Clipping = Inf;
    Raw = false;
elseif nargin < 6
    Raw = false;
end

%% Eccentricity range undefined?
if length(Ecc) == 1
    Ecc = [Ecc 0];
end
if length(Ecc) == 2
    Ecc = [Ecc 2*Ecc];
end

%% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

%% Use raw data?
if Raw
    if isfield(Srf, 'Raw_Data')
        Srf.Data = Srf.Raw_Data;
    end
end

%% Mask out non-ROI vertices
% If ROI input is a string
if ~isempty(Roi) && ischar(Roi)
    rver = samsrf_loadlabel(Roi);
    if isnan(rver)
        error(['ROI ' Roi ' not found!']);
    end
end
% Roi provided as vector
if ~ischar(Roi)
    rver = Roi;
    Roi = 'Vertices provided directly';
end
% If ROI undefined
if isempty(Roi) 
    rver = 1:size(Srf.Vertices,1);
    Roi = 'All vertices';
end
% Limit to ROI
Srf.Data = Srf.Data(:,rver);

%% Remove rubbish
sv = Srf.Data(1,:) > R2Thrsh & Srf.Data(4,:) > 0 & sqrt(Srf.Data(2,:).^2+Srf.Data(3,:).^2) > Ecc(2) & sqrt(Srf.Data(2,:).^2+Srf.Data(3,:).^2) < Ecc(3);
Srf.Data = Srf.Data(:,sv);

%% pRF centre positions
Xpos = Srf.Data(2,:);
Ypos = Srf.Data(3,:);

%% Initialise the image matrix
Mov = zeros(200, 200, size(Srf.Data,2));

%% Detect type of data (pRF, DoG) and normalise data by eccentricity
if sum(strcmpi(Srf.Values, 'Surround'))
    DataType = 'DoG';
    Srf.Data(2:6, :) = Srf.Data(2:6, :) / Ecc(1);
else
    DataType = 'pRF';
    Srf.Data(2:4, :) = Srf.Data(2:4, :) / Ecc(1);
end

%% Loop through all vertices
for v = 1:size(Srf.Data, 2)
    if strcmpi(DataType, 'pRF')
        Mov(:,:,v) = prf_gaussian_rf(Srf.Data(2,v), Srf.Data(3,v), Srf.Data(4,v));
    elseif strcmpi(DataType, 'DoG')
        Mov(:,:,v) = prf_dog_rf(Srf.Data(2,v), Srf.Data(3,v), Srf.Data(4,v), Srf.Data(5,v), Srf.Data(6,v));
    end
end
  
%% Display the images
figure('Name', 'Visual field coverage');
g = Ecc(3)/100;
EccGrd = [-Ecc(3):g:-g g:g:Ecc(3)];
[Xc,Yc] = meshgrid(EccGrd, EccGrd);
Yc = flipud(Yc); % Because Matlab is stupid...
Img = mean(Mov,3); 
if isinf(Clipping)
    Img = Img / max(Img(:));
    Clipping = 1;
else
    Img = Img / Clipping;
    Img(Img > 1) = 1;
end
contourf(Xc, Yc, Img, 100, 'linestyle', 'none');
hold on
scatter(Xpos, Ypos, '.k');
line(xlim, [0 0], 'color', 'k', 'linestyle', '--');
line([0 0], ylim, 'color', 'k', 'linestyle', '--');
set(gca, 'fontsize', 20);
xlabel('Horizontal coordinate (deg)');
ylabel('Vertical coordinate (deg)');
if ischar(Roi)
    [~,RoiName] = fileparts(Roi);
    title(RoiName);
end
axis square
cb = colorbar;
caxis([0 1]);
set(cb, 'ytick', 0:.25:1, 'yticklabel', 0:Clipping/4:Clipping);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
