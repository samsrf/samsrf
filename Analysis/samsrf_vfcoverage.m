function [Img, Mov] = samsrf_vfcoverage(Srf, Ecc, Roi, R2Thrsh, Clipping, Raw)
%
% [Img, Mov] = samsrf_vfcoverage(Srf, Ecc, [Roi='', R2Thrsh=0.05, Clipping, Raw])
%
% Creates an image showing the visual field coverage of a particular ROI. 
%
% ROI is a string to a ROI label. If it is undefined or empty, all vertices are displayed.
% If ROI is a scalar vector it instead directly refers to the ROI vertices.
%
% Ecc defines the eccentricity of the mapping stimulus. R2Thrsh defines the 
% R^2 threshold (defaults = 0.1). Clipping defines the value above which
% all pixels are coloured the same (default = Inf, i.e. no clipping).
% Finally, Raw toggles whether raw data are used.
%
% The optional output Mov contains all the pRF profiles as individual frames.
%
% 10/08/2018 - SamSrf 6 version (DSS)
% 11/08/2018 - Fixed bug with switch but need to add support 
%               for flexible function handles (DSS)
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
sv = Srf.Data(1,:) > R2Thrsh & Srf.Data(4,:) > 0 & sqrt(Srf.Data(2,:).^2+Srf.Data(3,:).^2) < 2*Ecc;
Srf.Data = Srf.Data(:,sv);

%% Initialise the image matrix
Mov = zeros(200, 200, size(Srf.Data,2));

%% Detect type of data (pRF, DoG) and normalise data by eccentricity
if sum(strcmpi(Srf.Values, 'Surround'))
    DataType = 'DoG';
    Srf.Data(2:6, :) = Srf.Data(2:6, :) / Ecc;
else
    DataType = 'pRF';
    Srf.Data(2:4, :) = Srf.Data(2:4, :) / Ecc;
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
Img = mean(Mov,3); 
Img = Img(26:175,26:175);
Img = flipud(Img);
if isinf(Clipping)
    Img = Img / max(Img(:));
else
    Img = Img / Clipping;
    Img(Img > 1) = 1;
end
contourf(Img, 50, 'linestyle', 'none');
colormap jet 
set(gca, 'fontsize', 12);
xlabel('Horizontal coordinate (deg)');
ylabel('Vertical coordinate (deg)');
if ischar(Roi)
    [~,RoiName] = fileparts(Roi);
    title(RoiName);
end
axis square
axis([1 150 1 150]);
set(gca, 'xtick', 0:25:150, 'xticklabel', -Ecc*1.5:Ecc/2:Ecc*1.5);
set(gca, 'ytick', 0:25:150, 'yticklabel', -Ecc*1.5:Ecc/2:Ecc*1.5);
colorbar
caxis([0 1]);
