function [D, gx, gy] = samsrf_cometplot(SrfDv, ValDv, SrfIv, ValIv, Limits, Roi, Threshold, Cmap)
%
% [D, gx, gy] = samsrf_cometplot(SrfDv, ValDv, SrfIv, ValIv, Limits, [Roi='', Threshold=0.01, Cmap='gray'])
%
% Creates a comet plot of the ValDv data in SrfDv against the ValIv data in SrfIv.
%  This can be used as an alternative to binned plots. 
%  The comet plots from individuals can be averaged at the group level.
%  One could also fit a curve to the maximum densities.
%
% SrfDv/SrfIv:  Srf structures. When plotting data within a map, use the -same- Srf.
%
% ValDv/ValIv:  The value name as given by Srf.Values. Returns an error if 
%                more than one entry of that name exists.
%               This can also be 'Eccentricity' or 'Polar', in which case it 
%                may throw up an error if SrfIv.Data(2:3,:) are not X0 and Y0.
%               Prefixing ValDv/ValIv by ':' uses unsmoothed data (if possible).
%               Finally, this can also refer to anatomical statistics such as
%                'Area', 'Thickness', or 'Curvature'. (These will always be
%                anatomical so don't have any Srf.Values of the same names!)
%
% Limits:       1x4 vector defining X-minimum, X-maximum, Y-minimum & Y-maximum, respectively.
%               The optional 5th & 6th components define the granularity of the X & Y dimension.
%                This defaults to 40 and 30, respectively.
%
% Cmap:         Optional, defines the colourmap to be used (defaults to 'gray').
%               If this is preceded by '-', the colour map is inverted.
%               Can also contain any 1*3 colour value or colour string (e.g. 'r') 
%                which defines the colour of a contour line plot instead.
%
% Returns in D the raw density matrix & in gx and gy the coordinates for making a contour plot.
%
% 29/04/2021 - Written (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 29/07/2022 - Changed to logarithmic density scale (DSS)
%              Changed default granularity (DSS)
% 04/08/2022 - Returns raw density matrix even through plot is logarithmic (DSS)
%

if length(Limits) < 4
    samsrf_error('Limits have not been defined!');
end
if length(Limits) < 5
    Limits(5) = 50;
end
if length(Limits) == 5
    Limits(6) = 30;
end
if nargin < 6
    Roi = '';
end
if nargin < 7
    Threshold = 0.01;
end
if nargin < 8
    Cmap = 'gray';
end

%% Expand Srfs if necessary
SrfDv = samsrf_expand_srf(SrfDv);    
SrfIv = samsrf_expand_srf(SrfIv);

%% Check compatibility
if size(SrfDv.Data,2) ~= size(SrfIv.Data,2)
    samsrf_error('SrfDv & SrfIv are not the same mesh!');
end

%% Retrieve data
% Loop thru dependent & independent variable
for i_var = 1:2
    if i_var == 1
        Srf = SrfDv;
        Val = ValDv;
        ValLab = 'ValDv';
    elseif i_var == 2
        Srf = SrfIv;
        Val = ValIv;
        ValLab = 'ValIv';
    end       
    
    % Use unsmoothed data?
    RawLabel = '';
    if Val(1) == ':'
        Val = Val(2:end);
        if isfield(Srf, 'Raw_Data')
            Srf.Data = Srf.Raw_Data;
            RawLabel = ' (unsmoothed)';
        end
    end

    % Which data to use?
    if strcmpi(Val, 'Eccentricity')
        % Eccentricity
        Data = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    elseif strcmpi(Val, 'Polar')
        % Polar angle
        Data = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    elseif strcmpi(Val, 'Area')
        % Surface area
        Data = Srf.Area;
    elseif strcmpi(Val, 'Thickness')
        % Cortical thickness
        Data = Srf.Thickness;
    elseif strcmpi(Val, 'Curvature')
        % Curvature
        Data = Srf.Curvature;
    else
        %% Anything else
        loc = strcmpi(Srf.Values, Val);
        if sum(loc) > 1
            samsrf_error([ValLab ' ' Val ' is ambiguous!']);
        elseif sum(loc) == 0
            samsrf_error([ValLab ' ' Val ' does not exist!'])
        end
        Data = Srf.Data(loc,:);
    end
        
    if i_var == 1
        DataDv = Data;
        ValDv = Val;
        RawLabelDv = RawLabel;
    elseif i_var == 2
        DataIv = Data;
        ValIv = Val;
        RawLabelIv = RawLabel;
    end
end

%% Filter data & restrict to ROI
GoF = true(1,size(SrfIv.Data,2)); % Goodness-of-fit placeholder
% Is 1st row R^2 (usually in pRF maps only)?
if strcmpi(SrfIv.Values{1}, 'R^2') || strcmpi(SrfIv.Values{1}, 'nR^2')
    GoF(SrfIv.Data(1,:) <= Threshold) = false; % Unlabel bad fits
else
    % Warn user if 1st row isn't R^2
    samsrf_disp('WARNING: 1st data row isn''t R^2 so using all good data... Do you really want this?');
end

% Load ROI label
if ~isempty(Roi)
    % ROI vertex indeces
    RoiVx = samsrf_loadlabel(Roi);
    if isnan(RoiVx)
        samsrf_error(['Could not load ROI ' Roi '!']);
    end
    % Label ROI vertices
    RoiLab = false(1,size(SrfIv.Data,2));
    RoiLab(RoiVx) = true;
else
    % No ROI so label all
    RoiLab = true(1,size(SrfIv.Data,2));
end
% Filtered & restricted vectors
DataDv = DataDv(GoF & RoiLab); % Dependent variable
DataIv = DataIv(GoF & RoiLab); % Independent variable

%% Grid parameters
cx = Limits(5);
cy = Limits(6);
rx = range(Limits(1:2))/(cx*2);
ry = range(Limits(3:4))/(cy*2);
[gx,gy] = meshgrid(Limits(1)+rx:rx*2:Limits(2)-rx, Limits(3)+ry:ry*2:Limits(4)-ry);

%% Determine density
D = NaN(cy,cx);
for ix = 1:cx
    for iy = 1:cy
        D(iy,ix) = sum(DataIv>=gx(iy,ix)-rx & DataIv<gx(iy,ix)+rx & DataDv>=gy(iy,ix)-ry & DataDv<gy(iy,ix)+ry);
    end
end
% Normalize
D = abs(log10(D / max(D(:))));
% Invert Y axis
D = flipud(D);
gy = flipud(gy);

%% Contour plot
if ischar(Cmap) && length(Cmap) > 1
    contourf(gx, gy, D, 100, 'edgecolor', 'none');
    if Cmap(1) == '-'
        Cmap = Cmap(2:end);
        cm = colormap(Cmap);
        colormap(flipud(cm));
    else
        colormap(Cmap);
    end
else
    contour(gx, gy, D, 10, 'color', Cmap);
end

%% Cosmetic changes
% Labels
ylabel([ValDv RawLabelDv]);
xlabel([ValIv RawLabelIv]);
% Reformat title
Roi(strfind(Roi,'_')) = '-';
Roi([strfind(Roi,'/'),  strfind(Roi,'\')]) = ' ';
title(Roi);

%% Return raw density matrix
D = 10.^-D;