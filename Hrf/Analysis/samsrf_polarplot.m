function [Res, h] = samsrf_polarplot(Srf, Roi, Val, Thr)
%
% [Res, h] = samsrf_polarplot(Srf, [Roi='', Val='Sigma', Thr=[0.01 0 Inf]])
%
% Plots the pRFs inside a ROI as a scatter plot projected to the visual field. 
%
% ROI can be a string pointing to a ROI label or a vector of vertices.
%
% Note: The X and Y coordinates are correct, but the exact size of the pRFs must 
% be calibrated because of the incredibly stupid way Matlab handles marker sizes. 
% This function does this automatically but you will need to do this again
% if you rescale the figure. To do so you simply need to run scatter_size(h) 
% where h is the handle to the scatter plot (second output of this function).
%
% By default the colour of each circle also codes for the sigma though and 
% can be read off the colourbar. The optional input argument Val defines 
% which of the data fields to display in the colour code. This can also 
% be 'Polar' or 'Eccentricity' although both of these are obviously not
% terribly meaningful in this plot... You can also plot the cortical
% surface 'Area', Thickness', or 'Curvature'. Finally, if the first letter 
% is ':' it plots the raw unsmoothed parameters (if available).
%
% The optional input Thr defines the R^2 threshold of the data. Optionally,
% you can also define an eccentricity range in Thr(2:3). 
%
% The function returns a matrix containing the data: [X0 Y0 Sigma Value]
%   where Value stands for the value you plotted in the colour code.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    Roi = '';
end
if nargin < 3
    Val = 'Sigma';
end
if nargin < 4
    Thr = 0.01;
end

if length(Thr) == 1
    Thr = [Thr 0 Inf];
elseif length(Thr) == 2
    Thr = [Thr Inf];    
end

%% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

%% Use raw data?
if Val(1) == ':'
    Val = Val(2:end);
    if isfield(Srf, 'Raw_Data')
        Srf.Data = Srf.Raw_Data;
    end
end

%% Mask out non-ROI vertices
% If ROI input is a string
if ~isempty(Roi) && ischar(Roi)
    rver = samsrf_loadlabel(Roi);
    if isnan(rver)
        samsrf_error(['ROI ' Roi ' not found!']);
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

%% Extract data
E = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); % All eccentricities
% Now get only desired data 
R = Srf.Data(1,:) > Thr(1) & Srf.Data(4,:) > 0 & E >= Thr(2) & E <= Thr(3);
E = sqrt(Srf.Data(2,R).^2 + Srf.Data(3,R).^2);
P = atan2(Srf.Data(3,R), Srf.Data(2,R));
S = Srf.Data(4,R);
% Cartesian space
[X,Y] = pol2cart(P,E);

% Value name
Val = lower(Val);
Val(1) = upper(Val(1));
% Determine data
if strcmpi(Val, 'Polar')
    D = P/pi*180;
elseif strcmpi(Val, 'Eccentricity')
    D = E;
elseif strcmpi(Val, 'Area')
    D = Srf.Area(R);
elseif strcmpi(Val, 'Thickness')
    D = Srf.Thickness(R);
elseif strcmpi(Val, 'Curvature')
    D = Srf.Curvature(R);
else
    D = Srf.Data(strcmpi(Srf.Values, Val),R);
    if isempty(D)
        samsrf_error(['Value ' Val ' does not exist!']);
    end
    if size(D,1) > 1
        samsrf_error(['Value ' Val ' is ambiguous!']);
    end
end

%% Sort by ascending absolute value
[~,sx] = sort(abs(D));
X = X(sx);
Y = Y(sx);
S = S(sx);
D = D(sx);

%% Polar plot
figure('name', Val);
h = scatter(X, Y, S, D, 'filled');
scatter_size(h);
alpha(h, 0.1);
Ecc = max([abs(X)+max(S) abs(Y)+max(S)]);
axis([-Ecc Ecc -Ecc Ecc]);
set(gca, 'fontsize', 12);
hold on
line(xlim, [0 0], 'color', 'k')
line([0 0], ylim, 'color', 'k')
axis square
xlabel('Horizontal coordinate (deg)');
ylabel('Vertical coordinate (deg)');
if ischar(Roi)
    [~,RoiName] = fileparts(Roi);
    title(RoiName);
end
colorbar

%% Return data
Res = [X; Y; S; D];