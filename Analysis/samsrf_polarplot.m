function Res = samsrf_polarplot(Srf, Roi, Scale, Val, Thr)
%
% Res = samsrf_polarplot(Srf, [Roi='', Scale=10, Val='Sigma', Thr=[0.01 0 Inf]])
%
% Plots the pRFs inside a ROI as a scatter plot projected to the visual field. 
%
% ROI can be a string pointing to a ROI label or a vector of vertices.
%
% The X and Y coordinates are correct, but the exact size of the pRFs must 
% be calibrated because of the incredibly stupid way Matlab handles marker sizes. 
%
% The optional input argument Scale defines the number of points for a
% symbol of size 1. This is necessary if you want the size of symbols to
% truely reflect the sigma of the vertex. This depends on the axis limits 
% you choose, on the size of the figure, and most likely also differs for
% different computers. To calibrate it, set this to a negative number:
%
%   samsrf_polarplot(Srf,'',-15,'Sigma');
%
% This will plot a red circle with the (absolute) width you defined. You
% need to change this number until the circle fits snugly into the small 
% square around the origin defined by the dotted lines. This is the 
% calibrated width. When using this as a positive number will plot the data 
% with this width. 
%
% If you subsequently want to rescale the axis limits, you will have to 
% multiply the calibrated width by that factor: for example, if the axes
% go up to 25 deg and you want to rescale the axis so they only go up to 
% 5 deg, you will need to multiply the width with 5.
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
% 19/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 2
    Roi = '';
end
if nargin < 3
    Scale = 10;
end
if nargin < 4
    Val = 'Sigma';
end
if nargin < 5
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
        error(['Value ' Val ' does not exist!']);
    end
    if size(D,1) > 1
        error(['Value ' Val ' is ambiguous!']);
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
h = scatter(X, Y, (S*abs(Scale)).^2, D, 'filled');
alpha(h, 0.1);
Ecc = max([abs(X)+max(S) abs(Y)+max(S)]);
if Scale < 0
    scatter(0, 0, abs(Scale).^2, 'r');
end
axis([-Ecc Ecc -Ecc Ecc]);
set(gca, 'fontsize', 12);
hold on
if Scale < 0
    line(xlim, +[1 1], 'color', 'k', 'linestyle', ':');
    line(xlim, -[1 1], 'color', 'k', 'linestyle', ':');
    line(+[1 1], ylim, 'color', 'k', 'linestyle', ':');
    line(-[1 1], ylim, 'color', 'k', 'linestyle', ':');
end
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