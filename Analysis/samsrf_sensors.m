function samsrf_sensors(Srf, Map, R2Thr, TimePt, PlotType)
%
% samsrf_sensors(Srf, Map, [R2Thr=0], [TimePt=1], [PlotType='Scalp'])
%
% Plots the data in Map (from Srf.Values) for each sensor. 
% This requires a Srf file with EEG/MEG data. 
%
% R2Thr is the R^2 threshold, which defaults to 0.
%
% TimePt is the time point to plot, which defaults to 1.
%
% PlotType defines the arrangement of the plot:
%       'Scalp': Two-dimensional scalp distribution plot 
%                (This is default but requires coordinates in Srf.Vertices) 
%       '2D':    Flat plot of 3D coordinates without depth dimension
%       '3D':    Three-dimensional plot of sensors/electrodes 
%
% 30/08/2023 - Written (DSS)
% 07/09/2023 - Axis now square in flat mode (DSS)
% 08/09/2023 - Corrected colour scheme for pRF sizes (DSS)
% 19/09/2023 - pRF size measures now zero-bound (DSS)
% 26/09/2023 - Added new default of scalp distribution plots (DSS)
%              Adapted colour schemes (DSS)
% 28/09/2023 - Changed Srf structure for scalp distribution plots (DSS)
%              Sub-threshold sensors are now plotted but set to 0 (DSS)
%

if nargin < 3
    R2Thr = 0;
end
if nargin < 4
    TimePt = 1;
end
if nargin < 5
    PlotType = 'Scalp';
end

%% Determine plot type
PlotType = upper(PlotType(1));
if PlotType == '3' && ~isfield(Srf, 'Sphere')
    warning('No sphere coordinates included - using 2D plot instead...');
    PlotType = '2';
end

%% Check if M/EEG data
if ~strcmpi(Srf.Hemisphere, 'eeg')
    error('This Srf does not contain MEG/EEG data!');
end

%% Default colour schemes
disp(['Using defaults in: ' which('SamSrf_defaults.mat')]);
load('SamSrf_defaults.mat');
% Ensure colour maps have sign
if def_cmap_angle(1) ~= '-' && def_cmap_angle(1) ~= '+'
    def_cmap_angle = ['+' def_cmap_angle];
end
if def_cmap_eccen(1) ~= '-' && def_cmap_eccen(1) ~= '+'
    def_cmap_eccen = ['+' def_cmap_eccen];
end
if def_cmap_other(1) ~= '-' && def_cmap_other(1) ~= '+'
    def_cmap_other = ['+' def_cmap_other];
end
if def_cmap_sigma(1) ~= '-' && def_cmap_sigma(1) ~= '+'
    def_cmap_sigma = ['+' def_cmap_sigma];
end

%% Extract data
% Limit to time points 
ts = Srf.TimePts == TimePt;
if ~isempty(ts)
    Srf.Data = Srf.Data(:,ts);
    Srf.Vertices = Srf.Vertices(ts,:);
    if isfield(Srf, 'Sphere')
        Srf.Sphere = Srf.Sphere(ts,:);
    end
end
% Which data to use?
if strcmpi(Map, 'Eccentricity')
    % Eccentricity
    Data = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);   
    Cmap = colormap(def_cmap_eccen(2:end));
    if def_cmap_eccen(1) == '-'
        Cmap = flipud(Cmap);
    end
elseif strcmpi(Map, 'Polar') || strcmpi(Map, 'Phase') || strcmpi(Map, 'Phi')
    % Polar angle
    Data = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    Cmap = colormap(def_cmap_angle(2:end));
    if def_cmap_angle(1) == '-'
        Cmap = flipud(Cmap);
    end
else
    % Is pRF size?
    if strcmpi(Map, 'Sigma') || strcmpi(Map, 'Fwhm') || strcmpi(Map, 'Centre') || strcmpi(Map, 'Surround') 
        Cmap = colormap(def_cmap_sigma(2:end));
        if def_cmap_sigma(1) == '-'
            Cmap = flipud(Cmap);
        end
    else
        Cmap = colormap(def_cmap_other(2:end));
        if def_cmap_other(1) == '-'
            Cmap = flipud(Cmap);
        end
    end
    % Anything else
    loc = strcmpi(Srf.Values, Map);
    if sum(loc) > 1
        error([Map ' is ambiguous!']);
    elseif sum(loc) == 0
        error([Map ' does not exist!'])
    end
    Data = Srf.Data(loc,:);
end

%% Threshold data
Data(Srf.Data(1,:) <= R2Thr) = NaN;

%% Plot values
if PlotType == '3'
    % 3D scatter plot
    scatter3(Srf.Sphere(:,1), Srf.Sphere(:,2), Srf.Sphere(:,3), 80, Data, 'filled', 'markeredgecolor', 'k');
elseif PlotType == 'S' || PlotType == '2' 
    if PlotType == '2'
        % 2D scatter plot
        scatter(Srf.Vertices(:,1), Srf.Vertices(:,2), 60, Data, 'filled', 'markeredgecolor', 'k');
    else
        % Scalp distribution plot
        PlotScalpDist(Srf.Vertices(:,1), Srf.Vertices(:,2), Data);
    end
    axis square
    axis([-.8 .8 -.8 .8]);
else
    error('Invalide plot type specified!');
end
% Colour scheme
colormap(Cmap);
if strcmpi(Map, 'R^2') || strcmpi(Map, 'Eccentricity') || strcmpi(Map, 'Sigma') ...
                       || strcmpi(Map, 'Fwhm') || strcmpi(Map, 'Centre') || strcmpi(Map, 'Surround') 
    caxis([0 1]*max(Data)); % Scale colour scheme
else
    caxis([-1 1]*max(abs(Data))); % Scale colour scheme
end
cb = colorbar;
set(get(cb,'label'), 'string', Map); % Colour bar label
title([Map ': t = ' num2str(TimePt)]);

%% Plot scalp distribution
function PlotScalpDist(sX, sY, Data)

k = 8; % Number of neighbours to average

% Grid for interpolation
xl = max(abs(sX)) * 1.1; % X limit 
yl = max(abs(sY)) * 1.1; % Y limit
rl = max(sqrt(sX.^2 + sY.^2)); % Rho limit
[X,Y] = meshgrid(-xl:.01:xl, fliplr(-xl:.01:xl)); % Position grid 
Rho = sqrt(X.^2 + Y.^2); % Distance from centre
Mu = NaN(size(X)); % Output matrix

% Loop thru X coordinates
for c = 1:size(X,2)
    % Loop thru Y coordinates
    for r = 1:size(Y,1)
        if Rho(r,c) < rl
            Euc = sqrt((X(r,c)-sX).^2 + (Y(r,c)-sY).^2); % Euclidean distance of sampling point from each sensor
            Wts = 0.3 - Euc; % Turn into weights
            Wts(Wts < 0) = 0; % No negative weights
            Wts(isnan(Data)) = NaN; % Set bad sensors to NaN
            [~,sx] = sort(Wts, 'descend'); % Determine ascending order
            Mu(r,c) = nansum(Data(sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k))); % Weighted average of k nearest neighbours            
        end
    end
end

% Plot head 
hold on
[hx,hy] = pol2cart((0:360)/180*pi, rl);
plot(hx, hy, 'k-', 'linewidth', 2);
% Plot nose
[nx,ny] = pol2cart([100 90 80]/180*pi, [1 1.25 1]*rl);
plot(nx, ny, 'k-', 'linewidth', 2);
% Plot ears
[ex,ey] = pol2cart((0:360)/180*pi, .1 + (.1*sind(0:360)));
plot(ex+rl, ey-.1, 'k-', 'linewidth', 2); % Left ear
plot(ex-rl, ey-.1, 'k-', 'linewidth', 2); % Right ear

% Plot scalp map
contourf(X, Y, Mu, 100, 'EdgeColor', 'none');
scatter(sX, sY, 60, Data, 'filled', 'markeredgecolor', 'w');
scatter(sX, sY, 30, Data, 'k', 'linewidth', 2);
axis off
set(gcf, 'color', 'w');












