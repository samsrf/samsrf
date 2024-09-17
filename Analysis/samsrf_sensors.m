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
%       '2D':    Two-dimensional plot of sensors only
%       'Ball':  Three-dimensional distribution plot on a sphere 
%       'Head':  Three-dimensional scalp distribution plot on head model
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
% 03/10/2023 - Removed overly verbose defaults message (DSS)
%              Added hold off to prevent overloading plots (DSS)
%              R^2 maps now use eccentricity colour scheme (DSS)
% 05/10/2023 - Scalp distribution function now takes flexible input (DSS)
%              Can now also plot sphere & head maps (DSS)
% 31/10/2023 - Polar angle now calculates circular mean as it should (DSS)
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
if (PlotType == '3' ||  PlotType == 'B' || PlotType == 'H') && ~isfield(Srf, 'Sphere')
    samsrf_disp('WARNING: No sphere coordinates included - using 2D plot instead...');
    PlotType = '2';
end

%% Check if M/EEG data
if ~strcmpi(Srf.Hemisphere, 'eeg')
    samsrf_error('This Srf does not contain MEG/EEG data!');
end

%% Default colour schemes
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
    Data = Srf.Data(2:3,:);
%     Data = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
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
    elseif strcmpi(Map, 'R^2') || strcmpi(Map, 'nR^2') 
        Cmap = colormap(def_cmap_eccen(2:end));
        if def_cmap_eccen(1) == '-'
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
        samsrf_error([Map ' is ambiguous!']);
    elseif sum(loc) == 0
        samsrf_error([Map ' does not exist!'])
    end
    Data = Srf.Data(loc,:);
end

%% Threshold data
Data(Srf.Data(1,:) <= R2Thr) = NaN;

%% Plot values
if PlotType == 'H' || PlotType == 'B' || PlotType == '3'
    if PlotType == '3'
        % 3D scatter plot
        if size(Data,1) > 1
            % Circular data
            Data = atan2(Data(2,:), Data(1,:)) / pi*180;
        end
        scatter3(Srf.Sphere(:,1), Srf.Sphere(:,2), Srf.Sphere(:,3), 80, Data, 'filled', 'markeredgecolor', 'k');
    elseif PlotType == 'B'
        % Sphere distribution plot
        PlotScalpDist(Srf.Sphere, Data);
    else
        % Head distribution plot
        Head = load('ScalpHeadModel');
        PlotScalpDist(Srf.Sphere, Data, Head);
    end
elseif PlotType == 'S' || PlotType == '2' 
    if PlotType == '2'
        % 2D scatter plot
        if size(Data,1) > 1
            % Circular data
            Data = atan2(Data(2,:), Data(1,:)) / pi*180;
        end
        scatter(Srf.Vertices(:,1), Srf.Vertices(:,2), 60, Data, 'filled', 'markeredgecolor', 'k');
    else
        % Scalp distribution plot
        PlotScalpDist(Srf.Vertices, Data);
    end
    axis square
    axis([-.8 .8 -.8 .8]);
else
    samsrf_error('Invalid plot type specified!');
end
% Colour scheme
colormap(Cmap);
if strcmpi(Map, 'R^2') || strcmpi(Map, 'Eccentricity') || strcmpi(Map, 'Sigma') ...
                       || strcmpi(Map, 'Fwhm') || strcmpi(Map, 'Centre') || strcmpi(Map, 'Surround') 
    caxis([0 1]*max(Data)); % Scale colour scheme
elseif strcmpi(Map, 'Polar') 
    caxis([-180 180]); % Scale colour scheme
else
    caxis([-1 1]*max(abs(Data))); % Scale colour scheme
end
cb = colorbar;
set(get(cb,'label'), 'string', Map); % Colour bar label
title([Map ': t = ' num2str(TimePt)]);

%% Plot scalp distribution
function PlotScalpDist(XYZ, Data, Head)

% Number of neighbours to average
k = 8; 

% Sensor coordinates
sX = XYZ(:,1);
sY = XYZ(:,2);
if size(XYZ,2) > 2
    sZ = XYZ(:,3);
else
    sZ = [];
end

if isempty(sZ)
    %% 2D scalp plot 

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
                Wts(isnan(Data(1,:))) = NaN; % Set bad sensors to NaN
                [~,sx] = sort(Wts, 'descend'); % Determine ascending order
                if size(Data,1) == 1
                    % Linear data
                    Mu(r,c) = nansum(Data(sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k))); % Weighted average of k nearest neighbours            
                else
                    % Circular data
                    Mu(r,c) = atan2(nansum(Data(2,sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k))), nansum(Data(1,sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k)))) / pi*180; % Weighted circular average of k nearest neighbours            
                    Data = atan2(Data(2,:), Data(1,:)) / pi*180; 
                end
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
    
    % Plot sensors
    scatter(sX, sY, 60, Data, 'filled', 'markeredgecolor', 'w');
    scatter(sX, sY, 30, 'k', 'linewidth', 2);
    axis off
    set(gcf, 'color', 'w');
else
    %% 3D scalp plot
    
    if nargin > 2
        % Use head model
        X = Head.Vertices(:,1);
        Y = Head.Vertices(:,2);
        Z = Head.Vertices(:,3);
        Tri = Head.Faces;
        Scale = 1.07;
    else
        % Create sphere
        [X,Y,Z] = sphere(100); % Vertices of sphere
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        warning off
        Tri = delaunay(X, Y, Z);
        warning on
        Scale = 0.95;
    end
    
    % Transparency 
    Alpha = Z*4 + 2; 
    Alpha(Alpha < 0) = 0;
    Alpha(Alpha > 1) = 1;
    
    % Grid for interpolation
    Mu = NaN(size(X)); % Output matrix

    % Loop thru vertices
    for i = 1:size(X,1)
        Euc = sqrt((X(i)-sX).^2 + (Y(i)-sY).^2 + (Z(i)-sZ).^2); % Euclidean distance of sampling point from each sensor
        Wts = 0.5 - Euc; % Turn into weights
        Wts(Wts < 0) = 0; % No negative weights
        Wts(isnan(Data(1,:))) = NaN; % Set bad sensors to NaN
        [~,sx] = sort(Wts, 'descend'); % Determine ascending order
        if size(Data,1) == 1
            % Linear data
            Mu(i) = nansum(Data(sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k))); % Weighted average of k nearest neighbours            
        else
            % Circular data
            Mu(i) = atan2(nansum(Data(2,sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k))), nansum(Data(1,sx(1:k)).*Wts(sx(1:k))') / nansum(Wts(sx(1:k)))) / pi*180; % Weighted circular average of k nearest neighbours            
            Data = atan2(Data(2,:), Data(1,:)) / pi*180; 
        end
    end

    % Plot scalp map
    patch('Vertices', [X Y Z]*Scale, 'Faces', Tri, 'FaceVertexCData', Mu, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceVertexAlphaData', Alpha, 'FaceAlpha', 'interp');         
    daspect([1 1 1]); % Correct aspect ratio
    rotate3d
    
    % Plot sensors
    hold on
    scatter3(sX, sY, sZ, 80, Data, 'filled', 'markeredgecolor', 'w');
    scatter3(sX, sY, sZ, 50, 'k', 'linewidth', 2);
    scatter3(0, 1.3, 0, 90, 'k^', 'markerfacecolor', [.5 .5 .5], 'linewidth', 2);
    axis off
    set(gcf, 'color', 'w');
end 
    
hold off











