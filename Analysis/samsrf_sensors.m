function samsrf_sensors(Srf, Map, R2Thr, TimePt, IsFlat)
%
% samsrf_sensors(Srf, Map, [R2Thr=0], [TimePt=1], [IsFlat=false])
%
% Plots the data in Map (from Srf.Values) for each sensor. 
% This requires a Srf file with EEG/MEG data. 
%
% R2Thr is the R^2 threshold, which defaults to 0.
%
% TimePt is the time point to plot, which defaults to 1.
%
% IsFlat toggles whether to plot a 3D (default) or 2D arrangement.
%
% 30/08/2023 - Written (DSS)
% 07/09/2023 - Axis now square in flat mode (DSS)
%

if nargin < 3
    R2Thr = 0;
end
if nargin < 4
    TimePt = 1;
end
if nargin < 5
    IsFlat = false;
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
end
% Which data to use?
if strcmpi(Map, 'Eccentricity')
    % Eccentricity
    Data = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);   
    Cmap = colormap(def_cmap_eccen(2:end));
    if def_cmap_eccen(1) == '-'
        Cmap = flipud(Cmap);
    end
elseif strcmpi(Map, 'Polar')
    % Polar angle
    Data = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    Cmap = colormap(def_cmap_angle(2:end));
    if def_cmap_angle(1) == '-'
        Cmap = flipud(Cmap);
    end
else
    % Is pRF size?
    if strcmpi(Map, 'Sigma')
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
Srf.Vertices = Srf.Vertices(Srf.Data(1,:) > R2Thr,:);
Data = Data(Srf.Data(1,:) > R2Thr);

%% Plot values
if IsFlat
    scatter(Srf.Vertices(:,1), Srf.Vertices(:,2), 120, Data, 'filled', 'markeredgecolor', 'k');
    axis([-1.2 1.2 -1.2 1.2]);
    axis square
else
    scatter3(Srf.Vertices(:,1), Srf.Vertices(:,2), Srf.Vertices(:,3), 120, Data, 'filled', 'markeredgecolor', 'k');
end
% Colour scheme
colormap(Cmap);
if strcmpi(Map, 'R^2') || strcmpi(Map, 'Eccentricity') || strcmpi(Map, 'Sigma')
    caxis([0 1]*max(Data)); % Scale colour scheme
else
    caxis([-1 1]*max(abs(Data))); % Scale colour scheme
end
cb = colorbar;
set(get(cb,'label'), 'string', Map); % Colour bar label
title([Map ': t = ' num2str(TimePt)]);