function samsrf_simvsfit(Srf, Thresholds, SearchSpace)
%
% samsrf_simvsfit(Srf, [Thresholds=[NaN -Inf], SearchSpace=[]])
%
% Plots a comparison of simulated ground truth pRFs & model fits.
% At present this function only works for standard 2D Gaussian pRFs.
%
% Srf contains the model fit of a simulated pRF data set, so this must 
%   contain a Srf.Ground_Truth field. It assumes that Srf.Data(2:3,:) 
%   contains the X and Y coordinates & Srf.Data(4,:) contains Sigma. 
%
%   Plots the position shifts as a quiver graph: the dot symbols denote 
%   the modelled positions and the lines denote the shifts from the ground 
%   truth positions. The dot colours denote the difference in Sigma. 
%   The colours of the lines denote the modelled Beta parameters.
%   A circular mask for where the apertures should be is also included - 
%   so if for some reason your stimulus was different you will need to
%   change this to match your design.
%
% Thresholds(1) defines the ground truth Sigma to restrict the plot to. 
%   A typical simulation would contain a range of spatial positions and a 
%   range of Sigmas. You can restrict the plot to only one simulated Sigma
%   to prevent clutter. Defaults to NaN for no restriction.
%
% Thresholds(2) defines the R^2 threshold of the model fits to include in 
%   the comparison. Defaults to -Inf (includes all).
%
% SearchSpace defines the position parameters of the search grid. 
%   These are plotted on top of the plot as crosses. This must be a 
%   k-by-m matrix where n is the number of parameters in the search grid 
%   and m the number of grid positions. Note, however, only the unique 
%   spatial positions are plotted so it ignores Sigmas, Thetas, or any 
%   other such parameters. You would usually find this search space matrix 
%   in the S variable of your src_*.mat file.
%
% 02/06/2020 - SamSrf 7 version (DSS) 
%

if nargin < 2
    Thresholds = [NaN -Inf];
end
if length(Thresholds) == 1
    Thresholds = [Thresholds -Inf];
end
if nargin < 3
    SearchSpace = [];
end

%% Retrieve data
% Modelled parameters
R2 = Srf.Data(1,:);
mX = Srf.Data(2,:);
mY = Srf.Data(3,:);
mS = Srf.Data(4,:);
mB = Srf.Data(5,:);

% Ground truth posiitonsti
tX = Srf.Ground_Truth(1,:);
tY = Srf.Ground_Truth(2,:);
tS = Srf.Ground_Truth(3,:);

% Difference in Sigma
Delta = mS-tS;

%% Filter data
% Remove bad fits
g = R2 > Thresholds(2);
mX = mX(g); tX = tX(g);
mY = mY(g); tY = tY(g);
mS = mS(g); tS = tS(g);
Delta = Delta(g); 
if isempty(Delta)
    error(['No data with R^2 > ' num2str(Thresholds(2))]);
end
TitlStr = ['R^2 > ' num2str(Thresholds(2))];
% Limit ground truths?
if ~isnan(Thresholds(1))
    % Only this Sigma 
    g = tS == Thresholds(1);
    mX = mX(g); tX = tX(g); % X position
    mY = mY(g); tY = tY(g); % Y position
    mS = mS(g); tS = tS(g); % Sigma (not used at present)
    mB = mB(g); % Beta 
    Delta = Delta(g); % Sigma difference
    TitlStr = {TitlStr; ['Ground truth \sigma = ' num2str(Thresholds(1))]};
else
    TitlStr = {TitlStr; ['All ground truth \sigma']};
end
% If this removed all data
if isempty(Delta)
    error('No data with this ground truth in file!');
end

%% Scales for plot
sS = max(abs([min(Delta) max(Delta)])); % Sigma difference scale
TitlStr{3} = [num2str(min(mB)) ' \leq \beta \geq ' num2str(max(mB))];
mB = mB - min(mB); % Baseline correct Betas
mB = round(mB / max(mB) * 99)+1; % Normalise Betas
mB(isnan(mB)) = 1; % In case poor betas

%% Plot mask
hold on
set(gca, 'color', [1 1 1]*.7);
[sx,sy] = pol2cart((0:360)/180*pi, 1);
h = fill(sx, sy, [1 1 1]*.9);
alpha(h, 0.3);

%% Plot lines
Cmap = hotcold(100); % Beta colour map
for v = 1:length(tX)
    line([tX(v) mX(v)], [tY(v) mY(v)], 'color', Cmap(mB(v),:), 'linewidth', 2);
end

%% Plot modelled positions 
scatter(mX, mY, 30, Delta, 'filled', 'markeredgecolor', 'k');

%% Plot search space
if ~isempty(SearchSpace)
    s3 = unique(SearchSpace(3,:));
    sX = SearchSpace(1,SearchSpace(3,:) == s3(1));
    sY = SearchSpace(2,SearchSpace(3,:) == s3(1));
    scatter(sX, sY, 100, 'k+');
end

%% Cosmetic changes
grid on
cb = colorbar;
colormap hotcold
axis([-2 2 -2 2]);
axis square
set(gca, 'fontsize', 20, 'Clim', [-1 1]*sS);
xlabel('Horizontal position');
ylabel('Vertical position');
title(TitlStr);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
set(get(cb, 'Label'), 'String', '\Delta_{\sigma}');