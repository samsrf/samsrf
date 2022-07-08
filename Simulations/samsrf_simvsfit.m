function samsrf_simvsfit(Srf, Thresholds, SearchSpace, PlotRsq)
%
% samsrf_simvsfit(Srf, [Thresholds=[NaN -Inf], SearchSpace=[], PlotRsq])
%
% Plots a comparison of simulated ground truth pRFs & model fits in Srf.
% At present this function only works for standard 2D Gaussian pRFs.
% These plots really only make sense when you have noise-free simulations.
% When you have noisy simulations use samsrf_simvsfithist instead.
%
% Srf contains the model fit of a simulated pRF data set, so this must 
%   contain a Srf.Ground_Truth field. It assumes that Srf.Data(2:3,:) 
%   contains the X and Y coordinates & Srf.Data(4,:) contains Sigma. 
%
% Plots the position shifts as a quiver graph: the dot symbols denote 
%   the modelled positions and the lines denote the shifts from the ground 
%   truth positions. The dot colours denote the modelled Beta or R^2.
%
%   A circular mask for where the apertures should be is also included - 
%   so if for some reason your stimulus was different (e.g. square?) 
%   you will need to change this to match your design.
%
% An alternative plot is a scatter graph comparing ground truth and 
%   modelled pRF size (sigma) parameters. Each symbol is one pRF. 
%   The colour denotes the beta amplitude parameter.
%
%
% Thresholds(1) defines the ground truth Sigma to restrict the plot to. 
%   A typical simulation would contain a range of spatial positions and a 
%   range of Sigmas. You can restrict the plot to only one simulated Sigma
%   to prevent clutter. This defaults to NaN, which means no restriction.
%   When this is NaN the function also plots a comparison of Sigmas between
%   ground truth and modelled pRFs.
%
% Thresholds(2) defines the R^2 threshold of the model fits to include in 
%   the comparison. Defaults to -Inf (includes all).
%
% SearchSpace defines the position parameters of the search grid. 
%   These are plotted on top of the plot as crosses. This must be a 
%   k-by-m matrix where k is the number of parameters in the search grid 
%   and m the number of grid positions. Note, however, only the unique 
%   spatial positions are plotted so it ignores Sigmas, Thetas, or any 
%   other such parameters. You would usually find this search space matrix 
%   in the S variable of your src_*.mat file.
%
% PlotRsq is a boolean that toogles whether dot colours denote R^2 (true) 
%   or the modelled Betas (false).
%
% 08/07/2022 - Minor bug fixes & adjusted for SamSrf 9(DSS) 
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
if nargin < 4
    PlotRsq = false;
end

%% Retrieve data
% Modelled parameters
R2 = Srf.Data(1,:);
mX = Srf.Data(2,:);
mY = Srf.Data(3,:);
mS = Srf.Data(4,:);
mB = Srf.Data(5,:);
% Plot betas or R^2?
if PlotRsq
    mB = R2;
end

% Ground truth parameters
tX = Srf.Ground_Truth(1,:);
tY = Srf.Ground_Truth(2,:);
tS = Srf.Ground_Truth(3,:);

% Plot dimensions
dims = max(max(abs(Srf.Ground_Truth(1:2,:))));

%% Remove bad fits
g = R2 > Thresholds(2);
mX = mX(g); tX = tX(g);
mY = mY(g); tY = tY(g);
mS = mS(g); tS = tS(g);
mB = mB(g);
if isempty(mB)
    error(['No data with R^2 > ' num2str(Thresholds(2))]);
end
TitlStr = ['R^2 > ' num2str(Thresholds(2))];

%% Scales for sigma plot
sB = max(abs([min(mB) max(mB)])); % Beta difference / R^2 scale
if sB == 0
    sB = 1;
end

%% Plot pRF sizes
if isnan(Thresholds(1))
    figure; hold on 
    line([0 max([tS mS])], [0 max([tS mS])], 'color', 'k');
    scatter(tS, mS, 200, mB, 'filled', 'markeredgecolor', 'k');
    axis square
    cb = colorbar;
    colormap hotcold
    set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
    set(gca, 'fontsize', 20, 'Clim', [-1 1]*sB);
    xlabel('Ground truth \sigma');
    ylabel('Modelled \sigma');
    set(get(cb, 'Label'), 'String', 'Modelled \beta amplitude');
end

%% Limit ground truths?
if ~isnan(Thresholds(1))
    % Only this Sigma 
    g = tS == Thresholds(1);
    mX = mX(g); tX = tX(g); % X position
    mY = mY(g); tY = tY(g); % Y position
    mS = mS(g); tS = tS(g); % Sigma (not used in following plots)
    mB = mB(g); % Beta 
    TitlStr = {TitlStr; ['Ground truth \sigma = ' num2str(Thresholds(1))]};
else
    TitlStr = {TitlStr; 'All ground truth \sigma'};
end
% If this removed all data
if isempty(mB)
    error('No data with this ground truth in file!');
end

%% Scales for quiver plot
sB = max(abs([min(mB) max(mB)])); % Beta difference / R^2 scale
if sB == 0
    sB = 1;
end

%% Plot mask
figure
hold on
set(gca, 'color', [1 1 1]*.7);
[sx,sy] = pol2cart((0:360)/180*pi, 1);
h = fill(sx, sy, [1 1 1]*.9);
alpha(h, 0.3);

%% Plot lines
Cmap = hotcold(100); % Beta colour map
for v = 1:length(tX)
    line([tX(v) mX(v)], [tY(v) mY(v)], 'color', [1 1 1]/2);
end

%% Plot modelled positions 
scatter(mX, mY, 60, mB, 'filled', 'markeredgecolor', 'k');

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
axis([-1.5 1.5 -1.5 1.5]*dims);
axis square
set(gca, 'fontsize', 20, 'Clim', [-1 1]*sB);
xlabel('Horizontal position');
ylabel('Vertical position');
title(TitlStr);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
if PlotRsq
    set(get(cb, 'Label'), 'String', 'Model fit R^2');
else
    set(get(cb, 'Label'), 'String', 'Modelled \beta amplitude');
end