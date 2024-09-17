function samsrf_simvsfit_hist(Srf, Thresholds, Nreps)
%
% samsrf_simvsfit_hist(Srf, [Thresholds=[NaN -Inf 1], Nreps=1])
%
% Plots a comparison of simulated ground truth pRFs & noisy model fits in Srf.
% It plots position shifts as a 2D density histogram.
% These plots really only make sense when you have noisey simulations.
% When you have noise-free simulations use samsrf_simvsfit instead.
%
% Srf contains the model fit of a simulated pRF data set, so this must 
%   contain a Srf.Ground_Truth field. It assumes that Srf.Data(2:3,:) 
%   contains the X and Y coordinates & Srf.Data(4,:) contains Sigma. 
%
% Plots the position shifts as a 2D density histogram. This is clipped by
%   the 99th percentile to avoid overexposure artifacts. The ground truth
%   positions are plotted as red dots on top of that.
%
%   A circular mask for where the apertures should be is also included - 
%   so if for some reason your stimulus was different (e.g. square?) 
%   you will need to change this to match your design.
%
%
% Thresholds(1) defines the ground truth Sigma to restrict the plot to. 
%   A typical simulation would contain a range of spatial positions and a 
%   range of Sigmas. You can restrict the plot to only one simulated Sigma
%   to prevent clutter. Defaults to NaN for no restriction.
%
% Thresholds(2) defines the R^2 threshold of the model fits to include in 
%   the comparison. Defaults to -Inf (includes all).
%
% Thresholds(3) defines the scaling factor (eccentricity) to use. 
%   For simulations, you might want to simply use aperture space so the 
%   maximum eccentricity could be set to 1, which is the default. However,
%   you can also use other stimulus designs using this parameter.
%
% Nreps contains the number of times the ground truth was repeated in the
%   simulation which defaults to 1. This only affects how often the
%   ground truths are being replotted so it's not very crucial.
%
% 08/07/2022 - Minor bug fixes & adjusted for SamSrf 9(DSS) 
% 20/07/2022 - Added option to define scaling factor/eccentricity (DSS)
%

if nargin < 2
    Thresholds = [NaN -Inf 1];
end
if length(Thresholds) == 1
    Thresholds = [Thresholds -Inf];
end
if length(Thresholds) == 2
    Thresholds = [Thresholds -Inf 1];
end
if nargin < 3
    Nreps = 1;
end

%% Vectorise data
% Modelled parameters
R2 = Srf.Data(1,:);
mX = Srf.Data(2,:);
mY = Srf.Data(3,:);
mS = Srf.Data(4,:);
mB = Srf.Data(5,:);

% Ground truth parameters
tS = Srf.Ground_Truth(3,:);
nt = size(Srf.Ground_Truth,2) / Nreps;
tX = Srf.Ground_Truth(1,1:nt);
tY = Srf.Ground_Truth(2,1:nt);

% Plot dimensions
dims = max(max(abs(Srf.Ground_Truth(1:2,:))));

%% Remove bad fits
g = R2 > Thresholds(2) | mS <= 0 | mB <= 0;

%% Filter to sigma ground truth
mX = mX(g); 
mY = mY(g); 
mS = mS(g); tS = tS(g);
mB = mB(g);
if isempty(mB)
    samsrf_error(['No data with R^2 > ' num2str(Thresholds(2))]);
end
TitlStr = ['R^2 > ' num2str(Thresholds(2))];
% Limit ground truths?
if ~isnan(Thresholds(1))
    % Only this Sigma 
    g = tS == Thresholds(1);
    mX = mX(g); % X position
    mY = mY(g); % Y position
    mS = mS(g); % Sigma (not used in following plots)
    mB = mB(g); % Beta 
    TitlStr = {TitlStr; ['Ground truth \sigma = ' num2str(Thresholds(1))]};
else
    TitlStr = {TitlStr; 'All ground truth \sigma'};
end
% If this removed all data
if isempty(mB)
    samsrf_error('No data with this ground truth in file!');
end

%% Map density 
w = (3*dims)/200; % Window size
[gx,gy] = meshgrid(-1.5*dims:w:1.5*dims, -1.5*dims:w:1.5*dims);
% Determine density
D = NaN(201,201);
for ix = 1:201
    for iy = 1:201
        D(ix,iy) = sum(mX>=gx(ix,iy)-w & mX<gx(ix,iy)+w & mY>=gy(ix,iy)-w & mY<gy(ix,iy)+w);
    end
end
% Normalize
D = D / prctile(D(:), 99);
D(D > 1) = 1;

%% Plot density 
figure; hold on
contourf(gx, gy, D, 100, 'linestyle', 'none');

% Create colour scheme
Cmap = [1 1 1; ... % White at zero
        linspace(.99,0,63)' linspace(.99,0,63)' linspace(.99,0,63)'; ... % Grey > Black
        linspace(0,0,64)' linspace(0,0,64)' linspace(0,.5,64)'; ... % Black > Blue
        linspace(0,0,64)' linspace(0,1,64)' linspace(.5,0,64)'; ... % Blue > Green
        linspace(0,1,64)' linspace(1,1,64)' linspace(0,0,64)']; % Green > Yellow
colormap(Cmap);

%% Plot ground truths
scatter(tX, tY, 50, 'r', 'filled');

%% Plot mask
set(gca, 'color', [1 1 1]*.7);
[sx,sy] = pol2cart((0:360)/180*pi, Thresholds(3));
h = fill(sx, sy, [1 1 1]*.9);
alpha(h, 0.1);

%% Cosmetic changes
axis([-1.5 1.5 -1.5 1.5]*dims);
axis square
colorbar
set(gca, 'fontsize', 20);
xlabel('Horizontal position');
ylabel('Vertical position');
title(TitlStr);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
