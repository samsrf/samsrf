function h = samsrf_simvsfit(Srf, Thresholds, SearchSpace, Value)
%
% h = samsrf_simvsfit(Srf, [Thresholds=[NaN -Inf 1], SearchSpace=[], Value='R^2'])
%
% Plots a comparison of simulated ground truth pRFs & model fits in Srf.
% At present this function assumes that Srf.Data(2:3,:) contain X0 & Y0.
% When you have noisy simulations, the function samsrf_simvsfithist may
% produce more useful results, but your mileage may vary. The output h is
% the handle for the scatter plot.
%
% Srf contains the model fit of a simulated pRF data set, so this must 
%   contain a Srf.Ground_Truth field with true parameters.  
%
% Plots the position shifts as a quiver graph: the dot symbols denote 
%   the modelled positions and the lines denote the shifts from the ground 
%   truth positions. The dot colours can denote any value of your choosing 
%   but this defaults to R^2.
%
%   A circular mask for where the apertures should be is also included - 
%   so if for some reason your stimulus was different (e.g. square?) 
%   you will need to change this to match your design.
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
% Thresholds(3) defines the scaling factor (eccentricity) to use. 
%   For simulations, you might want to simply use aperture space so the 
%   maximum eccentricity could be set to 1, which is the default. However,
%   you can also use other stimulus designs using this parameter.
%
% SearchSpace defines the position parameters of the search grid. 
%   These are plotted on top of the plot as crosses. This must be a 
%   k-by-m matrix where k is the number of parameters in the search grid 
%   and m the number of grid positions. Note, however, only the unique 
%   spatial positions are plotted so it ignores Sigmas, Thetas, or any 
%   other such parameters. You would usually find this search space matrix 
%   in the S variable of your src_*.mat file.
%
% Value is a string that determines what dot colours represent. This can be 
%   any of the strings in Srf.Values. If prefixed with a minus it plots the
%   difference between the modelled & ground truth value (e.g. '-Sigma'), 
%   assuming the corresponding row in Ground_Truth (e.g. for Sigma this is
%   usually Srf.Data(4,:)-Srf.Ground_Truth(3,:)). You can also use '%'
%   instead of '-' to plot the percentage difference. Obviously, this only
%   works for data rows which have a ground truth!
%
% 08/07/2022 - Minor bug fixes & adjusted for SamSrf 9(DSS) 
% 20/07/2022 - Added option to define scaling factor/eccentricity (DSS)
% 16/11/2023 - Dot colour can now denote any value you choose (DSS)
%              Removed sigma plot when setting Thresholds(1) to NaN (DSS)
%              Reduced the axis limits of the plot (DSS)
%              Colour schemes for zero-bounded values now start at 0 (DSS)
% 18/11/2023 - Dot colour can now also plot differences & percentages (DSS)
% 20/11/2023 - Now returns handle to the scatter dots (DSS)
%

if nargin < 2
    Thresholds = [NaN -Inf 1];
end
if length(Thresholds) == 1
    Thresholds = [Thresholds -Inf];
end
if length(Thresholds) == 2
    Thresholds = [Thresholds 1];
end
if nargin < 3
    SearchSpace = [];
end
if nargin < 4
    Value = 'R^2';
end

%% Calculate differences?
if Value(1) == '-' 
    Value = Value(2:end);
    CalcDiff = '-';
elseif Value(1) == '%'
    Value = Value(2:end);
    CalcDiff = '%';
else
    CalcDiff = ' ';
end

%% Retrieve data
% Modelled parameters
R2 = Srf.Data(1,:);
mX = Srf.Data(2,:);
mY = Srf.Data(3,:);
mS = Srf.Data(4,:);
r = find(strcmpi(Srf.Values, Value)); % Which row are these data?
mV = Srf.Data(r,:);
% Calculate difference?
if CalcDiff ~= ' ' && r > 1 && r <= size(Srf.Ground_Truth,1) + 1
    if CalcDiff == '-'
        mV = mV - Srf.Ground_Truth(r-1,:); % Difference from ground truth
    elseif CalcDiff == '%'
        mV = (mV - Srf.Ground_Truth(r-1,:)) ./ Srf.Ground_Truth(r-1,:) * 100; % Percentage difference from ground truth
    end
end
% If value is no good
if isempty(mV)
    samsrf_error('Invalid value specified!');
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
mV = mV(g);
if isempty(mV)
    samsrf_error(['No data with R^2 > ' num2str(Thresholds(2))]);
end
TitlStr = ['R^2 > ' num2str(Thresholds(2))];

%% Limit ground truths?
if ~isnan(Thresholds(1))
    % Only this Sigma 
    g = tS == Thresholds(1);
    mX = mX(g); tX = tX(g); % X position
    mY = mY(g); tY = tY(g); % Y position
    mS = mS(g); tS = tS(g); % Sigma (not used in following plots)
    mV = mV(g); % Value 
    TitlStr = {TitlStr; ['Ground truth \sigma = ' num2str(Thresholds(1))]};
else
    TitlStr = {TitlStr; 'All ground truth \sigma'};
end
% If this removed all data
if isempty(mV)
    samsrf_error('No data with this ground truth in file!');
end

%% Scales for quiver plot
sB = max(abs([min(mV) max(mV)])); % Value difference / R^2 scale
if sB == 0
    sB = 1;
end

%% Plot mask
figure
hold on
set(gca, 'color', [1 1 1]*.7);
[sx,sy] = pol2cart((0:360)/180*pi, Thresholds(3));
h = fill(sx, sy, [1 1 1]*.9);
alpha(h, 0.3);

%% Plot lines
for v = 1:length(tX)
    line([tX(v) mX(v)], [tY(v) mY(v)], 'color', [1 1 1]/2);
end

%% Plot modelled positions 
h = scatter(mX, mY, 60, mV, 'filled', 'markeredgecolor', 'k');

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
axis square
if CalcDiff == ' ' && (strcmpi(Value, 'R^2') || strcmpi(Value, 'nR^2') || strcmpi(Value, 'Sigma') || strcmpi(Value, 'Sigma1') || strcmpi(Value, 'Sigma2') || strcmpi(Value, 'Fwhm') || strcmpi(Value, 'Centre') || strcmpi(Value, 'Surround'))  
    set(gca, 'fontsize', 20, 'Clim', [0 1]*sB);
    if strcmpi(Value, 'R^2') || strcmpi(Value, 'nR^2') 
        colormap batlow
    else
        colormap hawaii
    end
else
    set(gca, 'fontsize', 20, 'Clim', [-1 1]*sB);
    colormap berlin
end
xlabel('Horizontal position');
ylabel('Vertical position');
title(TitlStr);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
if CalcDiff == '-'
    CbStr = ['\Delta ' Srf.Values{strcmpi(Srf.Values, Value)}];
elseif CalcDiff == '%'
    CbStr = ['\Delta ' Srf.Values{strcmpi(Srf.Values, Value)} '%'];
else
    CbStr = ['Model ' Srf.Values{strcmpi(Srf.Values, Value)}];
end
set(get(cb, 'Label'), 'String', CbStr);
axis([-1 1 -1 1]*dims*1.05);
