function samsrf_simvsfit(Srf, Value, R2Thr, ParVal, SrcSpc)
%
% samsrf_simvsfit(Srf, [Value='Sigma', R2Thr=0, ParVal=Inf, SrcSpc=[]])
%
% Plots a comparison of simulated ground truth pRFs & model fits.
%
% Srf contains the model fit of a simulated pRF data set, so this must 
% contain a Srf.Ground_Truth field. It assumes that Srf.Data(2:3,:) 
% contains the X and Y coordinates & plots the differences as a quiver plot.
% The dot symbols denote the modelled positions and the lines denote the
% shifts from the ground truth positions.
% 
% Value defines which field of non-position parameters is being compared.
% This is denoted by the colours of the dots in the scatter plot.
% If this is set to R^2 or Beta it plots the raw value because there is no
% ground truth for these parameters. For other parameters it plots the
% difference between the ground truth and the modelled parameter estimate.
% This is optional and defaults to 'Sigma'.
%
% R2Thr defines the R^2 threshold of the model fits to include in the 
% comparison. This is optional and defaults to zero.
%
% ParVal limits the plotted data to only those ground truth values (as
% defined by Value) that fall within the range defined by this vector. 
% A scalar limits it to only that single value. This is optional and
% defaults to Inf which doesn't limit the range at all. (Note that obviously
% this doesn't work for R^2, Beta, or Baseline).
%
% SrcSpc defines the position parameters of the search grid. These are
% plotted as crosses. It assumes a k-by-m matrix where n is the number
% of parameters in the search grid and m the number of grid position.
% However, only the unique spatial positions are plotted so it ignores
% sigmas, thetas, or any other such parameters. You would usually find this
% search space matrix in the S variable of your src_*.mat file.
%
% 28/05/2020 - Written (DSS)
% 29/05/2020 - Added option to plot search grid (DSS)
%

if nargin < 2
    Value = 'Sigma';
end
if nargin < 3
    R2Thr = 0;
end
if nargin < 4
    ParVal = Inf;
end
if nargin < 5
    SrcSpc = [];
end

%% Retrieve data
% Modelled positions
mX = Srf.Data(2,:);
mY = Srf.Data(3,:);

% Ground truth posiitonsti
tX = Srf.Ground_Truth(1,:);
tY = Srf.Ground_Truth(2,:);

% Extract scatter data
r = find(strcmpi(Srf.Values, Value));
mD = Srf.Data(r,:); 
if strcmpi(Value, 'R^2') || strcmpi(Value, 'Beta') || strcmpi(Value, 'Baseline')
    tD = mD; % No ground truth exists 
    Delta = mD; % Used for scaling
else
    tD = Srf.Ground_Truth(r-1,:);
    Value = ['\Delta_{' Value '}'];
    Delta = mD-tD; % Difference between modelled & true value
end
% Scale for plot
s = max(abs([min(Delta) max(Delta)])); 

%% Filter data
% Remove bad fits
g = Srf.Data(1,:) > R2Thr;
mX = mX(g); tX = tX(g);
mY = mY(g); tY = tY(g);
mD = mD(g); tD = tD(g);
Delta = Delta(g); 
% Limit ground truths?
LabStr = Value;
if ~isinf(ParVal)
    if length(ParVal) == 1
        ParVal = [1 1]*ParVal;
    end
    g = tD >= ParVal(1) & tD <= ParVal(2);
    mX = mX(g); tX = tX(g);
    mY = mY(g); tY = tY(g);
    mD = mD(g); tD = tD(g);
    Delta = Delta(g); 
    Value = [num2str(ParVal(1)) ' \leq ' Value ' \geq ' num2str(ParVal(2))]; % Add range if defined
end
Value = {Value; ['R^2 > ' num2str(R2Thr)]};

%% Plot mask
hold on
set(gca, 'color', [1 1 1]*.7);
[mx,my] = pol2cart((0:360)/180*pi, 1);
fill(mx, my, [1 1 1]*.9);

%% Plot lines
for v = 1:length(tX)
    line([tX(v) mX(v)], [tY(v) mY(v)], 'color', 'k');
end

%% Plot modelled positions 
scatter(mX, mY, 100, Delta, 'filled', 'markeredgecolor', 'k');

%% Plot search space
if ~isempty(SrcSpc)
    s3 = unique(SrcSpc(3,:));
    sX = SrcSpc(1,SrcSpc(3,:) == s3(1));
    sY = SrcSpc(2,SrcSpc(3,:) == s3(1));
    scatter(sX, sY, 100, 'k+');
end

%% Cosmetic changes
cb = colorbar;
colormap hotcold
axis([-1 1 -1 1]*2.5);
axis square
set(gca, 'fontsize', 20, 'Clim', [-1 1]*s);
xlabel('Horizontal position');
ylabel('Vertical position');
title(Value);
set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
set(get(cb, 'Label'), 'String', LabStr);