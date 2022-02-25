function R = samsrf_cfcorr(Y, X, S, TR, Hrf, Downsampling)
%
% R = samsrf_cfcorr(Y, X, S, [TR=1, Hrf=[], Downsampling=1])
%
% Plots the squared correlations of time course Y with the regressors in X 
%   just as the coarse fit does. S contains the search space parameters. 
% You can retrieve X and S from the search space file src_*.mat.
%
% TR is optional & defines the repetition time (TR) of your scans (default = 1).
%   This can also be the temporal resolution of the stimulus if that is faster.
%   In that case, use Downsampling to match with the true scanner TR.
% Hrf is optional & contains the HRF vector. Use empty for canonical (default).
% Downsampling is optional & defines the factor by which to match stimulus timing with true TR (default = 1).
%
% Returns the R^2 of Y with each search grid position.
% Also plots the R^2 for each X and Y position (i.e. S(1:2,:) using 
%   a separate plot for each unique Sigma (S(3,:). Please note that additional 
%   parameters are currently ignored & thus included in these plots...
%   The plot containing the maximum correlation is indicated by asterisks
%   and the peak correlation is shown by a green cross.
%
% 28/06/2020 - SamSrf 7 version (DSS) 
% 26/02/2022 - Fixed bug with HRF convolution not being used (DSS) 
%

if nargin < 4
    TR = 1;
end
if nargin < 5
    Hrf = [];
end
if nargin < 6
    Downsampling = 1;
end

% Canonical HRF?
if isempty(Hrf)
    Hrf = samsrf_hrf(TR);
end

% Convolve with HRF
for p = 1:size(X,2)
    X(:,p) = prf_convolve_hrf(X(:,p), Hrf, Downsampling); % Convolve & downsample if desired
end

% Calculate R^2
R = corr(Y,X).^2;
% Maximum correlation
M = S(:, R == max(R));

% Plot correlations
xU = unique(S(1,:)); % Unique X values
yU = unique(S(2,:)); % Unique Y values
sU = unique(S(3,:)); % Unique Sigma values
[Xc,Yc] = meshgrid(xU,yU); % Meshgrid of unique X & Y values
figure('Units', 'Normalized', 'Position', [0 0 1 1]);
p = 0; 
for s = sU
    p = p + 1;
    subplot(ceil(sqrt(length(sU))), ceil(sqrt(length(sU))), p); 
    Rc = NaN(length(xU), length(yU));
    for x = 1:length(xU)
        for y = 1:length(yU)
            Rc(x,y) = R(S(1,:)==Xc(x,y) & S(2,:)==Yc(x,y) & S(3,:)==s);
        end
    end
    contourf(Xc, Yc, Rc, 40, 'edgecolor', 'none'); 
    hold on
    colormap hot
    TitlStr = ['\sigma = ' num2str(s)];
    if M(3) == s
        scatter(M(1), M(2), 200, [0 .5 0], '+', 'linewidth', 2);
        TitlStr = ['*** ' TitlStr ' ***'];
    end
    axis([min(S(1,:)) max(S(1,:)) min(S(2,:)) max(S(2,:))]);
    axis square
    set(gca, 'Clim', [0 1]*max(R)); 
    title(TitlStr); 
    colorbar
end