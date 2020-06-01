function R = samsrf_cfcorr(Y, X, S, TR, Hrf)
%
% R = samsrf_cfcorr(Y, X, S, [TR=1, Hrf=[]])
%
% Plots the correlations of time course Y with the regressors in X just as
%   the coarse fit does. S contains the search space parameters. 
% You can retrieve X and S from the search space file src_*.mat.
%
% TR is optional & defines the repetition time (TR) of your scans (default = 1).
% Hrf is optional & contains the HRF vector. Use empty for canonical (default).
%
% Returns the Pearson correlation of Y with each search grid position.
% Also plots the correlation for each X and Y position (i.e. S(1:2,:) using 
%   a separate plot for each unique Sigma (S(3,:). Please note that additional 
%   parameters are currently ignored & thus included in these plots...
%   The plot containing the maximum correlation is indicated by asterisks
%   and the peak correlation is shown by a green cross.
%
% 02/06/2020 - SamSrf 7 version (DSS) 
%

if nargin < 4
    TR = 1;
end
if nargin < 5
    Hrf = [];
end

% Canonical HRF?
if isempty(Hrf)
    Hrf = samsrf_hrf(TR);
end

% Convolve with HRF
for p = 1:size(X,2)
    X(:,p) = prf_convolve_hrf(X(:,p), Hrf);
end

% Calculate correlation
R = corr(Y,X);
% Maximum correlation
M = S(:, R == max(R));

% Plot correlations
U = unique(S(3,:)); % Unique Sigma values
figure('Units', 'Normalized', 'Position', [0 0 1 1]);
p = 0; 
for s = U
    p = p + 1;
    subplot(ceil(sqrt(length(U))), ceil(sqrt(length(U))), p); 
    scatter(S(1,S(3,:)==s)', S(2,S(3,:)==s)', 100, R(S(3,:)==s), 'filled'); 
    hold on
    colormap hotcold
    TitlStr = ['\sigma = ' num2str(s)];
    if M(3) == s
        scatter(M(1), M(2), 200, [0 .5 0], '+', 'linewidth', 2);
        TitlStr = ['*** ' TitlStr ' ***'];
    end
    axis([min(S(1,:)) max(S(1,:)) min(S(2,:)) max(S(2,:))]);
    axis square
    set(gca, 'Clim', [-1 1]); 
    title(TitlStr); 
    colorbar
end