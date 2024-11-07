function R = samsrf_gsr2(Y, X, S, Model) 
%
% R = samsrf_gsr2(Y, X, S, [Model])
%
% Plots the squared correlations of time course Y with the regressors in X 
%   just as the grid search (coarse fit) does. S contains the search space 
%   parameters. You can retrieve X and S from the search space file src_*.mat.
%
% The third input Model is optional. This can be the model structure for
%   your pRF analysis and it defines the TR, Hrf, Downsampling constant,
%   and also defines whether the model is 2D or 1D, and if you used a polar
%   or a Cartesian search grid.
%
%   Model.TR is defines the repetition time (TR) of your scans.
%       This can also be the temporal resolution of the stimulus if that is faster.
%       In that case, use Downsampling to match with the true scanner TR.
%
%   Model.Hrf contains the HRF vector. Use empty for canonical or 1 for none.
%
%   Model.Downsampling defines the factor by which to match stimulus timing with true TR (default = 1).
%
%   Model.Polar_Search_Space toggles whether a polar or Cartesian search space is used (default = false).
%
%   Model.Param_Names is a cell array defining the names of free parameters in your pRF model
%       (excluding betas). This is used to determine if the model is 2D or 1D.
%
%   Look at the Model structure for forward-model pRFs for more information.
%
% Returns the R^2 of Y with each search grid position.
%
% Also plots the R^2 for each X and Y position (i.e. S(1:2,:) using a separate plot for each unique Sigma (S(3,:). 
%   If this is a 1D model, there is only one plot comparing the tuning (X-axis) to the sigma parameter (Y-axis).
% IMPORTANT: Any additional parameters are currently ignored & thus included in these plots...
%
% The plot containing the maximum correlation is indicated by asterisks
%   and the peak correlation is shown by a green cross.
%
% 02/03/2022 - Added support for polar search grids & 1D models (DSS) 
% 14/04/2022 - Renamed function to avoid confusion with CF analysis (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 03/12/2023 - Bugfix when using concurrent HRF fitting (DSS)
% 04/11/2024 - Fixed bug assuming wrong HRF when concurrent fitting or using SPM (DSS) 
%

if ~isfield(Model, 'Downsampling')
    Model.Downsampling = 1;
end
if ~isfield(Model, 'Polar_Search_Space')
    Model.Polar_Search_Space = false;
end

TR = Model.TR; % Repetition time (TR)
Hrf = Model.Hrf; % Haemodynamic response
Downsampling = Model.Downsampling; % Microtime resolution when downsampling
IsPolar = Model.Polar_Search_Space; % Are we using polar search grid (2D models only)?
if length(Model.Param_Names) == 2    
    Is1D = true; % This is a 1D model
    IsPolar = false; % No polar grid possible so ignore this
else
    Is1D = false; % This is a 2D model
end

% Canonical HRF or concurrent fitting?
if isempty(Hrf) 
    Hrf = samsrf_hrf(TR); % de Haas canonical HRF
elseif isinf(Hrf) || Hrf == 0
    Hrf = samsrf_doublegamma(TR); % SPM 12 canonical HRF
end

% Convolve with HRF
for p = 1:size(X,2)
    X(:,p) = prf_convolve_hrf(X(:,p), Hrf, Downsampling); % Convolve & downsample if desired
end

% Calculate R^2
R = corr(Y,X).^2;
% Maximum correlation
M = S(:, R == max(R));

% Search space parameters
xU = unique(S(1,:)); % Unique X values
yU = unique(S(2,:)); % Unique Y values
sU = unique(S(3,:)); % Unique Sigma values

% Polar or Cartesian search grid?
if IsPolar
    % Polar search grid is same as unique 
    Xc = S(1,:);
    Yc = S(2,:);
else
    [Xc,Yc] = meshgrid(xU,yU); % Meshgrid of unique X & Y values
end

% Plot correlations
figure('Units', 'Normalized', 'Position', [0 0 1 1]);
p = 0; 
for s = sU
    p = p + 1;
    subplot(ceil(sqrt(length(sU))), ceil(sqrt(length(sU))), p); 
    hold on
    if IsPolar
        polarpatch(Xc(S(3,:)==s), Yc(S(3,:)==s), R(S(3,:)==s));
    else
        Rc = NaN(length(yU), length(xU));
        for x = 1:length(xU)
            for y = 1:length(yU)
                Rc(y,x) = R(S(1,:)==Xc(y,x) & S(2,:)==Yc(y,x) & S(3,:)==s); % Store correlation
            end
        end
        contourf(Xc, Yc, Rc, 40, 'edgecolor', 'none'); 
        axis([min(S(1,:)) max(S(1,:)) min(S(2,:)) max(S(2,:))]);
        axis square
    end
    TitlStr = ['\sigma = ' num2str(s)];
    if M(3) == s
        scatter(M(1), M(2), 200, [0 .5 0], '+', 'linewidth', 2);
        TitlStr = ['*** ' TitlStr ' ***'];
    end
    colormap berlin
    caxis([-1 +1]*max(R)); 
    if Is1D
        title('1D pTC model');
        xlabel(Model.Param_Names{1});
        ylabel(Model.Param_Names{2});
    else
        title(TitlStr); 
    end
    colorbar
end