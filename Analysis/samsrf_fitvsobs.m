function [S, X, Y] = samsrf_fitvsobs(Srf, Model, v)
%
% [S, X, Y] = samsrf_fitvsobs(Srf, Model, v)
%
% Plots the predicted (X) and observed (Y) time course in pRF-fit Srf at vertex v.
% You must define Model to convolve the predicted time course with the HRF.
% Also reports the model parameters (S) in the figure.
%
% 12/02/2022 - Now allows plotting convolved predictions stored in coarse-fit data files (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 12/07/2023 - Fixed bug with downsampling predictions (DSS)
%

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Observed time course
Y = Srf.Y(:,v); 

% Predicted time course
if isfield(Srf, 'X_glm')
    % GLM design matrix in Srf
    Model.Downsample_Predictions = 1; % Dummy model parameter
    TR = 1;  % Dummy TR value but will display volumes 
    
    % Create full GLM predictor
    X = Srf.X_glm * Srf.Data(1:size(Srf.X_glm,2), v);
elseif isfield(Srf, 'X')
    % Predictions defined in Srf
    X = Srf.X(:,v); % Predicted time course
    % Backwards compatibility with versions prior to 7.5
    if ~isfield(Model, 'Downsample_Predictions')
        Model.Downsample_Predictions = 1; % In case downsampling undefined
    end
    % Convolve prediction with HRF?
    if isfield(Model, 'Hrf') && ~Model.Coarse_Fit_Only % Only needed for fine-fits
        X = prf_convolve_hrf(X, Model.Hrf, 1); % No downsampling to retain temporal resolution of prediction!
    end

    % If raw data exists use that
    if isfield(Srf, 'Raw_Data')
        S = Srf.Raw_Data(:,v);
    else
        S = Srf.Data(:,v);
    end

    % Scale predictor appropriately
    if sum(contains(Srf.Values, 'Beta')) && Srf.Version >= 6
        % Multiply by Beta
        br = strcmpi(Srf.Values,'Beta'); % Which value contains beta?
        cr = strcmpi(Srf.Values,'Baseline'); % Which value contains baseline?
        X = X * S(br) + S(cr); % Scale predictor
    else
        % Z-transform time series
        X = zscore(X);
        Y = zscore(Y);
    end
end

% Plot time courses
if isfield(Model, 'TR')
    TR = Model.TR;
else
    TR = 1;
end
hold off
plot((1:length(Y))*TR*Model.Downsample_Predictions, Y, 'color', [.5 .5 1], 'linewidth', 2);
hold on
plot((1:length(X))*TR, X, 'r', 'linewidth', 2);
xlim([1 length(Y)]*TR);
grid on
line(xlim, [0 0], 'color', [1 1 1]/2, 'linewidth', 2);
legend({'Observed' 'Predicted'});
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [.1 .1 .8 .4]);
set(gca, 'fontsize', 12);
if isfield(Model, 'TR')
    xlabel('Time (s)');
else
    xlabel('Volumes (#)');
end
ylabel('Response (z)');
ts = '';
if exist('S', 'var')
    for i = 1:size(S,1)
        ts = [ts Srf.Values{i} '=' num2str(round(S(i),2)) '   '];
    end
    title({['Vertex: ' num2str(v)]; ts});
else
    title(['Vertex: ' num2str(v) ', R^2=' num2str(corr(X,Y)^2)]);
end
