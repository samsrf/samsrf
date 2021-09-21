function [S, X, Y] = samsrf_fitvsobs(Srf, Model, v)
%
% [S, X, Y] = samsrf_fitvsobs(Srf, Model, v)
%
% Plots the predicted (X) and observed (Y) time course in pRF-fit Srf at vertex v.
% You must define Model to convolve the predicted time course with the HRF.
% Also reports the model parameters (S) in the figure.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 17/03/2021 - Removed redundant rounding function (DSS) 
% 22/09/2021 - Now supports downsampling if stimulus timing mismatches TR (DSS)
%              Removed dual-Y axis for tuning curves as no idea why that was there... (DSS) 
%              Changed predicted timeseries colour to red (DSS)
%

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Time courses
X = Srf.X(:,v); % Predicted time course
Y = Srf.Y(:,v); % Observed time course
% Convolve prediction with HRF
if ~isfield(Model, 'Downsample_Predictions')
    Model.Downsample_Predictions = 1; % Backwards compatibility with versions prior to 7.5
end
X = prf_convolve_hrf(X, Model.Hrf, Model.Downsample_Predictions);

% If raw data exists use that
if isfield(Srf, 'Raw_Data')
    S = Srf.Raw_Data(:,v);
else
    S = Srf.Data(:,v);
end

% Multiply by Beta
br = strcmpi(Srf.Values,'Beta'); % Which value contains beta?
cr = strcmpi(Srf.Values,'Baseline'); % Which value contains baseline?
X = X * S(br) + S(cr); % Scale predictor 

% Plot time courses
figure('name', ['Vertex: ' num2str(v)]);
plot(Y, 'color', [1 1 1]/2, 'linewidth', 2);
hold on
plot(1:Model.Downsample_Predictions:length(Y), X, 'r', 'linewidth', 2);

legend({'Observed' 'Predicted'});
set(gcf, 'Units', 'normalized');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-pos(3) pos(2) pos(3)*3 pos(4)]);
set(gca, 'fontsize', 15);
xlabel('Volumes (#)');
ylabel('Response (z)');
ts = '';
for i = 1:size(S,1)
    ts = [ts Srf.Values{i} '=' num2str(round(S(i),2)) '   '];
end
title(ts);
