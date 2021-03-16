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
%

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Time courses
X = Srf.X(:,v); % Predicted time course
Y = Srf.Y(:,v); % Observed time course
% Convolve prediction with HRF
X = prf_convolve_hrf(X, Model.Hrf);

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

if find(strcmpi(Srf.Values, 'Mu'))
    [Ax, line1, line2] = plotyy(1:length(X), X, 1:length(Y), Y);
    line1.Color = [1 1 1]/2;
    line1.LineWidth = 2;
    line2.Color = [0 0 1];
    line2.LineWidth = 2;
    Ax(2).FontSize = 15;
    Ax(2).YColor = [1 1 1]/2;
    ylabel(Ax(2), 'Observed');
else    
    plot(Y, 'color', [1 1 1]/2, 'linewidth', 2);
    hold on
    plot(X, 'b', 'linewidth', 2);
end

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
