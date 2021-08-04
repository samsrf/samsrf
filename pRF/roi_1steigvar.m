function [tY, r2] = roi_1steigvar(Ys)
%
% [tY, r2] = roi_1steigvar(Ys)
%
% Returns the first eigenvariate for the time series in Ys (row = time points, columns = individual vertices).
% The second output argument contains the variance explained by this component.
%
% 11/07/2021 - Written (DSS)
% 04/08/2021 - Now returns NaNs when input is empty (DSS)
%

if ~isempty(Ys)
    % Remove redundant time series
    Ys = rem_redcols(Ys);

    % Principal component analysis
    warning off
    [~,tY,~,~,r2] = pca(Ys(:,~isnan(sum(Ys))));
    warning on
    if ~isempty(tY)
        % Only take first component
        tY = tY(:,1); % 1st component
        r2 = r2(1); % Variance explained

        % Z-score time course
        tY = zscore(tY);
    end
else
    tY = NaN;
    r2 = NaN;
end