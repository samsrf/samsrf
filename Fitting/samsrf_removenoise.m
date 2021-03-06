function cSrf = samsrf_removenoise(Srf, X)
%
% cSrf = samsrf_removenoise(Srf, X)
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%  A global regressor is added automatically so do not include one.
%
% 29/06/2020 - SamSrf 7 version (DSS)
% 24/05/2021 - Improved speed & removed expansion/compression (DSS)
% 01/07/2021 - Removed utterly redundant expansion/compression again (DSS)
%

%% Regression analysis
cSrf = Srf; % Cleaned Srf
Y = Srf.Data; % Time courses
Nv = size(Srf.Data,2); % Number of vertices
X = zscore(X); % Standardise covariates 
X = [X ones(size(X,1),1)]; % Add global regressor
B = X \ Y; % Regression betas
mY = X * B; % Modelled time series
R = Y - mY; % Residuals

% Loop thru vertices
disp('Regressing out noise...'); 
Data = zeros(size(cSrf.Data));
clear Srf
parfor v = 1:Nv
    % Only do this for good data
    if ~isnan(sum(Y(:,v)))
        Data(:,v) = R(:,v);
    end
end
cSrf.Data = Data;
