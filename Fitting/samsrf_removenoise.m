function cSrf = samsrf_removenoise(Srf, X, globalcovar)
%
% cSrf = samsrf_removenoise(Srf, X, [globalcovar=true])
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%
%  A global regressor is added automatically so do not include one, 
%  unless you turn this off with the optional boolean input globalcovar.
%
% 29/06/2020 - SamSrf 7 version (DSS)
% 24/05/2021 - Improved speed & removed expansion/compression (DSS)
% 01/07/2021 - Removed utterly redundant expansion/compression again (DSS)
% 12/07/2021 - Added stand-by message since parallel progress reports are a pain (DSS)
% 08/11/2021 - Optional input now allows you to turn off automatical global covariate (DSS)
%

if nargin < 3
    globalcovar = true;
end

%% Regression analysis
cSrf = Srf; % Cleaned Srf
Y = Srf.Data; % Time courses
Nv = size(Srf.Data,2); % Number of vertices
X = zscore(X); % Standardise covariates 
if globalcovar
    X = [X ones(size(X,1),1)]; % Add global regressor
end
B = X \ Y; % Regression betas
mY = X * B; % Modelled time series
R = Y - mY; % Residuals

% Loop thru vertices
disp('Regressing out noise...'); 
disp(' Please stand by...');
Data = zeros(size(cSrf.Data));
clear Srf
parfor v = 1:Nv
    % Only do this for good data
    if ~isnan(sum(Y(:,v)))
        Data(:,v) = R(:,v);
    end
end
cSrf.Data = Data;
