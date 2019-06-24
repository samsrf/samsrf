function cSrf = samsrf_removenoise(Srf, X)
%
% cSrf = samsrf_removenoise(Srf, X)
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%  A global regressor is added automatically so do not include one.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

%% Check if waitbar is to be used
wb = samsrf_waitbarstatus;

%% Expand Srf if necessary
[Srf,vx] = samsrf_expand_srf(Srf);

%% Regression analysis
cSrf = Srf; % Cleaned Srf
Y = Srf.Data; % Time courses
Nv = size(Srf.Data,2); % Number of vertices
X = zscore(X); % Standardise covariates 
X = [X ones(size(X,1),1)]; % Add global regressor

% Loop thru vertices
if wb h = waitbar(0, 'Regressing...', 'Units', 'pixels', 'Position', [100 100 360 70]); end
for v = 1:Nv
    if ~isnan(sum(Y(:,v)))
        B = regress(Y(:,v), X);
        nY = X * B;
        cSrf.Data(:,v) = Y(:,v) - nY;
    end
    if wb waitbar(v/Nv, h); end
end
if wb close(h); end

%% Compress Srf again if needed
cSrf = samsrf_compress_srf(cSrf,vx);
