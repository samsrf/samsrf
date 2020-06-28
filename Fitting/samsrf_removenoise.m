function cSrf = samsrf_removenoise(Srf, X)
%
% cSrf = samsrf_removenoise(Srf, X)
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%  A global regressor is added automatically so do not include one.
%
% 29/06/2020 - SamSrf 7 version (DSS)
%

%% Expand Srf if necessary
[Srf,vx] = samsrf_expand_srf(Srf);

%% Regression analysis
cSrf = Srf; % Cleaned Srf
Y = Srf.Data; % Time courses
Nv = size(Srf.Data,2); % Number of vertices
X = zscore(X); % Standardise covariates 
X = [X ones(size(X,1),1)]; % Add global regressor

% Loop thru vertices
disp('Regressing out noise...'); 
Data = zeros(size(cSrf.Data));
parfor v = 1:Nv
    if ~isnan(sum(Y(:,v)))
        B = regress(Y(:,v), X);
        nY = X * B;
        Data(:,v) = Y(:,v) - nY;
    end
end
cSrf.Data = Data;

%% Compress Srf again if needed
cSrf = samsrf_compress_srf(cSrf,vx);
