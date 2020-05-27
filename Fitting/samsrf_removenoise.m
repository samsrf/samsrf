function cSrf = samsrf_removenoise(Srf, X)
%
% cSrf = samsrf_removenoise(Srf, X)
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%  A global regressor is added automatically so do not include one.
%
% 09/08/2018 - SamSrf 6 version (DSS)
% 27/05/2020 - Streamlined how waitbar is handled (DSS)
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
h = samsrf_waitbar('Regressing...'); 
for v = 1:Nv
    if ~isnan(sum(Y(:,v)))
        B = regress(Y(:,v), X);
        nY = X * B;
        cSrf.Data(:,v) = Y(:,v) - nY;
    end
    samsrf_waitbar(v/Nv, h); 
end
samsrf_waitbar('', h); 

%% Compress Srf again if needed
cSrf = samsrf_compress_srf(cSrf,vx);
