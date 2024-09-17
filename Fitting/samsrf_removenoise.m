function cSrf = samsrf_removenoise(Srf, X, mver, globalcovar)
%
% cSrf = samsrf_removenoise(Srf, X, [mver=[], globalcovar=true])
%
% Uses linear regression to remove the noise (covariates of no interest in X)
%  from the time series data in Srf.Data. Returns the "cleaned" Srf.
%
%  By default, this analysis is performed on all vertices in Srf. 
%  You can limit this to a ROI by providing a list of vertices in mver.
%
%  A global regressor is added automatically so do not include one, 
%  unless you turn this off with the optional boolean input globalcovar.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 15/09/2022 - Analysis can now be restricted to a ROI mask (DSS)
%

if nargin < 3
    mver = [];
end
if isempty(mver)
    mver = 1:size(Srf.Vertices,1);
end
if nargin < 4
    globalcovar = true;
end

%% Regression analysis
samsrf_disp('Regressing out noise...'); 
cSrf = Srf; % Cleaned Srf
Y = Srf.Data(:,mver); % Time courses
Nv = length(mver); % Number of vertices
X = zscore(X); % Standardise covariates 
if globalcovar
    X = [X ones(size(X,1),1)]; % Add global regressor
end
B = X \ Y; % Regression betas
mY = X * B; % Modelled time series
R = Y - mY; % Residuals

% Loop thru vertices
samsrf_disp(' Removing covariates...');
Data = zeros(size(cSrf.Data,1), Nv);
clear Srf
parfor v = 1:Nv
    % Only do this for good data
    if ~isnan(sum(Y(:,v)))
        Data(:,v) = R(:,v);
    end
end
cSrf.Data(:,mver) = Data;
