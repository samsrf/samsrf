function xysb = samsrf_backproj_sctr(Response, pRF_Data, Threshold)
%
% xysb = samsrf_backproj_sctr(Response, pRF_Data, [Threshold=[0, 85, 0, Inf]])
%
% Projects the activity values in Response back into visual space as a
% scatter plot. pRF size is determined by pRF_Data(4,:).
% Response values are sorted in ascending order so that the strongest
% response is shown on top.
%
% The optional input Threshold defines the minimal R^2 of the pRFs to be projected 
% and the percentile beyond which to plot responses (across all of Response).
% The R^2 threshold defaults to 0 and the percentile defaults to 85.
%
%   The third and fourth value define the inner and outer eccentricity to be used. 
%   This defaults to the most inclusive criteria, i.e. all R^2 and all eccentricities.
%
% Returns a cell array with as many components as there are rows in Response. 
% Each component contains a 4-column matrix which you can use for a scatter plot:
%
%   h = scatter(xysb{1}(:,1), xysb{1}(:,2), xysb{1}(:,3), xysb{1}(:,4), 'filled');
%   scatter_size(h);
%   alpha(h, 0.2);
%   axis square
%   colormap(flipud(gray));
%
% Because of the stupid way Matlab determines the size of symbols in scatter plots,
% you need to manually adjust the sizes with this utility function scatter_size.
% (The same problem also applies to samsrf_polarplot)
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 08/09/2021 - Finally fixed the issue with plotting pRF size properly (DSS)
% 13/10/2021 - Fixed misleading typo in the help section (DSS)
%

if nargin < 3
    Threshold = [];
end
if isempty(Threshold)
    Threshold = [0 85 0 Inf]; 
elseif length(Threshold) == 1
    Threshold = [Threshold 85 0 Inf];
elseif length(Threshold) == 2
    Threshold = [Threshold 0 Inf];
elseif length(Threshold) == 3
    Threshold = [Threshold Inf];
end

% pRF map data
gof = pRF_Data(1,:); % Goodness of fit
ecc = sqrt(pRF_Data(2,:).^2 + pRF_Data(3,:).^2); % Eccentricity
sigma = pRF_Data(4,:).^2; % pRF size squared

% Percentile level
Pct = prctile(Response(:), Threshold(2));

% Filter vertices
nanresp = isnan(Response(1,:)); % Determine NaNs 
Response = Response(:, gof > Threshold(1) & ecc > Threshold(3) & ecc < Threshold(4) & sigma > 0 & nanresp == 0); % Throw out NaNs 
pRF_Data = pRF_Data(:, gof > Threshold(1) & ecc > Threshold(3) & ecc < Threshold(4) & sigma > 0 & nanresp == 0); % Throw out NaNs 

% Loop thru volumes
xysb = {};
for i = 1:size(Response,1)
    % Sort current response & map
    CurrResp = Response(i,:);
    [CurrResp,sx] = sort(CurrResp, 'ascend');
    CurrMap = pRF_Data(:,sx);
    % Only vertices above percentile
    pv = CurrResp > Pct;    
    % Add to output
    xysb{i} = [CurrMap(2,pv)' CurrMap(3,pv)' CurrMap(4,pv)' CurrResp(pv)'];    
end