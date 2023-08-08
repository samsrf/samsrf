function Rmap = samsrf_revcor_rmap(Srf, Model, Roi)
%
% Rmap = samsrf_revcor_rmap(Srf, Model, Roi)
%
% Returns the reverse correlation profiles for the vertices in vector Roi
% for the time series data in the reverse correlation pRF map in Srf & 
% the corresponding analysis parameters in Model.
%
% Each row of Rmap is a pixel in the reverse correlation profile.
% Columns correspond to the vertices in Roi.
%
% 08/08/2023 - Written (DSS)
%

% Initialise output R matrix
Rmap = zeros(Model.Rdim^2, length(Roi));
% Dimensions of profiles
dim = sqrt(size(Srf.Regs,2));

% Loop thru ROI vertices
for v = 1:length(Roi)
    % Calculate r-map
    cY = Srf.Y(:,Roi(v));  % Time course of current vertex
    warning off
    cM = [cY ones(size(cY,1),1)] \ Srf.Regs; % Linear regression
    warning on
    cM = cM(1,:); % Remove intercept beta
    % Determine peak
    mM = max(cM); % Find peak activation in each map
    mM = mM(1); % Ensure only one value
    gp = cM > mM/2; % Full area at half maximum
    m = find(cM==mM,1); % Find peak coordinate
    mR = corr(cY, Srf.Regs(:,m)); % Peak correlation with stimulus design
    % Create profile & store parameters
    if ~isempty(mR) && ~isnan(mR) 
        cM = reshape(cM,dim,dim); % Reshape into a map
        if dim ~= Model.Rdim
            cM = imresize(cM,[Model.Rdim Model.Rdim]); % Down-sample r-map
        end
        cM = cM(:); % Vectorise again
        % Store pRF profile
        Rmap(:,v) = cM; % Activation map as vector 
    end
end
