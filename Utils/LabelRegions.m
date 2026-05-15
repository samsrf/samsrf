function Lab = LabelRegions(Srf, Tmp, Rad, Ecc)
%
% Lab = LabelRegions(Srf, Tmp, [Rad=10, Ecc=[0 Inf]])
%
% Simple decoder that uses a set of spatially normalised maps warped into
% the subject's native space to determine the most likely ROI each vertex
% in a retinotopic pRF map belongs to. Saves the labels with the suffix '_pred'.
%
% Goes through all vertices in the template maps where at least one individual 
% has been labelled & creates a Euclidean searchlight (in sphere space) around it. 
% It uses the mean correlation for both the x0 and y0 maps in the searchlight 
% between the empirical data and the comparison subjects to calculate weighted mean 
% of ROI labels to predict the ROI (rounded to categorise).
%
%   Srf:    Empirical pRF map (assumes Srf.Data(2:3) are x0 and y0 coordinates)
%   Tmp:    Comparison maps warped from template into native space (using Template2NativeMaps)
%   Rad:    Radius in mm of search light around each vertex to use (default = 10)
%   Ecc:    Eccentricity range of empirical map to be considered (default = [0 Inf])
%
% Also returns the predicted label for each vertex in Lab(1,:).
% The second row Lab(2,:) contains the proportion of consistent labels across the training data.
% You can use this as a confidence probability for the label.
%
% 04/05/2026 - Written (DSS)
%

if nargin < 3
    Rad = 10; % Searchlight radius
end
if nargin < 4
    Ecc = [0 Inf]; % Eccentricity range
end

%% Mask out eccentricity range
Eccen = sqrt(sum(Srf.Data(2:3,:).^2)); % Empirical eccentricity map
Emsk = Eccen < Ecc(1) | Eccen > Ecc(2); % Vertices out of eccentiricy range
Srf.Data(:,Emsk) = 0;
Tmp.Data(:,Emsk) = 0;

%% Loop thru mask vertices
Lab = zeros(2,size(Srf.Vertices,1)); % Labelling empirical map
mver = find(Tmp.Data(6,:) > 0); % Only vertices labelled in any individual map
samsrf_progbar(0);
for i = 1:length(mver)
    % Current vertex info
    v = mver(i);
    L = squeeze(Tmp.Raw_Data(5,v,:)); % ROI labels of this vertex   
    sl = find(sqrt(sum((Srf.Sphere-Srf.Sphere(v,:)).^2,2)) < Rad); % Searchlight around this vertex   
    
    % Closest matching individual map
    Rx = corr(Srf.Data(2,sl)', squeeze(Tmp.Raw_Data(2,sl,:)))'; % Correlations with X-coordinates
    Ry = corr(Srf.Data(3,sl)', squeeze(Tmp.Raw_Data(3,sl,:)))'; % Correlations with Y-coordinates
    W = (Rx+Ry)/2; % Mean correlation to be used as weight
    W(W < 0) = 0; % Negative correlations weighted as zero
    W(L == 0) = 0; % Any individuals were this vertex is unlabelled don't matter    
    Pred = wmode(L, W); % Predicted label for this vertex
    Lab(1,v) = Pred; % Rounded weighted mean of ROI label
    Lab(2,v) = mean(L==Pred); % Proportion of consistent labels

    samsrf_progbar(i/length(mver));
end

%% Save out ROI labels
if Srf.Hemisphere(1) == 'l'
    h = 'lh_';
elseif Srf.Hemisphere(1) == 'r'
    h = 'rh_';
else
    h = '';
end
% Loop thru ROI labels
for r = 1:nanmax(Lab(1,:))
    samsrf_srf2label(Srf, [h Tmp.Roi_Names{r} '_pred'], 1, find(Lab(1,:)==r));
end

%% Weighted mode of labels
function wL = wmode(L, W)
% Returns mode of labels L weighted by W

uL = unique(L); % Unique labels
uL = uL(uL > 0); % Don't even bother including zero
WbL = []; % Summed weights by unique label
% Loop thru unique labels
for u = uL'
    WbL = [WbL; u sum(W(L==u))];
end
wL = WbL(find(WbL(:,2)==max(WbL(:,2)),1),1); % Return label with highest weight
% Ensure not empty
if isempty(wL)
    wL = 0;
end