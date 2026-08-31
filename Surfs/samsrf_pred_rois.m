function Srf = samsrf_pred_rois(SrfOrName, Tmp, Kernels, Threshs)
%
% Srf = samsrf_pred_rois(SrfOrName, Tmp, [Kernels=[10 0], Threshs=[0.1 0.01 0 Inf]])
%
% Simple decoder that uses a set of spatially normalised maps warped into
% the subject's native space to determine the most likely ROI each vertex
% in a retinotopic pRF map belongs to. Saves the labels with the suffix '_pred'.
%
% Goes through all vertices in the template maps where at least one individual 
% has been labelled & creates a Euclidean searchlight (in sphere space) around it. 
% It uses the mean correlation for both the x0 and y0 maps in the searchlight 
% between the empirical data and the comparison subjects to calculate weighted mean 
% of ROI labels to predict the ROI (rounded to categorise). The final ROI map is 
% further smoothed with the mode function within the searchlight radius.
%
% Optionally, you can adjust the searchlight radius, restrict the eccentricity range, 
% and also smooth the empirical map (only within the search region from the template).
%
%   SrfOrName:  Empirical pRF map (assumes Srf.Data(2:3) are x0 and y0 coordinates)
%   Tmp:        Comparison maps warped from template into native space (using Template2NativeMaps)
%   Kernels:    1 x 2 vector with kernels for analysis:
%                   Radius in mm of search light around each vertex to use (default = 10)
%                   Smoothing kernel (on surface) to apply to empirical map (default = 0, no smoothing)
%   Threshs:    1 x 4 vector with thresholds for analysis:
%                   Threshold of confidence map to be used for labelling ROIs (default = 0.1)
%                   Goodness-of-fit threshold of empirical map (default = 0.01) - applied before smoothing if used
%                   Minimum eccentricity of empirical map to be considered (default = 0)
%                   Maximum eccentricity of empirical map to be considered (default = Inf)
%
% Also returns a Srf with the predicted label for each vertex in the "ROI" value & 
% the proportion of consistent labels across the training data in "Confidence".
%
% 04/05/2026 - Written (DSS)
% 17/06/2026 - Added option to pre-smooth empirical maps (DSS)
% 18/06/2026 - Can now take Srf file as input (DSS)
% 22/06/2026 - Changed input arguments (DSS)
%              Added confidence thresholding (DSS)
% 25/06/2026 - ROI map is now smoothed with mode function (DSS)
% 29/07/2026 - Added R^2 threshold for empirical map (DSS)
% 26/08/2026 - Fixed mistake in help section (DSS)
% 27/08/2026 - Still work in progress... (DSS)
%

% Empirical map
if ischar(SrfOrName)
    % Load Srf
    load(SrfOrName);
else
    % Directly supplied
    Srf = SrfOrName;
end
% Expand if needed
Srf = samsrf_expand_srf(Srf);
Tmp = samsrf_expand_srf(Tmp);

if nargin < 3
    Kernels = [10 0]; % Searchlight radius & smoothing kernel
end
if isscalar(Kernels)
    Kernels = [Kernels 0];
end
if nargin < 4
    Threshs = [0.1 0.01 0 Inf]; % Confidence & R^2 threshold & eccentricity range
end
if isscalar(Threshs)
    Threshs = [Threshs 0.01 0 Inf];
elseif length(Threshs) == 2
    Threshs = [Threshs 0 Inf];
elseif length(Threshs) == 3
    Threshs = [Threshs Inf];
end
% Assign to named variables
Rad = Kernels(1); % Searchlight radius
Smo = Kernels(2); % Smoothing kernel
Con = Threshs(1); % Confidence threshold
GoF = Threshs(2); % Goodness-of-fit threshold
Ecc = Threshs(3:4); % Eccentricity range

%% Prepare maps
mver = find(Tmp.Data(6,:) > 0); % Only vertices labelled in any individual map
% Update version number in Srf
Srf.Version = samsrf_version;
% Filter out bad vertices
Srf.Data(:, Srf.Data(1,:) < GoF) = 0; 
% Smooth map
Srf = samsrf_smooth_sphere(Srf, Smo, mver);

%% Mask out eccentricity range
Eccen = sqrt(sum(Srf.Data(2:3,:).^2)); % Empirical eccentricity map
Emsk = Eccen < Ecc(1) | Eccen > Ecc(2); % Vertices out of eccentricity range
Srf.Data(:,Emsk) = 0;
Tmp.Data(:,Emsk) = 0;
mver = find(sum(Srf.Data) ~= 0); % Create new mask
% % Create meridian maps
% Mer = cart2pol(abs(Srf.Data(2,:)), abs(Srf.Data(3,:)));
% TmpMer = cart2pol(abs(Tmp.Raw_Data(2,:,:)), abs(Tmp.Raw_Data(3,:,:)));

%% Loop thru mask vertices
samsrf_disp('Labelling ROIs...'); 
samsrf_disp(' Please stand by...');
LabMap = zeros(1,length(mver)); % Labelling empirical map
ConMap = zeros(1,length(mver)); % Confidence of label 
parfor i = 1:length(mver)
    % Current vertex info
    v = mver(i);
    L = squeeze(Tmp.Raw_Data(5,v,:)); % ROI labels of this vertex   
    sl = find(sqrt(sum((Srf.Sphere-Srf.Sphere(v,:)).^2,2)) < Rad); % Searchlight around this vertex   
    
    % Closest matching individual map
    Rx = corr(Srf.Data(2,sl)', squeeze(Tmp.Raw_Data(2,sl,:)))'; % Correlations with X-coordinates
    Ry = corr(Srf.Data(3,sl)', squeeze(Tmp.Raw_Data(3,sl,:)))'; % Correlations with Y-coordinates
    W = (Rx+Ry)/2; % Mean correlation to be used as weight
    % W = corr(Mer(sl)', squeeze(TmpMer(1,sl,:)))'; % Correlations with meridian map
    % W(W < 0) = 0; % Negative correlations weighted as zero
    % W(L == 0) = 0; % Any individuals were this vertex is unlabelled don't matter    
    % Pred = wmode(L, W); % Predicted label for this vertex
    Pred = L(W==max(W));
    LabMap(i) = Pred; % Rounded weighted mean of ROI label
    ConMap(i) = mean(L==Pred); % Proportion of consistent labels
end

%% Add prediction to Srf
Srf.Values{end+1} = 'ROI';
Srf.Values{end+1} = 'Confidence';
Srf.Data = [Srf.Data; zeros(2,size(Srf.Vertices,1))];
Srf.Data(end-1,mver) = LabMap; 
Srf.Data(end,mver) = ConMap; 

%% Smooth ROI map?
sLabMap = zeros(1,size(Srf.Vertices,1)); % Smoothed ROI map
samsrf_progbar(0);
for i = 1:length(mver)
    % Current vertex info
    v = mver(i);
    sl = sqrt(sum((Srf.Sphere-Srf.Sphere(v,:)).^2,2)) < Rad; % Searchlight around this vertex   
    sLabMap(v) = mode(Srf.Data(end-1,sl)); % Most frequently occurring ROI label in searchlight
    samsrf_progbar(i/length(mver));
end
Srf.Data(end-1,:) = sLabMap; % Replace with smoothed ROI map

%% Save out ROI labels
if Srf.Hemisphere(1) == 'l'
    h = 'lh_';
elseif Srf.Hemisphere(1) == 'r'
    h = 'rh_';
else
    h = '';
end
% Loop thru ROI labels thresholded by confidence
for r = 1:nanmax(Srf.Data(end-1,:))
    samsrf_srf2label(Srf, [h Tmp.Roi_Names{r} '_pred'], 1, find(Srf.Data(end-1,:)==r & Srf.Data(end,:)>=Con));
end

%% Save maps?
if ischar(SrfOrName)    
    if ~exist('Model','var')
        Model = NaN; % In case model wasn't saved in map file
    end
    save([SrfOrName '_pred'], 'Srf', 'Model');
end


%% Weighted mode of labels
function wL = wmode(L, W)
% Returns mode of labels L weighted by W

uL = unique(L); % Unique labels
uL = uL(uL > 0); % Don't even bother including zero
WbL = []; % Summed weights by unique label
if ~isempty(uL)
    % Loop thru unique labels
    for u = uL'
        WbL = [WbL; u sum(W(L==u))];
    end
    wL = WbL(find(WbL(:,2)==max(WbL(:,2)),1),1); % Return label with highest weight
    % Ensure not empty
    if isempty(wL)
        wL = 0;
    end
else
    wL = 0;
end