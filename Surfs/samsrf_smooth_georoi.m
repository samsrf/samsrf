function Srf = samsrf_smooth_georoi(InSrf, fwhm, roi, thrsh)
%
% Srf = samsrf_smooth_georoi(InSrf, fwhm, [roi='', thrsh=0.01])
%
% Surface-based smoothing of the data in surface data InSrf using geodesic 
% steps along the white matter surface. The R^2 threshold of vertices to include is defined by thrsh.
%
% This can and probably should be limited to a ROI.
%
% If InSrf is a data set that has already been smoothed, the function
% automatically selects the unsmoothed data to be smoothed again. This way
% you won't end up smoothing smoothed data repeatedly. By specifying a
% smoothing kernel of 0 you can restore the unsmoothed data set.
%
% Smoothing is only performed if at least half the vertices in the search
% radius are above threshold. This avoids smoothing in meaningless data from
% small clusters containing mostly rubbish.
%
% Stores the smoothed data in Srf.Data. The original raw data are stored 
% inside Srf.Raw_Data.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 26/09/2022 - Fixed bug with NaNs in coordinates (DSS) 
%

%% Default parameters
if nargin == 2
    roi = '';
end
if nargin <= 3
    thrsh = 0.01;
end
t0 = tic;

% Expand Srf if necessary
Srf = samsrf_expand_srf(InSrf);
clear InSrf

% Load data
if isfield(Srf, 'Raw_Data')
    Data = Srf.Raw_Data;
else
    Data = Srf.Data;
end
nver = size(Data,2);
aVs = 1:nver; % All vertex indices
% Remove smoothed value labels
if size(Data,1) < length(Srf.Values) % Only if more values than data rows!
    Srf.Values = Srf.Values(1:size(Data,1));
end

% Remove smoothing string if it exists
if iscellstr(Srf.Functional)
    for iStr=1:length(Srf.Functional)
        ss = strfind(Srf.Functional{iStr}, ' (Smoothed with');
        if isempty(ss)
            ss = length(Srf.Functional{iStr})+1;
        end
        Srf.Functional{iStr} = Srf.Functional{iStr}(1:ss-1);      
    end
else
    ss = strfind(Srf.Functional, ' (Smoothed with');
    if isempty(ss)
        ss = length(Srf.Functional)+1;
    end
    Srf.Functional = Srf.Functional(1:ss-1);
end

% No smoothing?
if fwhm == 0
    if isfield(Srf, 'Raw_Data')
        Srf.Data = Srf.Raw_Data;
        Srf = rmfield(Srf, 'Raw_Data');
    end
    return
else
    % Convert FWHM into standard deviation
    stdev = fwhm / (2*sqrt(2*log(2)));
    % Clear Srf.Data field
    Srf.Data = zeros(size(Data));
end

% Is R^2 present in data?
if isfield(Srf, 'Values')
    if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
        Rsq = Data(1,:);
    else
        Rsq = ones(1,size(Data,2));
    end
else
    Rsq = ones(1,size(Data,2));
end

% Load region of interest
if ~isempty(roi)
    Vs = samsrf_loadlabel(roi);
    si = 1;
    roistr = [' in RoI: ' roi ')'];
    % Set R^2 of vertices outside ROI to 0
    RoiVx = zeros(1, size(Data,2));
    RoiVx(Vs) = 1;
    Rsq = Rsq .* RoiVx;
else
    si = ceil(nver/50000);  % Smoothing iterations so we don't run out of memory
    roistr = ')';
end

samsrf_disp('Smoothing vertices...');
for j = 1:si
    if isempty(roi)
        if j == si
            Vs = ((j-1)*50000+1:nver)';
        else
            Vs = ((j-1)*50000+1:j*50000)';
        end
        samsrf_disp([' Smoothing vertices... (Block #' num2str(j) ')']); 
    end
    samsrf_disp(' Please stand by...');
    % Smooth data for each mask vertex
    SmoothedData = zeros(size(Data,1), length(Vs));
    parfor vi = 1:length(Vs)
        % Current vertex
        v = Vs(vi);
        
        % Deteremine geodesic ROI
        [Nv Dd] = samsrf_georoi(v, fwhm, Srf.Vertices, Srf.Faces);
        
        % Only if in a good neighbourhood
        if mean(Rsq(Nv) > thrsh) >= 0.5            
            % Remove rubbish vertices
            Dd = Dd(Rsq(Nv) > thrsh)';
            Nv = Nv(Rsq(Nv) > thrsh)';

            % Weight vertices by distance
            W = exp(-(Dd.^2)/(2*stdev.^2));  
            SmoothedData(:,vi) = nansum(repmat(W,size(Data,1),1) .* Data(:,Nv),2) / nansum(W);
        end
    end
    % Store smoothed data
    Srf.Data(:,Vs) = SmoothedData;
end

% Store unsmoothed raw data
Srf.Raw_Data = Data;
if iscellstr(Srf.Functional)
    for iStr = 1:length(Srf.Functional)
        Srf.Functional{iStr} = [Srf.Functional{iStr} ' (Smoothed with geodesic FWHM=' num2str(fwhm) roistr];
    end
else
    Srf.Functional = [Srf.Functional ' (Smoothed with geodesic FWHM=' num2str(fwhm) roistr];
end
samsrf_disp(['Smoothing finished after ' num2str(toc(t0)/60) ' minutes.']);
samsrf_newline;
