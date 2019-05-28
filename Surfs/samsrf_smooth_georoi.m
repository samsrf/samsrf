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
% 27/05/2019 - Fixed bug with thresholding when using already smoothed data (DSS)
%

%% Default parameters
if nargin == 2
    roi = '';
end
if nargin <= 3
    thrsh = 0.01;
end
% Check if waitbar is to be used
wb = samsrf_waitbarstatus;

% Expand Srf if necessary
InSrf = samsrf_expand_srf(InSrf);

% Load data
Srf = InSrf;
if isfield(Srf, 'Raw_Data')
    Data = Srf.Raw_Data;
else
    Data = Srf.Data;
end
SmoothedData = zeros(size(Data));
nver = size(Data,2);
aVs = 1:nver; % All vertex indices

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
end

% Is R^2 present in data?
if isfield(Srf, 'Values')
    if strcmpi(Srf.Values{1}, 'R^2')
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

if wb h = waitbar(0, 'Smoothing vertices...', 'Units', 'pixels', 'Position', [100 100 360 70]); end
disp('  Smoothing vertices...');
for j = 1:si
    if wb waitbar(0, h); end
    i = 0;
    if isempty(roi)
        if j == si
            Vs = ((j-1)*50000+1:nver)';
        else
            Vs = ((j-1)*50000+1:j*50000)';
        end
        if wb waitbar(0, h, ['Smoothing vertices... (Block #' num2str(j) ')']); end
    end
    for v = Vs'
        i = i + 1;        
        % Deteremine geodesic ROI
        [Nv Dd] = samsrf_georoi(v, fwhm, Srf.Vertices, Srf.Faces);
        
        % Only if in a good neighbourhood
        if mean(Rsq(Nv) > thrsh) >= 0.5            
            % Remove rubbish vertices
            Dd = Dd(Rsq(Nv) > thrsh)';
            Nv = Nv(Rsq(Nv) > thrsh)';

            % Weight vertices by distance
            W = exp(-(Dd.^2)/(2*stdev.^2));  
            SmoothedData(:,v) = sum(repmat(W,size(Data,1),1) .* Data(:,Nv),2) / sum(W);
        end
        if wb waitbar(i/length(Vs), h); end
    end
end
if wb close(h); end;

% Store smoothed data
Srf.Data = SmoothedData;
Srf.Raw_Data = Data;
if iscellstr(Srf.Functional)
    for iStr = 1:length(Srf.Functional)
        Srf.Functional{iStr} = [Srf.Functional{iStr} ' (Smoothed with geodesic FWHM=' num2str(fwhm) roistr];
    end
else
    Srf.Functional = [Srf.Functional ' (Smoothed with geodesic FWHM=' num2str(fwhm) roistr];
end
disp('  Smoothing complete.');
