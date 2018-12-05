function Srf = samsrf_smooth_sphere(InSrf, fwhm, roi, thrsh)
%
% Srf = samsrf_smooth_sphere(InSrf, fwhm, [roi='', thrsh=0.01])
%
% Surface-based smoothing of the data in surface data InSrf using the 
% sphere surface in surfdir. So be aware that the smoothing kernel becomes 
% increasingly less geodesic the further away we are from the centre and 
% there may be wrapping artifacts if you set the kernel size really large).
% The kernel is the full width half maximum (fwhm) of the spherical Gaussian 
% smoothing kernel used to determine Euclidian distances of each vertex from 
% its neighbours. For the sake of speed, only vertices within a radius of 
% 3*stdev are used for this smoothing, so it is not completely invertible.
% The R^2 threshold of vertices to include is defined by thrsh. 
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
% 09/08/2018 - SamSrf 6 version (DSS)
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
    D = Srf.Raw_Data;
    if isfield(Srf, 'Values')
        Srf.Values = Srf.Values(1:size(Srf.Raw_Data,1));
    end
else
    D = Srf.Data;
end
sD = zeros(size(Srf.Data));
nver = size(D,2);
aVs = 1:nver;

% Remove smoothing string if it exists
if iscellstr(Srf.Functional)
    for iStr = 1:length(Srf.Functional)
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
        Rsq = Srf.Data(1,:);
    else
        Rsq = ones(1,size(Srf.Data,2));
    end
else
    Rsq = ones(1,size(Srf.Data,2));
end

% Load region of interest
if ~isempty(roi)
    Vs = samsrf_loadlabel(roi);
    si = 1;
    roistr = [' in RoI: ' roi ')'];
    % Set R^2 of vertices outside ROI to 0
    RoiVx = zeros(1, size(Srf.Data,2));
    RoiVx(Vs) = 1;
    Rsq = Rsq .* RoiVx;
else
    si = ceil(nver/50000);  % Smoothing iterations so we don't run out of memory
    roistr = ')';
end

if isfield(Srf, 'Sphere')
    % Load sphere surface
    sphV = Srf.Sphere;

    % Smoothing
    radius = stdev * 3;
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
            % Vertices in a spherical ROI within radius 
            Dm = sphV - repmat(sphV(v,:), size(sphV,1), 1);
            Nd = sqrt(Dm(:,1).^2 + Dm(:,2).^2 + Dm(:,3).^2);
            Nv = aVs(Nd <= radius);
            Nd = Nd(Nd <= radius);

            % Only if in a good neighbourhood
            if mean(Rsq(Nv) > thrsh) >= 0.5            
                % Remove rubbish vertices
                Nd = Nd(Rsq(Nv) > thrsh);
                Nv = Nv(Rsq(Nv) > thrsh);

                % Weight vertices by distance
                W = exp(-(Nd.^2)/(2*stdev.^2))';  
                sD(:,v) = sum(repmat(W,size(D,1),1) .* D(:,Nv),2) / sum(W);
            end
            if wb waitbar(i/length(Vs), h); end
        end
    end
    if wb close(h); end

    % Store smoothed data
    Srf.Data = sD;
    Srf.Raw_Data = D;
    if iscellstr(Srf.Functional)
        for iStr = 1:length(Srf.Functional)
            Srf.Functional{iStr} = [Srf.Functional{iStr} ' (Smoothed with spherical FWHM=' num2str(fwhm) roistr];
        end
    else
        Srf.Functional = [Srf.Functional ' (Smoothed with spherical FWHM=' num2str(fwhm) roistr];
    end
    disp('  Smoothing complete.'); 
else
    warning('Skipping smoothing: no sphere data in Srf');
end    
