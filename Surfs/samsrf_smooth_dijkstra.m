function Srf = samsrf_smooth_dijkstra(InSrf, fwhm, roi, thrsh)
%
% Srf = samsrf_smooth_dijkstra(InSrf, fwhm, [roi='', thrsh=0.01])
%
% Surface-based smoothing of the data in surface data InSrf using Dijkstr's
% geodesic distance algorithm (Dijkstra 1959). The kernel is the full width 
% half maximum (fwhm) of the Gaussian smoothing kernel. For the sake of speed, 
% only vertices within a radius of 4*stdev are used for this smoothing, so it 
% is not completely invertible. It is still very slow though. It is also
% important to note that a smaller smoothing kernel must be used compared
% to the spherical smoothing approach because the geodesic distances are
% much closer. The R^2 threshold of vertices to include is defined by thrsh.
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
% IMPORTANT: Requires the external dijkstra.m function from MatLab Central:
%            https://au.mathworks.com/matlabcentral/fileexchange/20025-dijkstra-s-minimum-cost-path-algorithm
%
% 18/07/2020 - SamSrf 7 version (DSS)
%

if ~exist('dijkstra.m', 'file')
    error('Joseph Kirk''s Dijkstra algorithm is not on your path! (type ''help samsrf_smooth_dijkstra'')');
end

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
Srf.Values = Srf.Values(1:size(Data,1));

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

if isfield(Srf, 'Sphere')
    % Load sphere surface
    sphV = Srf.Sphere;
    
    % Smoothing
    radius = stdev * 4;
    disp('Smoothing vertices...');
    for j = 1:si
        if isempty(roi)
            if j == si
                Vs = ((j-1)*50000+1:nver)';
            else
                Vs = ((j-1)*50000+1:j*50000)';
            end
            disp([' Smoothing vertices... (Block #' num2str(j) ')']); 
        end
        % Smooth data for each mask vertex
        SmoothedData = zeros(size(Data,1), length(Vs));
        parfor vi = 1:length(Vs) 
            % Current vertex
            v = Vs(vi);
            
            % Vertices in a spherical ROI within large radius 
            Dm = sphV - repmat(sphV(v,:), size(sphV,1), 1); % Euclidian distance vectors 
            Nd = sqrt(Dm(:,1).^2 + Dm(:,2).^2 + Dm(:,3).^2); % Euclidian distance in spherical mesh   
            Nv = aVs(Nd <= radius); % Vertices within a sphere on the spherical mesh
            % Calculate geodesic distance matrix
            Vxy = Srf.Vertices(Nv,:);
            Faces = Srf.Faces(sum(ismember(Srf.Faces,Nv),2) == 3,:); % Faces containing vertices 
            nF = NaN(size(Faces)); 
            for fv=1:length(Nv) 
                nF(Faces == Nv(fv)) = fv; % Change face indeces to reflect selected vertices only
            end
            nA = nF(:,[2 3 1]); % Reordered to compute edges
            Dd = dijkstra(Vxy, [nF(:) nA(:)], 1); % Compute Dijkstra's distance

            % Only if in a good neighbourhood
            if mean(Rsq(Nv) > thrsh) >= 0.5            
                % Remove rubbish vertices
                Dd = Dd(Rsq(Nv) > thrsh);
                Nv = Nv(Rsq(Nv) > thrsh);

                % Weight vertices by distance
                W = exp(-(Dd.^2)/(2*stdev.^2));  
                SmoothedData(:,vi) = sum(repmat(W,size(Data,1),1) .* Data(:,Nv),2) / sum(W);
            end
        end
        % Store smoothed data
        Srf.Data(:,Vs) = SmoothedData;
    end

    % Store unsmoothed raw data
    Srf.Raw_Data = Data;
    if iscellstr(Srf.Functional)
        for iStr = 1:length(Srf.Functional)
            Srf.Functional{iStr} = [Srf.Functional{iStr} ' (Smoothed with Dijkstra FWHM=' num2str(fwhm) roistr];
        end
    else
        Srf.Functional = [Srf.Functional ' (Smoothed with Dijkstra FWHM=' num2str(fwhm) roistr];
    end
    disp(['Smoothing finished after ' num2str(toc(t0)/60) ' minutes.']);
    new_line;
else
    warning('Skipping smoothing: no sphere data in Srf');
end