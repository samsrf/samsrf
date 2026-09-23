function Srf = samsrf_fieldsign(InSrf, radius, roi, thrsh)
%
% Srf = samsrf_fieldsign(InSrf, [radius=[-5 0], roi='', thrsh=0.1])
%
% Calculates a field sign map from the retinotopic maps in InSrf. 
% You can either calculate the traditional field sign ala Marty Sereno,
%  the sign of the summed cross product of the polar & eccentricity gradients,
%  or a new algorithm we call "Poor Man's Field Sign". 
% 
% Don't be fooled by the name! As far as we can tell this can actually produce
% rather precise delineations. It is simply the normalised gradient (correlation) 
% between sphere Z-coordinates & the polar angle map. It is -not- rectified
% but instead ranges between -1 and 1.
%
% The function determines the neighbours (and thus the gradient) for each vertex. 
% All vertices within radius of a given vertex are considered a neighbour. 
% This is then used to calculate the gradient & thus the local field sign.
%
% If radius(1) is positive, this is using geodesic steps along the surface.
%
% If radius(1) is negative, it is using the Euclidean distance in sphere space.
%   This may be a lot faster although perhaps less accurate.
%
% If radius(2) is zero or undefined the Poor Man's Field Sign algorithm is used.
% If radius(2) is negative, Sereno gradient angle is calculated. This is the 
%  Sereno field sign calculation but not rectified, so it is the angular difference 
%  between the polar angle & eccentricity gradients.           
%
% Only vertices with R^2>=thrsh (default = 0.1) are included. 
%
% This function can (and probably should) be limited to a ROI to reduce computation time.
%  The roi input can be the filename of a label (without extension) or a vector of vertex numbers.
%
% Adds the field sign map to Srf.Data as the bottom row.
% The corresponding Srf.Values entry will either be 'Field Sign', 'Gradient Angle', or 'PMFS'
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 14/03/2024 - Can now use vector of vertices as ROI (DSS)
% 23/04/2026 - Added option for spherical Euclidean instead of geodesic distance (DSS)
% 01/09/2026 - Changed algorithm to increase speed & use correlation (DSS)
%              Fixed issue with wrapping artifacts in right hemisphere (DSS)
% 02/09/2026 - New option for Poor Man's Field Sign algorithm (DSS)
%              Changed default radius to -5 using Euclidean radius (DSS)
% 03/09/2026 - PMFS algorithm now uses polar angles instead of meridians (DSS)
% 04/09/2026 - Now uses samsrf_sphere_roi for Euclidean search lights (DSS)
%              Can also calculate gradient angle from field sign now (DSS)
% 10/09/2026 - Now reports more information (DSS)  
% 11/09/2026 - New default is gradient angle as needed for fast delineation (DSS)
%

%% Default parameters
if nargin < 2
    radius = -5;
    roi = '';
    thrsh = 0.1;
elseif nargin < 3
    roi = '';
    thrsh = 0.1;
elseif nargin < 4
    thrsh = 0.1;
end
% Default to using Gradient Angle
if isscalar(radius)
    radius = [radius -1];
end

% Expand Srf if necessary
InSrf = samsrf_expand_srf(InSrf);

% Load data
Srf = InSrf;
D = Srf.Data;
nver = size(D,2);

% Calculate retinotopy
E = sqrt(D(2,:).^2 + D(3,:).^2); % Eccentricity map
P = atan2(D(3,:), abs(D(2,:))); % Polar angle map (collapsed to right hemifield to avoid wrapping artifacts)

% Load region of interest
if ~isempty(roi)
    if ischar(roi)
        % Load ROI label
        Vs = samsrf_loadlabel(roi);
    else
        % Use list of vertices
        Vs = roi;
    end
    si = 1;
else
    si = ceil(nver/50000);  % Smoothing iterations so we don't run out of memory
end

if isfield(Srf, 'Sphere')
    % Add new data row
    Srf.Data(end+1,:) = 0;
    if radius(2)
        % Traditional Sereno Field Sign
        if radius(2) > 0
            Srf.Values{end+1} = 'Field Sign';
        else
            Srf.Values{end+1} = 'Gradient Angle';
        end
    else
        % Poor Man's Field Sign
        Srf.Values{end+1} = 'PMFS';
    end
    
    % Load sphere surface
    sphV = Srf.Sphere;
    
    % Calculate field sign
    samsrf_disp(['Calculating ' Srf.Values{end} '...']);
    for j = 1:si
        if isempty(roi)
            if j == si
                Vs = ((j-1)*50000+1:nver)';
            else
                Vs = ((j-1)*50000+1:j*50000)';
            end
            samsrf_disp([' Calculating field signs... (Block #' num2str(j) ')']); 
        end
        fsD = zeros(1,length(Vs));
        fsE = zeros(1,length(Vs));
        samsrf_disp(' Please stand by...');
        parfor vi = 1:length(Vs)
            v = Vs(vi);
            % Vertices in ROI within radius
            if radius(1) < 0
                % Euclidean distance in sphere space
                Nv = samsrf_sphere_roi(v, -radius(1), sphV);  
            else
                % Geodesic distance in steps
                Nv = samsrf_georoi(v, radius(1), Srf.Vertices, Srf.Faces);
            end
            Nv = Nv(Srf.Data(1,Nv) >= thrsh);
            if length(Nv) > 3 % No point doing this for a handful of vertices 
                if radius(2)
                    % Traditional Sereno Field Sign
                    R = corr([sphV(Nv,:) P(Nv)' E(Nv)']); % Correlate polar & eccentricity with patch coordinates
                    grdang = -sum(cross(R(5,1:3), R(4,1:3))); % Gradient angle
                    if radius(2) < 0
                        fsD(1,vi) = grdang; % Angle between polar & eccentricity gradients
                    else
                        fsD(1,vi) = sign(grdang); % Field sign based on polar & eccentricity gradients
                    end
                else
                    % Poor Man's Field Sign
                    Slt = [sphV(Nv,3) P(Nv)' E(Nv)']; % Sphere Z-coordinates & polar map in searchlight
                    Slt = Slt(sum(isnan(Slt),2)==0, :); % Remove NaNs because they screw up the edges of the ROI
                    fsD(1,vi) = corr(Slt(:,1), Slt(:,2)); % Normalised gradient of polar map vs Z-coordinates
                end
            else
                fsD(1,vi) = 0;
            end
        end
    end

    % Save field sign data
    Srf.Data(end,Vs) = fsD;
    % Ensure sign is correct
    if radius(2) % Only need to do this in classical field sign
        Srf.Data(end,:) = Srf.Data(end,:) .* sign(Srf.Data(2,:));
    end
else
    % No sphere data in Srf
    samsrf_disp('WARNING: Skipping field sign calculation: no sphere data in Srf'); 
end
samsrf_newline;
    