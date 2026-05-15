function Srf = samsrf_fieldsign(InSrf, radius, roi, thrsh)
%
% Srf = samsrf_fieldsign(InSrf, [radius=2, roi='', thrsh=0.1])
%
% Calculates a field sign map from the retinotopic maps in InSrf.
%
% The function determines the neighbours (and thus the gradient) for each vertex. 
% All vertices within radius of a given vertex are considered a neighbour. 
%
% If radius is positive, this is using geodesic steps along the surface.
%
% If radius is negative, it is using the Euclidean distance in sphere space.
%   This may be a lot faster although perhaps less accurate.
%
% Only vertices with R^2>=thrsh (default = 0.1) are included. 
%
% This can (and probably should) be limited to a ROI.
%
% Adds the field sign map to Srf.Data as the bottom row. 
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 14/03/2024 - Can now use vector of vertices as ROI (DSS)
% 23/04/2026 - Added option for spherical Euclidean instead of geodesic distance (DSS)
%

%% Default parameters
if nargin < 2
    radius = 2;
    roi = '';
    thrsh = 0.1;
elseif nargin < 3
    roi = '';
    thrsh = 0.1;
elseif nargin < 4
    thrsh = 0.1;
end

% Expand Srf if necessary
InSrf = samsrf_expand_srf(InSrf);

% Load data
Srf = InSrf;
D = Srf.Data;
nver = size(D,2);

% Calculate retinotopy
E = sqrt(D(2,:).^2 + D(3,:).^2);
P = atan2(D(3,:), D(2,:));

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
    Srf.Values{end+1} = 'Field Sign';
    
    % Load sphere surface
    sphV = Srf.Sphere;
    
    % Calculate field sign
    samsrf_disp('  Calculating field signs...');
    for j = 1:si
        if isempty(roi)
            if j == si
                Vs = ((j-1)*50000+1:nver)';
            else
                Vs = ((j-1)*50000+1:j*50000)';
            end
            samsrf_disp(['   Calculating field signs... (Block #' num2str(j) ')']); 
        end
        fsD = zeros(1,length(Vs));
        samsrf_disp('    Please stand by...');
        parfor vi = 1:length(Vs)
            v = Vs(vi);
            % Vertices in ROI within radius
            if radius < 0
                % Euclidean distance in sphere space
                Nv = find(sqrt(sum((sphV-sphV(v,:)).^2,2)) < -radius);
            else
                % Geodesic distance in steps
                Nv = samsrf_georoi(v, radius, Srf.Vertices, Srf.Faces);
            end
            Nv = Nv(Srf.Data(1,Nv) >= thrsh);
            cE = E(Nv);
            cP = P(Nv);
            cX = sphV(Nv,1)';
            cZ = sphV(Nv,3)';

            % Only if there are viable neighbours 
            if length(Nv) > 1 && ~isnan(mean(cE))
                % Determine gradients
                Efx = polyfit(cX, cE, 1); % Slope of eccentricity in left-right axis
                Efz = polyfit(cZ, cE, 1); % Slope of eccentricity in inferior-superior axis
                Pfx = polyfit(cX, cP, 1); % Slope of polar angle in left-right axis
                Pfz = polyfit(cZ, cP, 1); % Slope of polar angle in inferior-superior axis
                E_vec = [Efx(1) Efz(1) 0]; % Eccentricity vector
                P_vec = [Pfx(1) Pfz(1) 0]; % Polar angle vector

                % Field sign
                fsD(1,vi) = sign(sum(cross(E_vec, P_vec)));
            end
        end
    end
    % Save field sign data
    Srf.Data(end,Vs) = fsD;
else
    % No sphere data in Srf
    samsrf_disp('WARNING: Skipping field sign calculation: no sphere data in Srf'); 
end
    