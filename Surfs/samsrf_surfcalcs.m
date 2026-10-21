function Srf = samsrf_surfcalcs(InSrf, Roi, R2Thrsh, Eccens, Kernels, Method, GenCmf)
%
% Srf = samsrf_surfcalcs(InSrf, [Roi='', R2Thrsh=0.01, Eccens=[0 Inf], Kernels=[8 -5 -1], Method='S', GenCmf=[false]])
% 
% Runs a range of surface calculations on a pRF map:
% First, it smoothes the map with a large kernel and calculates the field sign.
% Next, it smoothes the map again with a smaller kernel and calculates the CMF.
% Finally, it rectifies the now smoothed field sign map and returns the Srf.
%
%   InSrf:      Input Srf
%   Roi:        ROI label to restrict analysis
%   R2Thrsh:    R^2 threshold
%   Eccens:     Eccentricity range (1x2 vector for minimum and maximum)
%   Kernels:    Smoothing kernels / search light size (vector: 1: smoothing FWHM, 2:3: field sign radius)
%   Method:     Smoothing method (only first letter needed) from
%                'Sphere':      Distance determined on the spherical surface (SamSrf standard)
%                'Geodeosic':   Distance determined by geodesic steps (FreeSurfer standard)
%                'Dijkstra':    Distance determined by Dijkstra's geodesic distance (best method, but slow)
%   GenCmf:     Toggles whether cortical magnification factor is being calculated 
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 03/09/2026 - Changed default kernels (now use Euclidean distance & Poor Man's Field Sign)
%              Removed second smoothing pass (DSS)
%              Made CMF calculation optional (DSS)
%              Fixed incorrect help section (DSS)
% 04/09/2026 - Changed name of Kernels & fixed bug with defaults (DSS)
% 07/09/2026 - Second pass of eccentricity filtering after calculations (DSS)
%              Changed order of input arguments (DSS)
% 08/09/2026 - Now uses Sereno field sign & wider default kernel (DSS)
% 28/09/2026 - Fixed error in help section (DSS)
%

if nargin < 2
    Roi = '';
end
if nargin < 3
    R2Thrsh = 0.01;
end
if nargin < 4
    Eccens = [0 Inf];
end
if nargin < 5
    Kernels = [8 -5 -1];
end
if isscalar(Kernels) 
    Kernels = [Kernels -5 -1];
elseif length(Kernels) < 3
    Kernels = [Kernels -1];
end
if nargin < 6
    Method = 'S';
end
if nargin < 7
    GenCmf = false;
end

% Expand if needed & store raw data
[Srf, RoiVx] = samsrf_expand_srf(InSrf);
Srf = samsrf_smooth_sphere(Srf, 0); % Remove previous smoothing

% Filter eccentricity 
Ecc = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); % Eccentricities
Srf.Data(:,Ecc < Eccens(1) | Ecc > Eccens(2)) = 0; % Remove anything outside eccentricity range

% Smooth map 
switch upper(Method(1))
    case 'S'
        Srf = samsrf_smooth_sphere(Srf, Kernels(1), Roi, R2Thrsh);
    case 'G'
        Srf = samsrf_smooth_georoi(Srf, Kernels(1), Roi, R2Thrsh);
    case 'D'
        Srf = samsrf_smooth_dijkstra(Srf, Kernels(1), Roi, R2Thrsh);
end

% Calculate field sign
Srf = samsrf_fieldsign(Srf, Kernels(2:3), Roi, R2Thrsh);

% Calculate CMF?
if GenCmf
    Srf = samsrf_cortmagn(Srf, Roi); % Calculate CMF
end

% Filter eccentricity again
Ecc = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); % Eccentricities
Srf.Data(:,Ecc < Eccens(1) | Ecc > Eccens(2)) = 0; % Remove anything outside eccentricity range

% Compress data if needed
Srf = samsrf_compress_srf(Srf, RoiVx);
